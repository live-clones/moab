#include "moab/HiReconstruction.hpp"
#include "moab/Solvers.hpp"

#include <math.h>
#include <deque>
#include <iostream>

namespace moab
{
	HiReconstruction::HiReconstruction(Core *impl, ParallelComm *comm, EntityHandle meshIn, int minpnts, bool recwhole)
	 : mbImpl(impl),pcomm(comm),_mesh2rec(meshIn),_MINPNTS(minpnts){
	 	assert(NULL!=impl);
	 	ErrorCode error;

	 #ifdef MOAB_HAVE_MPI
	 	if(!pcomm){
	 		pcomm = moab::ParallelComm::get_pcomm(mbImpl,0);
	 	}
	 #endif

	 	error = initialize(recwhole);
	 	if(MB_SUCCESS!=error){
	 		std::cout << "Error initializing HiReconstruction\n" << std::endl;
	 		exit(1);
	 	}
	 }

	HiReconstruction::~HiReconstruction(){
	#ifdef MOAB_HAVE_AHF
	 	ahf = NULL;
	#else
	 	delete ahf;
	#endif
	}

	ErrorCode HiReconstruction::initialize(bool recwhole){
		ErrorCode error;

		ahf = new HalfFacetRep(mbImpl,pcomm,_mesh2rec);
		if(!ahf){
			return MB_MEMORY_ALLOCATION_FAILED;
		}

		error = ahf->initialize(); MB_CHK_ERR(error);
		if(HalfFacetRep::CURVE==ahf->thismeshtype){
			_dim = 1; _MAXPNTS = 13;
		}else if(HalfFacetRep::SURFACE==ahf->thismeshtype){
			_dim = 2; _MAXPNTS = 128;
		}else{
			MB_SET_ERR(MB_FAILURE,"Encountered a non-manifold mesh or a mesh with volume elements")
		}

		error = ahf->get_entity_ranges(_inverts,_inedges,_infaces,_incells); MB_CHK_ERR(error);

		//get locally hosted vertices by filtering pstatus
	#ifdef MOAB_HAVE_MPI
		error = pcomm->filter_pstatus(_inverts,PSTATUS_GHOST,PSTATUS_NOT,-1,&_verts2rec); MB_CHK_ERR(error);
	#else
		_verts2rec = _inverts;
	#endif
		_nv2rec = _verts2rec.size();

		if(recwhole){
			//compute normals(surface) or tangent vector(curve) for all locally hosted vertices
			if(2==_dim){
				compute_average_vertex_normals_surf();
			}else if(1==_dim{
				compute_average_vertex_tangents_curve();
			}else{
				MB_SET_ERR(MB_FAILURE,"Unknow space dimension");
			}
			_hasderiv = true;
		}
		return error;
	}

	/***************************************************
	 *  User Interface for Reconstruction of Geometry  *
	 ***************************************************/

	ErrorCode HiReconstruction::reconstruct3D_surf_geom(int degree, bool interp, bool safeguard, bool reset){
		assert(2==_dim);
		if(_hasfittings&&!reset){
			//This object has precomputed fitting results and user don't want to reset
			return MB_SUCCESS;
		}else{
			_initfittings = _hasfittings = false;
		}
		//initialize for geometric information
		initialize_surf_geom(degree);
		ErrorCode error;
		double *coeffs,*coords;
		int *degree_out;
		for(Range::iterator ivert=_verts2rec.begin();ivert!=_verts2rec.end();++ivert){
			int index = _verts2rec.index(*ivert);
			size_t istr = _vertID2coeffID[index];
			coords = &(_local_coords[9*index]);
			coeffs = &(_local_fit_coeffs[istr]);
			degree_out = &(_degrees_out[index]);
			_interps[index] = interp;
			error  = polyfit3d_walf_surf_vertex(*ivert,interp,degree,_MINPNTS,safeguard,coords,degree_out,coeffs);
			MB_CHK_ERR(error);
		}
		_geom = HISURFACE;
		_hasfittings = true;
		return error;
	}

	ErrorCode HiReconstruction::reconstruct3D_surf_geom(size_t npts, int* degrees, bool* interps, bool safeguard, bool reset){

	}

	ErrorCode HiReconstruction::reconstruct3D_curve_geom(int degree, bool interp, bool safeguard, bool reset){

	}

	ErrorCode HiReconstruction::reconstruct3D_curve_geom(size_t npts, int* degrees, bool* interps, bool safeguard, bool reset){

	}

	ErrorCode HiReconstruction::polyfit3d_walf_surf_vertex(const EntityHandle vid, const bool interp, int degree, int minpnts, const bool safeguard, double* coords, int* degree_out, double* coeffs){
		assert(_dim==2);
		ErrorCode error;
		int ring = estimate_num_rings(degree,interp);
		//get n-ring neighbors
		Range ngbvs;
		error = obtain_nring_ngbvs(vid,ring,minpnts,ngbvs);
		//get coordinates;
		int nverts = ngbvs.size(); assert(nverts);
		double *ngbcoords = new double[nverts*3];
		error = mbImpl->get_coords(ngbvs,ngbcoords); MB_CHK_ERR(error);
		//get normals
		double *ngbnrms = new double[nverts*3];
		error = get_normals_surf(ngbvs,ngbnrms); MB_CHK_ERR(error);
		//local WLS fitting
		int ncoeffspv = (degree+2)*(degree+1)/2;
		int degree_pnt,degree_qr;
		polyfit3d_surf_get_coeff(nverts,ngbcoords,ngbnrms,degree,interp,safeguard,9,coords,nceoffspv,coeffs,degree_out,&degree_pnt,&degree_qr);
		delete [] ngbcoords; delete [] ngbnrms;
		return error;
	}

	ErrorCode HiReconstruction::polyfit3d_walf_curve_vertex(const EntityHandle vid, const bool interp, int degree, int minpnts, const bool safeguard, double* coords, int* degree_out, double* coeffs){

	}

	/**************************************************************
	 *  User Interface for Evaluation via Reconstructed Geometry  *
	 **************************************************************/

	ErrorCode HiReconstruction::hiproj_walf_in_element(EntityHandle elem, const int nvpe, const int npts2fit, const double* naturalcoords2fit, double* newcoords){
		assert(newcoords);
		ErrorCode error;
		if(!_hasfittings){
			MB_SET_ERR(MB_FAILURE,"There is no existing fitting results");
		}
		//check correctness of input
		for(int i=0;i<npts2fit;++i){
			if(!check_barycentric_coords(nvpe,naturalcoords2fit+i*nvpe)){
				MB_SET_ERR(MB_FAILURE,"Wrong barycentric coordinates");
			}
		}

		//get connectivity table
		std::vector<EntityHandle> elemconn;
		error = mbImpl->get_connectivity(&elem,1,elemconn); MB_CHK_ERR(error);
		if(nvpe!=elemconn.size()){
			MB_SET_ERR(MB_FAILURE,"element connectivity table size doesn't match input size");
		}

		double *elemcoords = new double[nvpe*3];
		error = mbImpl->get_coords(elemconn,nvpe,elemcoords); MB_CHK_ERR(error);

		double *coords2fit = new double[3*npts2fit]();
		for(int i=0;i<npts2fit;++i){
			for(int j=0;j<nvpe;++j){
				coords2fit[3*i] += naturalcoords2fit[i*nvpe+j]*elemcoords[3*j];
				coords2fit[3*i+1] += naturalcoords2fit[i*nvpe+j]*elemcoords[3*j+1];
				coords2fit[3*i+2] += naturalcoords2fit[i*nvpe+j]*elemcoords[3*j+2];
			}
		}

		double *hiproj_new = new double[3*npts2fit];
		//initialize output
		for(int i=0;i<npts2fit;++i){
			newcoords[3*i] = newcoords[3*i+1] = newcoords[3*i+2] = 0;
		}
		//for each input vertex, call nvpe fittings and take average
		for(int j=0;j<nvpe;++j){
			error = hiproj_walf_around_vertex(elemconn[j],npts2fit,coords2fit,hiproj_new); MB_CHK_ERR(error);
			for(int i=0;i<npts;++i){
				newcoords[3*i] += naturalcoords2fit[i*nvpe+j]*hiproj_new[3*i];
				newcoords[3*i+1] += naturalcoords2fit[i*nvpe+j]*hiproj_new[3*i+1];
				newcoords[3*i+2] += naturalcoords2fit[i*nvpe+j]*hiproj_new[3*i+2];
			}
		}
		delete [] elemcoords; delete [] coords2fit; delete [] hiproj_new;
	}

	ErrorCode HiReconstruction::hiproj_walf_around_vertex(EntityHandle vid, const int npts2fit, const double* coords2fit, double* hiproj_new){
		if(!_hasfittings){
			MB_SET_ERR(MB_FAILURE,"There is no existing fitting results");
		}
		//get center of local coordinates system
		double local_origin[3];
		error = mbImpl->get_coords(&vid,1,local_origin); MB_CHK_ERR(error);
		//get local fitting parameters
		int index = _verts2rec.index(vid);
		bool interp = _interps[index];
		int local_deg = _degrees_out[index];
		double *uvw_coords,*local_coeffs;
		if(_geom==HISURFACE){
			uvw_coords = &(_local_coords[9*index]);
			int ncoeffs = (local_deg+2)*(local_deg+1)>>1;
			size_t istr = _vertID2coeffID[index];
			local_coeffs = &(_local_fit_coeffs[istr]);
			walf3d_surf_vertex_eval(local_origin,uvw_coords,local_deg,local_coeffs,interp,npts2fit,coords2fit,hiproj_new);
		}else if(_geom==HI3DCURVE){
			uvw_coords = &(_local_coords[3*index]);
			size_t istr = _vertID2coeffID[index];
			local_coeffs = &(_local_fit_coeffs[istr]);
			walf3d_curve_vertex_eval(local_origin,uvw_coords,local_deg,local_coeffs,interp,npts2fit,coords2fit,hiproj_new);
		}
		return error;
	}

	void HiReconstruction::walf3d_surf_vertex_eval(const double* local_origin, const double* local_coords, const int local_deg, const double* local_coeffs, const bool interp, const int npts2fit, const double* coords2fit, double* hiproj_new){

	}

	void HiReconstruction::walf3d_curve_vertex_eval(const double* local_origin, const double* local_coords, const int local_deg, const double* local_coeffs, const bool interp, const int npts2fit, const double* coords2fit, double* hiproj_new){

	}

	/****************************************************************
	 *  Basic Internal Routines to initialize and set fitting data  *
	 ****************************************************************/

	 int HiReconstruction::estimate_num_rings(int degree, bool interp){
	 	return interp?(degree+1)/2:(degree+2)/2;
	 }

	 ErrorCode HiReconstruction::obtain_nring_ngbvs(const EntityHandle vid, int ring, const int minpnts, Range& ngbvs){
	 	ErrorCode error;
	 	std::deque<EntityHandle> todo;
	 	todo.push_back(vid); ngbvs.insert(vid);
	 	EntityHandle pre,nxt;
	 	for(int i=1;i<=ring;++i){
	 		int count = todo.size();
	 		while(count){
	 			EntityHandle center = todo.front();
	 			todo.pop_front(); --count;
	 			std::vector<EntityHandle> adjents;
	 			error = ahf->get_up_adjacencies(center,_dim,adjents); MB_CHK_ERR(error);
	 			for(int j=0;j<adjents.size();++j){
	 				std::vector<EntityHandle> elemconn;
	 				error = mbImpl->get_connectivity(&adjents[j],1,elemconn); MB_CHK_ERR(error);
	 				int nvpe = elemconn.size();
	 				for(int k=0;k<nvpe;++k){
	 					if(elemconn[k]==center){
	 						pre = k==0?elemconn[nvpe-1]:elemconn[k-1];
	 						nxt = elemconn[(k+1)%nvpe];
	 						if(ngbvs.find(pre)==ngbvs.end()){
	 							ngbvs.insert(pre);
	 							todo.push_back(pre);
	 						}
	 						if(ngbvs.find(nxt)==ngbvs.end()){
	 							ngbvs.insert(nxt);
	 							todo.push_back(nxt);
	 						}
	 						break;
	 					}
	 				}
	 			}
	 		}
	 		if(_MAXPNTS<=ngvs.size()){
	 			//obtain enough points
	 			return error;
	 		}
	 		if(!todo.size()){
	 			//current ring cannot introduce any points, return incase deadlock
	 			return error;
	 		}
	 		if(i==ring&&minpnts>ngbvs.size()){
	 			//reach maximum ring but not enough points
	 			++ring;
	 		}
	 	}
	 	return error;
	 }

	 void HiReconstruction::initialize_surf_geom(const int degree){
	 	if(!_hasderiv){
	 		compute_average_vertex_normals_surf();
	 		_hasderiv = true;
	 	}
	 	if(!_initfittings){
	 		int nceoffspv = (degree+2)*(degree+1)/2;
	 		_degrees_out.assign(_nv2rec,0);
	 		_interps.assign(_nv2rec,false);
	 		_vertID2coeffID.reserve(_nv2rec);
	 		_local_fit_coeffs.assign(_nv2rec*nceoffspv,0);
	 		for(size_t i=0;i<_nv2rec;++i){
	 			_vertID2coeffID.push_back(i*ncoeffspv);
	 		}
	 		_initfittings = true;
	 	}
	 }

	 void HiReconstruction::initialize_surf_geom(const size_t npts, const int* degrees){

	 }

	 void HiReconstruction::initialize_3Dcurve_geom(const int degree){

	 }

	 void HiReconstruction::initialize_3Dcurve_geom(const size_t npts, const int* degrees){
	 	
	 }

	 ErrorCode HiReconstruction::set_geom_data_surf(const EntityHandle vid, const double* coords, const double degree_out, const double* coeffs, bool interp){

	 }

	 ErrorCode HiReconstruction::set_geom_data_3Dcurve(const EntityHandle vid, const double* coords, const double degree_out, const double* coeffs, bool interp){

	 }

	 ErrorCode HiReconstruction::get_geom_data_surf(const EntityHandle vid, double* coords, double& degree_out, double* coeffs, bool& interp){

	 }

	 ErrorCode HiReconstruction::get_geom_data_3Dcurve(const EntityHandle vid, double* coords, double& degree_out, double* coeffs, bool& interp){

	 }

	 ErrorCode HiReconstruction::average_vertex_normal(const EntityHandle vid, double* nrm){
	 	ErrorCode error;
	 	std::vector<EntityHandle> adjfaces;
	 	error = ahf->get_up_adjacencies(vid,2,adjfaces); MB_CHK_ERR(error);
	 	int npolys = adjfaces.size();
	 	if(!nploys){
	 		MB_SET_ERR(MB_FAILURE,"Vertex has no incident 2D entities");
	 	}else{
	 		double v1[3],v2[3],v3[3],a[3],b[3],c[3];
	 		for(int i=0;i<npolys;++i){
	 			//get incident "triangles"
	 			std::vector<EntityHandle> elemconn;
	 			error = mbImpl->get_connectivity(&adjfaces[i],1,elemconn); MB_CHK_ERR(error);
	 			EntityHandle pre,nxt;
	 			int nvpe = elemconn.size();
	 			for(int j=0;j<nvpe;++j){
	 				if(vid==elemconn[j]){
	 					pre = j==0?elemconn[nvpe-1]:elemconn[j-1];
	 					nxt = elemconn[(j+1)%nvpe];
	 					break;
	 				}
	 			}
	 			//compute area weighted normals
	 			error = mbImpl->get_coords(&pre,1,a); MB_CHK_ERR(error);
	 			error = mbImpl->get_coords(&vid,1,b); MB_CHK_ERR(error);
	 			error = mbImpl->get_coords(&nxt,1,c); MB_CHK_ERR(error);
	 			vec_linear_operation(3,1,c,-1,b,v1);
	 			vec_linear_operation(3,1,a,-1,b,v2);
	 			vec_crossprod(v1,v2,v3);
	 			vec_linear_operation(3,1,nrm,1,v3,nrm);
	 		}
	 		double len=vec_normalize(3,nrm,nrm); assert(len);
	 	}
	 	return error;
	 }

	 ErrorCode HiReconstruction::compute_average_vertex_normals_surf(){
	 	if(_hasderiv){
	 		return MB_SUCCESS;
	 	}
	 	ErrorCode error;
	 	_local_coords.assign(9*_nv2rec,0);
	 	size_t index=0;
	 	for(Range::iterator ivert=_verts2rec.begin();ivert!=_verts2rec.end();++ivert,++index){
	 		error = average_vertex_normal(*ivert,&(_local_coords[9*index+6])); MB_CHK_ERR(error);
	 	}
	 	return error;
	 }

	 ErrorCode HiReconstruction::get_normals_surf(const Range& vertsh, double* nrms){
	 	ErrorCode error = MB_SUCCESS;
	 	if(_hasderiv){
	 		size_t id=0;
	 		for(Range::iterator ivert=vertsh.begin();ivert!=vertsh.end();++ivert,++id){
	 			int index = _verts2rec.index(*ivert);
	 		#ifdef MOAB_HAVE_MPI
	 			if(-1==index){
	 				//ghost vertex
	 				error = average_vertex_normal(*ivert,nrms+3*id); MB_CHK_ERR(error);
	 			}else{
	 				nrms[3*id] = _local_coords[9*index+6];
	 				nrms[3*id+1] = _local_coords[9*index+7];
	 				nrms[3*id+2] = _local_coords[9*index+8];
	 			}
	 		#else
	 			assert(1!=index);
	 			nrms[3*id] = _local_coords[9*index+6];
	 			nrms[3*id+1] = _local_coords[9*index+7];
	 			nrms[3*id+2] = _local_coords[9*index+8];
	 		#endif
	 		}
	 	}else{
	 		size_t id=0;
	 		for(Range::iterator ivert=vertsh.begin();ivert!=vertsh.end();++ivert,++id){
	 			error = average_vertex_normal(*ivert,nrms+3*id); MB_CHK_ERR(error);
	 		}
	 	}
	 	return error;
	 }

	 ErrorCode HiReconstruction::average_vertex_tangent(const EntityHandle vid, double* tang){
	 	ErrorCode error;
	 	std::vector<EntityHandle> adjedges;
	 	error = ahf->get_up_adjacencies_1d(vid,adjedges); MB_CHK_ERR(error);
	 	int nedges = adjedges.size();
	 	if(!nedges){
	 		MB_SET_ERR(MB_FAILURE,"Vertex has no incident edges");
	 	}else{
	 		assert(nedges<=2);
	 		tang[0] = tang[1] = tang[2] = 0;
	 		for(int i=0;i<nedges;++i){
	 			std::vector<EntityHandle> edgeconn;
	 			error = mbImpl->get_connectivity(&adjedges[i],1,edgeconn);
	 			double istr[3],iend[3],t[3];
	 			error = mbImpl->get_coords(&(edgeconn[0]),1,istr);
	 			error = mbImpl->get_coords(&(edgeconn[1]),1,iend);
	 			vec_linear_operation(3,1,iend,-1,istr,t);
	 			vec_linear_operation(3,1,tang,1,t,tang);
	 		}
	 		double len=vec_normalize(3,tang,tang); assert(len);
	 	}
	 }

	 ErrorCode HiReconstruction::compute_average_vertex_tangents_curve(){
	 	if(_hasderiv){
	 		return MB_SUCCESS;
	 	}
	 	ErrorCode error;
	 	_local_coords.assign(3*_nv2rec,0);
	 	size_t index=0;
	 	for(Range::iterator ivert=_verts2rec.begin();ivert!=_verts2rec.end();++ivert,++index){
	 		error = average_vertex_tangent(*ivert,&(_local_coords[3*index])); MB_CHK_ERR(error);
	 	}
	 	return error;
	 }

	 ErrorCode HiReconstruction::get_tangents_curve(const Range& vertsh, double* tangs){

	 }

	/************************************************
	*	Internal Routines for local WLS fittings	*
	*************************************************/

	 void HiReconstruction::polyfit3d_surf_get_coeff(const int nverts, const double* ngbcors, const double* ngbnrms, int degree, const bool interp, const bool safeguard, const int ncoords, double* coords, const int ncoeffs, double* coeffs, double* degree_out, double* degree_pnt, double* degree_qr){
	 	if(nverts<=0){
	 		return;
	 	}
	 	//step 1. copmute local coordinate system
	 	double nrm[3] = {ngbnrms[0],ngbnrms[1],ngbnrms[2]}, tang1 = {0,0,0}, tang2 = {0,0,0};
	 	if(abs(nrm[0])>abs(nrms[1])&&abs(nrm[0])>abs(nrm[2])){
	 		tang1[1] = 1.0;
	 	}else{
	 		tang1[0] = 1.0;
	 	}

	 	vec_projoff(3,tang1,nrm,tang1);
	 	double len1 = vec_normalize(3,tang1,tang1); assert(len1);
	 	vec_crossprod(nrm,tang1,tang2);
	 	if(9==ncoords&&coords){
	 		coords[0] = tang1[0]; coords[1] = tang1[1]; coords[2] = tang1[2];
	 		coords[3] = tang2[0]; coords[4] = tang2[1]; coords[5] = tang2[2];
	 		coords[6] = nrm[0]; coords[7] = nrm[1]; coords[8] = nrm[2];
	 	}

	 	//step 2. project onto local coordinates system
	 	int npts2fit = nverts-interp;
	 	if(0==npts2fit){
	 		*degree_out = *degree_pnt = *degree_qr = 0;
	 		coeffs[0] = 0;
	 		return;
	 	}
	 	double *us = new double[npts2fit*2];
	 	double *bs = new double[npts2fit];
	 	for(int i=interp;i<nverts;++i){
	 		int k = i-interp;
	 		double uu[3];
	 		vec_linear_operation(3,1,ngbcoords[3*i],-1,ngbcoords[0],uu);
	 		us[k*2] = vec_innerprod(3,tang1,uu); us[k*2+1] = vec_innerprod(3,tang2,uu);
	 		bs[k] = vec_innerprod(3,nrm,uu);
	 	}

	 	//step 3. compute weights
	 	double *ws = new double[npts2fit];
	 	int nzeros = compute_weights(npts2fit,2,us,nverts,ngbnrms,degree,_MINEPS,ws);

	 	//step 4. adjust according to zero-weights
	 	if(nzeros){
	 		if(nzeros==npts2fit){
	 			*degree_out = *degree_pnt = *degree_qr = 0;
	 			coeffs[0] = 0;
	 			return;
	 		}
	 		int index=0;
	 		for(int i=0;i<npts2fit;++i){
	 			if(ws[i]){
	 				if(i>index){
	 					us[index*2] = us[i*2]; us[index*2+1] = us[i*2+1];
	 					bs[index] = bs[i]; ws[index] = ws[i];
	 				}
	 				++index;
	 			}
	 		}
	 		npts2fit -= nzeros; assert(index==npts2fit);
	 		/*us = realloc(us,npts2fit*2*sizeof(double));
	 		bs = realloc(bs,npts2fit*sizeof(double));
	 		ws = realloc(ws,npts2fit*sizeof(double));*/
	 	}

	 	//step 5. fitting
	 	eval_vander_bivar_cmf(npts2fit,us,1,bs,degree,ws,interp,safeguard,degree_out,degree_pnt,degree_qr);

	 	//step 6. organize output
	 	int ncoeffs_out = (*degree_out+2)*(*degree_out+1)/2;
	 	assert(ncoeffs_out<=ncoeffs);
	 	coeffs[0] = 0;
	 	for(int j=0;j<ncoeffs_out-interp;++j){
	 		coeffs[j+interp] =  bs[j];
	 	}
	 	delete [] us; delete [] bs; delete [] ws;
	 }

	 void HiReconstruction::eval_vander_bivar_cmf(const int npts2fit, const double* us, const int ndim, double* bs, int degree, const double* ws, const bool interp, const bool safeguard, int* degree_out, int* degree_pnt, int* degree_qr){
	 	//step 1. adjust the degree according to number of points to fit
	 	int ncols = ((degree+2)*(degree+1))>>1-interp;
	 	while(1<degree&&npts2fit<ncols){
	 		--degree;
	 		ncols = ((degree+2)*(degree+1))>>1-interp;
	 	}
	 	*degree_out = degree;

	 	//step 2. construct Vandermonde matrix, stored in columnwise
	 	double *V_init = new double[npts2fit*(ncols+interp)];
	 	gen_vander_bivar(npts2fit,us,degree,V_init);
	 	//remove the first column of 1s if interpolation
	 	double* V;
	 	if(interp){
	 		V = new double[npts2fit*ncols];
	 		std::memcpy(V,V_init+npts2fit,ncols*npts2fit*sizeof(double));
	 		delete [] V_init; V_init = 0;
	 	}else{
	 		V = V_init;
	 	}

	 	//step 3. Scale rows to assign different weights to different points
	 	for(int i=0;i<npts2fit;++i){
	 		for(int j=o;j<ncols;++j){
	 			V[j*npts2fit+i] *= ws[i];
	 		}
	 		for(int k=0;k<ndim;++k){
	 			bs[k*npts2fit+i] *= ws[i];
	 		}
	 	}

	 	//step 4. scale columns to reduce condition number
	 	double *ts = new double[ncols];
	 	rescale_matrix(npts2fit,ncols,V,ts);

	 	//step 5. Perform Householder QR factorization
	 	double *D = new double[ncols];
	 	int rank;
	 	qr_polyfit_safeguarded(V,npts2fit,ncols,D,rank);

	 	//step 6. adjust degree of fitting according to rank of Vandermonde matrix
	 	int ncols_sub = ncols;
	 	while(rank<ncols_sub){
	 		--degree;
	 		if(degree==0){
	 			//surface is flat, return 0
	 			*degree_out = *degree_qr = degree;
	 			for(int i=0;i<npts2fit;++i){
	 				for(int k=0;k<ndim;++k){
	 					bs[k*npts2fit+i] = 0;
	 				}
	 			}
	 			return;
	 		}else{
	 			ncols_sub = ((degree+2)*(degree+1))>>1-interp;
	 		}
	 	}
	 	*degree_qr = degree;

	 	//step 7. compute Q'b
	 	compute_qtransposeB(npts2fit,ncols_sub,V,ndim,bs);

	 	//step 8. perform backward substitution and scale the solution
	 	for(int i=0;i<ncols_sub;++i){
	 		V[i*npts2fit+i] = D[i];
	 	}

	 	//backsolve
	 	if(safeguard){
	 		backsolve_polyfit_safeguarded();
	 	}else{
	 		backsolve(npt2fit,ncols_sub,V,1,bs,ts);
	 		*degree_out = degree;
	 	}
	 	if(V_init){
	 		delete [] V_init;
	 	}else{
	 		delete [] V;
	 	}
	 }

	 void HiReconstruction::polyfit3d_curve_get_coeff(const int nverts, const double* ngbcors, const double* ngbtangs, int degree, const bool interp, const bool safeguard, const int ncoords, double* coords, const int ncoeffs, double* coeffs, double* degree_out){

	 }

	 void HiReconstruction::eval_vander_univar_cmf(const int npts2fit, const double* us, const int ndim, double* bs, int degree, const double* ws, const bool interp, const bool safeguard, int* degree_out){

	 }

	 int HiReconstruction::compute_weights(const int nrows, const int ncols, double* us, const int nngbs, double* ngbnrms, const int degree, const double toler, double* ws){
	 	assert(nrows<=_MAXPNTS)&&ws;
	 	bool interp=false;
	 	if(nngbs!=nrows){
	 		assert(nngbs=nrows+1);
	 		interp = true;
	 	}
	 	double epsilon = 1e-2;

	 	//First, compute squared distance from each input piont to the center
	 	for(int i=0;i<nrows;++i){
	 		ws[i] = vec_2norm(ncols,us+i*ncols);
	 	}

	 	//Second, compute a small correction termt o guard against zero
	 	double h=0;
	 	for(int i=0;i<nrows;++i){
	 		h += ws[i];
	 	}
	 	h /= (double) nrows;

	 	//Finally, compute the weights for each vertex
	 	int nzeros = 0;
	 	for(int i=0;i<nrows;++i){
	 		double costheta = vec_innerprod(3,ngbnrms,ngbnrms+i+interp);
	 		if(costheta>toler){
	 			ws[i] = costheta*pow(ws[i]/h+epsilon,-1*(double) degree/2.0);
	 		}else{
	 			ws[i] = 0;
	 			++nzeros;
	 		}
	 	}
	 	return nzeros;
	 }
	 bool HiReconstruction::check_barycentric_coords(const int nws, const double* naturalcoords){
		double sum=0;
		for(int i=0;i<nws;++i){
			if(naturalcoords[i]<0){
				return false;
			}
			sum += naturalcoords[i];
		}
		if(abs(1-sum)>_MINEPS){
			return false;
		}else{
			return true;
		}
	}
}//namespace moab