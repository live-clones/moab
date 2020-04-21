///////////////////////////////////////////////////////////////////////////////
///
///	\file    TempestOnlineMap.h
///	\author  Vijay Mahadevan
///	\version November 20, 2017
///

#ifndef _TEMPESTONLINEMAP_H_
#define _TEMPESTONLINEMAP_H_

#include "moab/MOABConfig.h"

// Tempest includes
#ifdef MOAB_HAVE_TEMPESTREMAP
#include "moab/Remapping/TempestRemapper.hpp"
#include "OfflineMap.h"
#endif

#include <string>
#include <vector>

#ifdef MOAB_HAVE_EIGEN
#include <Eigen/Sparse>
#endif

///////////////////////////////////////////////////////////////////////////////

#define RECTANGULAR_TRUNCATION
// #define TRIANGULAR_TRUNCATION

///////////////////////////////////////////////////////////////////////////////

// Forward declarations
class Mesh;

///////////////////////////////////////////////////////////////////////////////

namespace moab
{

///	<summary>
///		An offline map between two Meshes.
///	</summary>
class TempestOnlineMap : public OfflineMap {

public:

	///	<summary>
	///		Generate the metadata associated with the offline map.
	///	</summary>
	TempestOnlineMap(moab::TempestRemapper* remapper);

	///	<summary>
	///		Define a virtual destructor.
	///	</summary>
	virtual ~TempestOnlineMap();

public:

    // Input / Output types
    enum DiscretizationType
    {
        DiscretizationType_FV,
        DiscretizationType_CGLL,
        DiscretizationType_DGLL,
        DiscretizationType_PCLOUD
    };

	///	<summary>
	///		Generate the offline map, given the source and target mesh and discretization details.
	///     This method generates the mapping between the two meshes based on the overlap and stores
	///     the result in the SparseMatrix.
	///	</summary>
	moab::ErrorCode GenerateRemappingWeights( std::string strInputType="fv", std::string strOutputType="fv",
                                        const int nPin=1, const int nPout=1,
                                        bool fBubble=false, int fMonotoneTypeID=0,
                                        bool fVolumetric=false, bool fNoConservation=false, bool fNoCheck=false,
                                        const std::string srcDofTagName="GLOBAL_ID", const std::string tgtDofTagName="GLOBAL_ID",
                                        const std::string strVariables="",
                                        const std::string strInputData="", const std::string strOutputData="",
                                        const std::string strNColName="", const bool fOutputDouble=false,
                                        const std::string strPreserveVariables="", const bool fPreserveAll=false, const double dFillValueOverride=0.0,
                                        const bool fInputConcave = false, const bool fOutputConcave = false );

	///	<summary>
	///		Generate the metadata associated with the offline map.
	///	</summary>
	// moab::ErrorCode GenerateMetaData();

public:

	///	<summary>
	///		Read the OfflineMap from a NetCDF file.
	///	</summary>
	// virtual void Read(
	// 	const std::string & strSource
	// );

	///	<summary>
	///		Write the TempestOnlineMap to a parallel NetCDF file.
	///	</summary>
	// virtual void Write(
	// 	const std::string & strTarget
	// );

public:
	///	<summary>
	///		Determine if the map is first-order accurate.
	///	</summary>
	virtual int IsConsistent(double dTolerance);

	///	<summary>
	///		Determine if the map is conservative.
	///	</summary>
	virtual int IsConservative(double dTolerance);

	///	<summary>
	///		Determine if the map is monotone.
	///	</summary>
	virtual int IsMonotone(double dTolerance);

	///	<summary>
	///		If we computed the reduction, get the vector representing the source areas for all entities in the mesh
	///	</summary>
	const DataArray1D<double>& GetGlobalSourceAreas() const;

	///	<summary>
	///		If we computed the reduction, get the vector representing the target areas for all entities in the mesh
	///	</summary>
	const DataArray1D<double>& GetGlobalTargetAreas() const;

private:

	///	<summary>
	///		Compute the remapping weights as a permutation matrix that relates DoFs on the source mesh
	///     to DoFs on the target mesh.
	///	</summary>
	moab::ErrorCode LinearRemapNN_MOAB( bool use_GID_matching=false, bool strict_check=false );

	///	<summary>
	///		Compute the remapping weights for a FV field defined on the source to a 
	///     FV field defined on the target mesh.
	///	</summary>
	void LinearRemapFVtoFV_Tempest_MOAB( int nOrder );

	///	<summary>
	///		Generate the OfflineMap for linear conserative element-average
	///		spectral element to element average remapping.
	///	</summary>
	void LinearRemapSE0_Tempest_MOAB(
		const DataArray3D<int> & dataGLLNodes,
		const DataArray3D<double> & dataGLLJacobian);

	///	<summary>
	///		Generate the OfflineMap for cubic conserative element-average
	///		spectral element to element average remapping.
	///	</summary>
	void LinearRemapSE4_Tempest_MOAB(
	    const DataArray3D<int> & dataGLLNodes,
	    const DataArray3D<double> & dataGLLJacobian,
	    int nMonotoneType,
	    bool fContinuousIn,
	    bool fNoConservation);


	///	<summary>
	///		Generate the OfflineMap for remapping from finite volumes to finite
	///		elements.
	///	</summary>
	void LinearRemapFVtoGLL_MOAB(
		const DataArray3D<int> & dataGLLNodes,
		const DataArray3D<double> & dataGLLJacobian,
		const DataArray1D<double> & dataGLLNodalArea,
		int nOrder,
		int nMonotoneType,
		bool fContinuous,
		bool fNoConservation
	);

	///	<summary>
	///		Generate the OfflineMap for remapping from finite elements to finite
	///		elements.
	///	</summary>
	void LinearRemapGLLtoGLL2_MOAB(
		const DataArray3D<int> & dataGLLNodesIn,
		const DataArray3D<double> & dataGLLJacobianIn,
		const DataArray3D<int> & dataGLLNodesOut,
		const DataArray3D<double> & dataGLLJacobianOut,
		const DataArray1D<double> & dataNodalAreaOut,
		int nPin,
		int nPout,
		int nMonotoneType,
		bool fContinuousIn,
		bool fContinuousOut,
		bool fNoConservation
	);

	///	<summary>
	///		Generate the OfflineMap for remapping from finite elements to finite
	///		elements (pointwise interpolation).
	///	</summary>
	void LinearRemapGLLtoGLL2_Pointwise_MOAB(
		const DataArray3D<int> & dataGLLNodesIn,
		const DataArray3D<double> & dataGLLJacobianIn,
		const DataArray3D<int> & dataGLLNodesOut,
		const DataArray3D<double> & dataGLLJacobianOut,
		const DataArray1D<double> & dataNodalAreaOut,
		int nPin,
		int nPout,
		int nMonotoneType,
		bool fContinuousIn,
		bool fContinuousOut
	);

	///	<summary>
	///		Copy the local matrix from Tempest SparseMatrix representation (ELL)
	///		to the parallel CSR Eigen Matrix for scalable application of matvec
	///     needed for projections.
	///	</summary>
#ifdef MOAB_HAVE_EIGEN
	void copy_tempest_sparsemat_to_eigen3();
#endif

public:

	///	<summary>
	///		Store the tag names associated with global DoF ids for source and target meshes
	///	</summary>
	moab::ErrorCode SetDOFmapTags(const std::string srcDofTagName, 
								  const std::string tgtDofTagName);

	///	<summary>
	///		Compute the association between the solution tag global DoF numbering and
	///		the local matrix numbering so that matvec operations can be performed
	///     consistently.
	///	</summary>
	moab::ErrorCode SetDOFmapAssociation(DiscretizationType srcType, bool isSrcContinuous, 
		DataArray3D<int>* srcdataGLLNodes, DataArray3D<int>* srcdataGLLNodesSrc,
		DiscretizationType destType, bool isDestContinuous, 
		DataArray3D<int>* tgtdataGLLNodes);

#ifdef MOAB_HAVE_EIGEN

	typedef Eigen::Matrix< double, 1, Eigen::Dynamic > WeightDRowVector;
	typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > WeightDColVector;
	typedef Eigen::SparseVector<double> WeightSVector;
	typedef Eigen::SparseMatrix<double, Eigen::RowMajor> WeightRMatrix;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor> WeightCMatrix;

	typedef WeightDRowVector WeightRowVector;
	typedef WeightDColVector WeightColVector;
	typedef WeightRMatrix WeightMatrix;

	///	<summary>
	///		Get the raw reference to the Eigen weight matrix representing the projection from source to destination mesh.
	///	</summary>
	WeightMatrix& GetWeightMatrix();

	///	<summary>
	///		Get the row vector that is amenable for application of A*x operation.
	///	</summary>
	WeightRowVector& GetRowVector();

	///	<summary>
	///		Get the column vector that is amenable for application of A^T*x operation.
	///	</summary>
	WeightColVector& GetColVector();

#endif

	///	<summary>
	///		Get the number of total Degrees-Of-Freedom defined on the source mesh.
	///	</summary>
	int GetSourceGlobalNDofs();

	///	<summary>
	///		Get the number of total Degrees-Of-Freedom defined on the destination mesh.
	///	</summary>
	int GetDestinationGlobalNDofs();

	///	<summary>
	///		Get the number of local Degrees-Of-Freedom defined on the source mesh.
	///	</summary>
	int GetSourceLocalNDofs();

	///	<summary>
	///		Get the number of local Degrees-Of-Freedom defined on the destination mesh.
	///	</summary>
	int GetDestinationLocalNDofs();

	///	<summary>
	///		Get the number of Degrees-Of-Freedom per element on the source mesh.
	///	</summary>
	int GetSourceNDofsPerElement();

	///	<summary>
	///		Get the number of Degrees-Of-Freedom per element on the destination mesh.
	///	</summary>
	int GetDestinationNDofsPerElement();

	///	<summary>
	///		Get the global Degrees-Of-Freedom ID on the destination mesh.
	///	</summary>
    int GetRowGlobalDoF(int localID) const;

	///	<summary>
	///		Get the global Degrees-Of-Freedom ID on the source mesh.
	///	</summary>
    int GetColGlobalDoF(int localID) const;

	///	<summary>
	///		Apply the weight matrix onto the source vector provided as input, and return the column vector (solution projection) after the map application
	///     Compute:        \p tgtVals = A(S->T) * \srcVals, or
	///     if (transpose)  \p tgtVals = [A(T->S)]^T * \srcVals
	///	</summary>
	moab::ErrorCode ApplyWeights (std::vector<double>& srcVals, std::vector<double>& tgtVals, bool transpose=false);

	///	<summary>
	///		Apply the weight matrix onto the source vector (tag) provided as input, and return the column vector (solution projection) in a tag, after the map application
	///     Compute:        \p tgtVals = A(S->T) * \srcVals, or
	///     if (transpose)  \p tgtVals = [A(T->S)]^T * \srcVals
	///	</summary>
	moab::ErrorCode ApplyWeights (  moab::Tag srcSolutionTag, moab::Tag tgtSolutionTag, bool transpose=false );

	///	<summary>
	///		Parallel I/O with NetCDF to write out the SCRIP file from multiple processors.
	///	</summary>
	void WriteParallelWeightsToFile(std::string filename);

	///	<summary>
	///		Parallel I/O with HDF5 to write out the remapping weights from multiple processors.
	///	</summary>
	moab::ErrorCode WriteParallelMap (std::string strOutputFile);

    typedef double (*sample_function)(double, double);

    /// <summary>
    ///     Define an analytical solution over the given (source or target) mesh, as specificed in the context.
    ///     This routine will define a tag that is compatible with the specified discretization method type and order
    ///     and sample the solution exactly using the analytical function provided by the user.
    /// </summary>
    moab::ErrorCode DefineAnalyticalSolution ( moab::Tag& exactSolnTag, const std::string& solnName, Remapper::IntersectionContext ctx, 
                                                sample_function testFunction,
                                                moab::Tag* clonedSolnTag=NULL,
                                                std::string cloneSolnName="");

    /// <summary>
    ///     Compute the error between a sampled (exact) solution and a projected solution in various error norms.
    /// </summary>
	moab::ErrorCode ComputeMetrics ( Remapper::IntersectionContext ctx, moab::Tag& exactTag, moab::Tag& approxTag, 
									 std::map<std::string, double>& metrics, bool verbose = true );

public:

	///	<summary>
	///		The fundamental remapping operator object.
	///	</summary>
	moab::TempestRemapper* m_remapper;

#ifdef MOAB_HAVE_EIGEN

	WeightMatrix m_weightMatrix;
	WeightRowVector m_rowVector;
	WeightColVector m_colVector;

#endif

	///	<summary>
	///		The reference to the moab::Core object that contains source/target and overlap sets.
	///	</summary>
	moab::Interface* m_interface;

#ifdef MOAB_HAVE_MPI
	///	<summary>
	///		The reference to the parallel communicator object used by the Core object.
	///	</summary>
	moab::ParallelComm* m_pcomm;
#endif

	///	<summary>
	///		The original tag data and local to global DoF mapping to associate matrix values to solution
	///	<summary>
	moab::Tag m_dofTagSrc, m_dofTagDest;
	std::vector<unsigned> row_gdofmap, col_gdofmap, srccol_gdofmap;
	std::vector<unsigned> row_dtoc_dofmap, col_dtoc_dofmap, srccol_dtoc_dofmap;

	DataArray3D<int> dataGLLNodesSrc, dataGLLNodesSrcCov, dataGLLNodesDest;
	DiscretizationType m_srcDiscType, m_destDiscType;
	int m_nTotDofs_Src, m_nTotDofs_SrcCov, m_nTotDofs_Dest;
	
	// Key details about the current map
	int m_nDofsPEl_Src, m_nDofsPEl_Dest;
	DiscretizationType m_eInputType, m_eOutputType;
	bool m_bConserved;
	int m_iMonotonicity;

	Mesh* m_meshInput;
	Mesh* m_meshInputCov;
	Mesh* m_meshOutput;
	Mesh* m_meshOverlap;
	
	bool is_parallel, is_root;
	int rank, size;
};


///////////////////////////////////////////////////////////////////////////////

inline
int moab::TempestOnlineMap::GetRowGlobalDoF(int localRowID) const
{
    return row_gdofmap [ localRowID ];
}

///////////////////////////////////////////////////////////////////////////////

inline
int moab::TempestOnlineMap::GetColGlobalDoF(int localColID) const
{
    return col_gdofmap [ localColID ];
}


///////////////////////////////////////////////////////////////////////////////
#ifdef MOAB_HAVE_EIGEN

inline
int moab::TempestOnlineMap::GetSourceGlobalNDofs()
{
    return m_weightMatrix.cols(); // return the global number of rows from the weight matrix
}

// ///////////////////////////////////////////////////////////////////////////////

inline
int moab::TempestOnlineMap::GetDestinationGlobalNDofs()
{
    return m_weightMatrix.rows(); // return the global number of columns from the weight matrix
}

///////////////////////////////////////////////////////////////////////////////

inline
int moab::TempestOnlineMap::GetSourceLocalNDofs()
{
    return m_weightMatrix.cols(); // return the local number of rows from the weight matrix
}

///////////////////////////////////////////////////////////////////////////////

inline
int moab::TempestOnlineMap::GetDestinationLocalNDofs()
{
    return m_weightMatrix.rows(); // return the local number of columns from the weight matrix
}

///////////////////////////////////////////////////////////////////////////////

inline
int moab::TempestOnlineMap::GetSourceNDofsPerElement()
{
    return m_nDofsPEl_Src;
}

///////////////////////////////////////////////////////////////////////////////

inline
int moab::TempestOnlineMap::GetDestinationNDofsPerElement()
{
    return m_nDofsPEl_Dest;
}

///////////////////////////////////////////////////////////////////////////////

inline
moab::TempestOnlineMap::WeightMatrix& moab::TempestOnlineMap::GetWeightMatrix()
{
    assert(m_weightMatrix.rows() != 0 && m_weightMatrix.cols() != 0);
    return m_weightMatrix;
}

///////////////////////////////////////////////////////////////////////////////

inline
moab::TempestOnlineMap::WeightRowVector& moab::TempestOnlineMap::GetRowVector()
{
    assert(m_rowVector.size() != 0);
    return m_rowVector;
}

///////////////////////////////////////////////////////////////////////////////

inline
moab::TempestOnlineMap::WeightColVector& moab::TempestOnlineMap::GetColVector()
{
    assert(m_colVector.size() != 0);
    return m_colVector;
}

///////////////////////////////////////////////////////////////////////////////

#endif

}

#endif
