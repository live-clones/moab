/*
 * Usage: MOAB intersection check tool, for intersection on a sphere or in plane
 *
 * mpiexec -np n ./mbintx_check -s <source mesh> -t <target mesh> -i <intx mesh>  -o <source verif> -q <target verif>
 *
 * after a run with mbtempest in parallel, that computes intersection, we can use this tool to verify areas
 * of intersection polygons, compared to areas of source, target elements
 *
 */

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>


#include "moab/Core.hpp"
#include "moab/IntxMesh/IntxUtils.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CpuTimer.hpp"
#include "DebugOutput.hpp"

#ifdef MOAB_HAVE_MPI
// MPI includes
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

using namespace moab;

int main ( int argc, char* argv[] )
{
    std::stringstream sstr;

    int rank = 0, size = 1;
#ifdef MOAB_HAVE_MPI
    MPI_Init ( &argc, &argv );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &size );
#endif

    std::string sourceFile, targetFile, intxFile;
    std::string source_verif("outS.h5m"), target_verif("outt.h5m");
    int sphere = 1;
    int oldNamesParents = 1;
    double areaErrSource = -1;
    double areaErrTarget = -1;
    ProgOptions opts;

    opts.addOpt<std::string> ( "source,s", "source file ", &sourceFile );
    opts.addOpt<std::string> ( "target,t", "target file ", &targetFile );
    opts.addOpt<std::string> ( "intersection,i", "intersection file ", &intxFile );
    opts.addOpt<std::string> ( "verif_source,v", "output source verification ", &source_verif );
    opts.addOpt<std::string> ( "verif_target,w", "output target verification ", &target_verif );
    opts.addOpt<double> ( "threshold_source,m", "error source threshold ", &areaErrSource );
    opts.addOpt<double> ( "threshold_target,q", "error target threshold ", &areaErrTarget );

    opts.addOpt<int> ( "sphere,p", "mesh on a sphere", &sphere );
    opts.addOpt<int> ( "old_convention,n", "old names for parent tags", &oldNamesParents );

    opts.parseCommandLine ( argc, argv );
    // load meshes in parallel if needed
    std::string opts_read = (size == 1 ? "" : std::string("PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION")+
          std::string(";PARALLEL_RESOLVE_SHARED_ENTS"));

    // read meshes in 3 file sets
    ErrorCode rval;
    Core moab;
    Interface * mb = &moab; // global
    EntityHandle sset, tset, ixset;

    // create meshsets and load files
    if (0==rank)
      std::cout << "Creating mesh sets\n";
    rval = mb->create_meshset(MESHSET_SET, sset);MB_CHK_ERR(rval);
    rval = mb->create_meshset(MESHSET_SET, tset);MB_CHK_ERR(rval);
    rval = mb->create_meshset(MESHSET_SET, ixset);MB_CHK_ERR(rval);
    if (0==rank)
      std::cout << "Loading source file \n";
    rval = mb->load_file(sourceFile.c_str(), &sset, opts_read.c_str()); MB_CHK_SET_ERR(rval, "failed reading source file");
    if (0==rank)
      std::cout << "Loading target file \n";
    rval = mb->load_file(targetFile.c_str(), &tset, opts_read.c_str()); MB_CHK_SET_ERR(rval, "failed reading target file");

    if (0==rank)
      std::cout << "Loading intersection file \n";
    rval = mb->load_file(intxFile.c_str(), &ixset, opts_read.c_str()); MB_CHK_SET_ERR(rval, "failed reading intersection file");
    double R = 1.;
    if (sphere){
      // fix radius of both meshes, to be consistent with radius 1
      rval = ScaleToRadius(mb, sset, R); MB_CHK_ERR(rval);
      rval = ScaleToRadius(mb, tset, R); MB_CHK_ERR(rval);
    }
    Range intxCells;
    rval = mb->get_entities_by_dimension(ixset, 2, intxCells);MB_CHK_ERR(rval);

    Range sourceCells;
    rval = mb->get_entities_by_dimension(sset, 2, sourceCells);MB_CHK_ERR(rval);

    Range targetCells;
    rval = mb->get_entities_by_dimension(tset, 2, targetCells);MB_CHK_ERR(rval);

    Tag sourceParentTag;
    Tag targetParentTag;
    if (oldNamesParents)
    {
      rval = mb->tag_get_handle("RedParent", targetParentTag);MB_CHK_SET_ERR(rval, "can't find target parent tag");
      rval =  mb->tag_get_handle("BlueParent", sourceParentTag);MB_CHK_SET_ERR(rval, "can't find source parent tag");
    }
    else {
      rval = mb->tag_get_handle("TargetParent", targetParentTag);MB_CHK_SET_ERR(rval, "can't find target parent tag");
      rval =  mb->tag_get_handle("SourceParent", sourceParentTag);MB_CHK_SET_ERR(rval, "can't find source parent tag");
    }


    // error sets, for better visualization
    EntityHandle errorSourceSet;
    rval = mb->create_meshset(MESHSET_SET, errorSourceSet);MB_CHK_ERR(rval);
    EntityHandle errorTargetSet;
    rval = mb->create_meshset(MESHSET_SET, errorTargetSet);MB_CHK_ERR(rval);

    std::map<int, double> sourceAreas;
    std::map<int, double> targetAreas;

    std::map<int, double> sourceAreasIntx;
    std::map<int, double> targetAreasIntx;

    Tag gidTag = mb->globalId_tag();

    Tag areaTag;
    rval = mb->tag_get_handle("OrigArea", 1, MB_TYPE_DOUBLE, areaTag, MB_TAG_DENSE | MB_TAG_CREAT);MB_CHK_ERR(rval);

    for (Range::iterator eit = sourceCells.begin(); eit != sourceCells.end(); ++eit)
    {
      EntityHandle cell = *eit;
      const EntityHandle * verts;
      int num_nodes;
      rval = mb->get_connectivity(cell, verts, num_nodes);MB_CHK_ERR(rval);
      if (MB_SUCCESS != rval)
        return -1;
      std::vector<double> coords(3 * num_nodes);
      // get coordinates
      rval = mb->get_coords(verts, num_nodes, &coords[0]);
      if (MB_SUCCESS != rval)
        return -1;
      double area = area_spherical_polygon_lHuiller(&coords[0], num_nodes, R);
      int sourceID;
      rval = mb->tag_get_data(gidTag, &cell, 1, &sourceID);MB_CHK_ERR(rval);
      sourceAreas[sourceID] = area;
      rval = mb->tag_set_data(areaTag, &cell, 1, &area);MB_CHK_ERR(rval);
    }
    for (Range::iterator eit = targetCells.begin(); eit != targetCells.end(); ++eit)
    {
      EntityHandle cell = *eit;
      const EntityHandle * verts;
      int num_nodes;
      rval = mb->get_connectivity(cell, verts, num_nodes);MB_CHK_ERR(rval);
      if (MB_SUCCESS != rval)
        return -1;
      std::vector<double> coords(3 * num_nodes);
      // get coordinates
      rval = mb->get_coords(verts, num_nodes, &coords[0]);
      if (MB_SUCCESS != rval)
        return -1;
      double area = area_spherical_polygon_lHuiller(&coords[0], num_nodes, R);
      int targetID;
      rval = mb->tag_get_data(gidTag, &cell, 1, &targetID);MB_CHK_ERR(rval);
      targetAreas[targetID] = area;
      rval = mb->tag_set_data(areaTag, &cell, 1, &area);MB_CHK_ERR(rval);
    }

    for (Range::iterator eit = intxCells.begin(); eit != intxCells.end(); ++eit)
    {
      EntityHandle cell = *eit;
      const EntityHandle * verts;
      int num_nodes;
      rval = mb->get_connectivity(cell, verts, num_nodes);MB_CHK_ERR(rval);
      if (MB_SUCCESS != rval)
        return -1;
      std::vector<double> coords(3 * num_nodes);
      // get coordinates
      rval = mb->get_coords(verts, num_nodes, &coords[0]);
      if (MB_SUCCESS != rval)
        return -1;
      double intx_area = area_spherical_polygon_lHuiller(&coords[0], num_nodes, R);
      int sourceID, targetID;
      rval = mb->tag_get_data(sourceParentTag, &cell, 1, &sourceID);MB_CHK_ERR(rval);
      rval = mb->tag_get_data(targetParentTag, &cell, 1, &targetID);MB_CHK_ERR(rval);

      std::map<int, double>::iterator sit=sourceAreasIntx.find(sourceID);
      if (sit==sourceAreasIntx.end())
      {
        sourceAreasIntx[sourceID] = intx_area;
      }
      else
        sourceAreasIntx[sourceID] += intx_area;

      std::map<int, double>::iterator tit=targetAreasIntx.find(targetID);
      if (tit==targetAreasIntx.end())
      {
        targetAreasIntx[targetID] = intx_area;
      }
      else
        targetAreasIntx[targetID] += intx_area;
    }
    Tag diffTag;
    rval = mb->tag_get_handle("AreaDiff", 1, MB_TYPE_DOUBLE, diffTag, MB_TAG_DENSE | MB_TAG_CREAT);MB_CHK_ERR(rval);


    for (Range::iterator eit = sourceCells.begin(); eit != sourceCells.end(); ++eit)
    {
      EntityHandle cell = *eit;

      int sourceID;
      rval = mb->tag_get_data(gidTag, &cell, 1, &sourceID);MB_CHK_ERR(rval);
      double areaDiff = sourceAreas[sourceID];
      std::map<int, double>::iterator sit=sourceAreasIntx.find(sourceID);
      if (sit!=sourceAreasIntx.end())
      {
        areaDiff -= sourceAreasIntx[sourceID];
      }
      rval = mb->tag_set_data(diffTag, &cell, 1, &areaDiff);
      // add to errorSourceSet set if needed
      if ( (areaErrSource > 0) && (fabs(areaDiff) > areaErrSource))
      {
        rval = mb->add_entities(errorSourceSet, &cell, 1);MB_CHK_ERR(rval);
      }
    }
    rval = mb->write_file( source_verif.c_str(),0, 0,&sset, 1);MB_CHK_ERR(rval);
    if (areaErrSource > 0)
    {
      std::string filterSource = std::string("filt_")+source_verif;
      rval = mb->write_file( filterSource.c_str(),0, 0,&errorSourceSet, 1);
    }


    for (Range::iterator eit = targetCells.begin(); eit != targetCells.end(); ++eit)
    {
      EntityHandle cell = *eit;

      int targetID;
      rval = mb->tag_get_data(gidTag, &cell, 1, &targetID);MB_CHK_ERR(rval);
      double areaDiff = targetAreas[targetID];
      std::map<int, double>::iterator sit=targetAreasIntx.find(targetID);
      if (sit!=targetAreasIntx.end())
      {
        areaDiff -= targetAreasIntx[targetID];
      }

      rval = mb->tag_set_data(diffTag, &cell, 1, &areaDiff);
      // add to errorTargetSet set if needed
      if ( (areaErrTarget > 0) && (fabs(areaDiff) > areaErrTarget))
      {
        rval = mb->add_entities(errorTargetSet, &cell, 1);MB_CHK_ERR(rval);
      }

    }
    rval = mb->write_file(target_verif.c_str(), 0, 0, &tset, 1);MB_CHK_ERR(rval);
    if (areaErrTarget > 0)
    {
      std::string filterTarget = std::string("filt_")+target_verif;
      rval = mb->write_file( filterTarget.c_str(),0, 0,&errorTargetSet, 1);
    }

    return 0;
}
