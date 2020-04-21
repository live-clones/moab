/*
 * Usage: MOAB-Tempest tool
 *
 * Generate a Cubed-Sphere mesh: ./mbtempest -t 0 -res 25 -f cubed_sphere_mesh.exo
 * Generate a RLL mesh: ./mbtempest -t 1 -res 25 -f rll_mesh.exo
 * Generate a Icosahedral-Sphere mesh: ./mbtempest -t 2 -res 25 <-dual> -f icosa_mesh.exo
 *
 * Now you can compute the intersections between the meshes too!
 *
 * Generate the overlap mesh: ./mbtempest -t 3 -l cubed_sphere_mesh.exo -l rll_mesh.exo -f overlap_mesh.exo
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
#include "moab/Remapping/TempestRemapper.hpp"
#include "moab/Remapping/TempestOnlineMap.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CpuTimer.hpp"
#include "DebugOutput.hpp"

//#ifndef MOAB_HAVE_MPI
//    #error mbtempest tool requires MPI configuration
//#endif

#ifdef MOAB_HAVE_MPI
// MPI includes
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

struct ToolContext
{
        moab::Interface* mbcore;
#ifdef MOAB_HAVE_MPI
        moab::ParallelComm* pcomm;
#endif
        const int proc_id, n_procs;
        moab::DebugOutput outputFormatter;
        int blockSize;
        std::vector<std::string> inFilenames;
        std::vector<Mesh*> meshes;
        std::vector<moab::EntityHandle> meshsets;
        std::vector<int> disc_orders;
        std::vector<std::string> disc_methods;
        std::vector<std::string> doftag_names;
        std::string outFilename;
        std::string intxFilename;
        moab::TempestRemapper::TempestMeshType meshType;
        bool computeDual;
        bool computeWeights;
        bool verifyWeights;
        int ensureMonotonicity;
        bool fNoConservation;
        bool fVolumetric;
        bool rrmGrids;
        bool kdtreeSearch;
        bool fNoBubble, fInputConcave, fOutputConcave, fNoCheck;

#ifdef MOAB_HAVE_MPI
        ToolContext ( moab::Interface* icore, moab::ParallelComm* p_pcomm ) :
            mbcore(icore), pcomm(p_pcomm),
            proc_id ( pcomm->rank() ), n_procs ( pcomm->size() ),
            outputFormatter ( std::cout, pcomm->rank(), 0 ),
#else
        ToolContext ( moab::Interface* icore ) :
            mbcore(icore),
            proc_id ( 0 ), n_procs ( 1 ),
            outputFormatter ( std::cout, 0, 0 )
#endif
            blockSize ( 5 ), outFilename ( "output.exo" ), intxFilename ( "" ), meshType ( moab::TempestRemapper::DEFAULT ),
            computeDual ( false ), computeWeights ( false ), verifyWeights ( false ), ensureMonotonicity ( 0 ), 
            fNoConservation ( false ), fVolumetric ( false ), rrmGrids ( false ), kdtreeSearch ( true ),
            fNoBubble(false), fInputConcave(false), fOutputConcave(false), fNoCheck(false)
        {
            inFilenames.resize ( 2 );
            doftag_names.resize( 2 );
            timer = new moab::CpuTimer();

            outputFormatter.set_prefix("[MBTempest]: ");
        }

        ~ToolContext()
        {
            // for (unsigned i=0; i < meshes.size(); ++i) delete meshes[i];
            meshes.clear();
            inFilenames.clear();
            disc_orders.clear();
            disc_methods.clear();
            doftag_names.clear();
            outFilename.clear();
            intxFilename.clear();
            meshsets.clear();
            delete timer;
        }

        void timer_push ( std::string operation )
        {
            timer_ops = timer->time_since_birth();
            opName = operation;
        }

        void timer_pop()
        {
            double locElapsed=timer->time_since_birth() - timer_ops, avgElapsed=0, maxElapsed=0;
#ifdef MOAB_HAVE_MPI
            MPI_Reduce(&locElapsed, &maxElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, pcomm->comm());
            MPI_Reduce(&locElapsed, &avgElapsed, 1, MPI_DOUBLE, MPI_SUM, 0, pcomm->comm());
#else
            maxElapsed = locElapsed;
            avgElapsed = locElapsed;
#endif
            if (!proc_id) {
                avgElapsed /= n_procs;
                std::cout << "[LOG] Time taken to " << opName.c_str() << ": max = " << maxElapsed << ", avg = " << avgElapsed << "\n";
            }
            // std::cout << "\n[LOG" << proc_id << "] Time taken to " << opName << " = " << timer->time_since_birth() - timer_ops << std::endl;
            opName.clear();
        }

        void ParseCLOptions ( int argc, char* argv[] )
        {
            ProgOptions opts;
            int imeshType = 0;
            std::string expectedFName = "output.exo";
            std::string expectedMethod = "fv";
            std::string expectedDofTagName = "GLOBAL_ID";
            int expectedOrder = 1;

            if (!proc_id)
            {
                std::cout << "Command line options provided to mbtempest:\n  ";
                for(int iarg=0; iarg < argc; ++iarg)
                    std::cout << argv[iarg] << " ";
                std::cout << std::endl << std::endl;
            }

            opts.addOpt<int> ( "type,t", "Type of mesh (default=CS; Choose from [CS=0, RLL=1, ICO=2, OVERLAP_FILES=3, OVERLAP_MEMORY=4, OVERLAP_MOAB=5])", &imeshType );
            opts.addOpt<int> ( "res,r", "Resolution of the mesh (default=5)", &blockSize );
            opts.addOpt<void> ( "dual,d", "Output the dual of the mesh (generally relevant only for ICO mesh)", &computeDual );
            opts.addOpt<void> ( "weights,w", "Compute and output the weights using the overlap mesh (generally relevant only for OVERLAP mesh)", &computeWeights );
            opts.addOpt<void> ( "noconserve,c", "Do not apply conservation to the resultant weights (relevant only when computing weights)", &fNoConservation );
            opts.addOpt<void> ( "volumetric,v", "Apply a volumetric projection to compute the weights (relevant only when computing weights)", &fVolumetric );
            opts.addOpt<void> ( "rrmgrids", "At least one of the meshes is a regionally refined grid (relevant to accelerate intersection computation)", &rrmGrids );
            opts.addOpt<void> ( "nocheck", "Do not check the generated map for conservation and consistency", &fNoCheck );
            opts.addOpt<void> ( "advfront,a", "Use the advancing front intersection instead of the Kd-tree based algorithm to compute mesh intersections" );
            opts.addOpt<void> ( "verify", "Verify the accuracy of the maps by projecting analytical functions from source to target grid by applying the maps", &verifyWeights );
            opts.addOpt<int> ( "monotonic,n", "Ensure monotonicity in the weight generation", &ensureMonotonicity );
            opts.addOpt<std::string> ( "load,l", "Input mesh filenames (a source and target mesh)", &expectedFName );
            opts.addOpt<int> ( "order,o", "Discretization orders for the source and target solution fields", &expectedOrder );
            opts.addOpt<std::string> ( "method,m", "Discretization method for the source and target solution fields", &expectedMethod );
            opts.addOpt<std::string> ( "global_id,g", "Tag name that contains the global DoF IDs for source and target solution fields", &expectedDofTagName );
            opts.addOpt<std::string> ( "file,f", "Output remapping weights filename", &outFilename );
            opts.addOpt<std::string> ( "intx,i", "Output TempestRemap intersection mesh filename", &intxFilename );

            opts.parseCommandLine ( argc, argv );

            // By default - use Kd-tree based search; if user asks for advancing front, disable Kd-tree algorithm
            kdtreeSearch = opts.numOptSet("advfront,a") == 0;

            switch ( imeshType )
            {
                case 0:
                    meshType = moab::TempestRemapper::CS;
                    break;

                case 1:
                    meshType = moab::TempestRemapper::RLL;
                    break;

                case 2:
                    meshType = moab::TempestRemapper::ICO;
                    break;

                case 3:
                    meshType = moab::TempestRemapper::OVERLAP_FILES;
                    break;

                case 4:
                    meshType = moab::TempestRemapper::OVERLAP_MEMORY;
                    break;

                case 5:
                    meshType = moab::TempestRemapper::OVERLAP_MOAB;
                    break;

                default:
                    meshType = moab::TempestRemapper::DEFAULT;
                    break;
            }

            if ( meshType > moab::TempestRemapper::ICO )
            {
                opts.getOptAllArgs ( "load,l", inFilenames );
                opts.getOptAllArgs ( "order,o", disc_orders );
                opts.getOptAllArgs ( "method,m", disc_methods );
                opts.getOptAllArgs ( "global_id,i", doftag_names );

                if ( disc_orders.size() == 0 )
                { disc_orders.resize ( 2, 1 ); }

                if ( disc_orders.size() == 1 )
                { disc_orders.push_back ( 1 ); }

                if ( disc_methods.size() == 0 )
                { disc_methods.resize ( 2, "fv" ); }

                if ( disc_methods.size() == 1 )
                { disc_methods.push_back ( "fv" ); }

                if ( doftag_names.size() == 0 )
                { doftag_names.resize ( 2, "GLOBAL_ID" ); }

                if ( doftag_names.size() == 1 )
                { doftag_names.push_back ( "GLOBAL_ID" ); }

                assert ( inFilenames.size() == 2 );
                assert ( disc_orders.size() == 2 );
                assert ( disc_methods.size() == 2 );
            }

            expectedFName.clear();
        }
    private:
        moab::CpuTimer* timer;
        double timer_ops;
        std::string opName;
};

// Forward declare some methods
static moab::ErrorCode CreateTempestMesh ( ToolContext&, moab::TempestRemapper& remapper, Mesh* );
inline double sample_slow_harmonic ( double dLon, double dLat );
inline double sample_fast_harmonic ( double dLon, double dLat );
inline double sample_constant ( double dLon, double dLat );
inline double sample_stationary_vortex ( double dLon, double dLat );

int main ( int argc, char* argv[] )
{
    moab::ErrorCode rval;
    NcError error ( NcError::verbose_nonfatal );
    std::stringstream sstr;

    int proc_id = 0, nprocs = 1;
#ifdef MOAB_HAVE_MPI
    MPI_Init ( &argc, &argv );
    MPI_Comm_rank ( MPI_COMM_WORLD, &proc_id );
    MPI_Comm_size ( MPI_COMM_WORLD, &nprocs );
#endif

    moab::Interface* mbCore = new ( std::nothrow ) moab::Core;

    if ( NULL == mbCore ) { return 1; }

    ToolContext *runCtx;
#ifdef MOAB_HAVE_MPI
    moab::ParallelComm* pcomm = new moab::ParallelComm ( mbCore, MPI_COMM_WORLD, 0 );

    runCtx = new ToolContext( mbCore, pcomm );
    const char *writeOptions = (nprocs > 1 ? "PARALLEL=WRITE_PART" : "");
#else
    runCtx = new ToolContext( mbCore );
    const char *writeOptions = "";
#endif
    runCtx->ParseCLOptions ( argc, argv );

    const double radius_src = 1.0 /*2.0*acos(-1.0)*/;
    const double radius_dest = 1.0 /*2.0*acos(-1.0)*/;

    moab::DebugOutput& outputFormatter = runCtx->outputFormatter;

#ifdef MOAB_HAVE_MPI
    moab::TempestRemapper remapper ( mbCore, pcomm );
#else
    moab::TempestRemapper remapper ( mbCore );
#endif
    remapper.meshValidate = true;
    remapper.constructEdgeMap = false;
    remapper.initialize();

    Mesh* tempest_mesh = new Mesh();
    runCtx->timer_push ( "create Tempest mesh" );
    rval = CreateTempestMesh ( *runCtx, remapper, tempest_mesh ); MB_CHK_ERR ( rval );
    runCtx->timer_pop();

    // Some constant parameters
    const double epsrel = 1.e-15;
    const double boxeps = 1e-6;

    if ( runCtx->meshType == moab::TempestRemapper::OVERLAP_MEMORY )
    {
        // Compute intersections with MOAB
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        assert ( runCtx->meshes.size() == 3 );

#ifdef MOAB_HAVE_MPI
        rval = pcomm->check_all_shared_handles(); MB_CHK_ERR ( rval );
#endif

        // Load the meshes and validate
        rval = remapper.ConvertTempestMesh ( moab::Remapper::SourceMesh ); MB_CHK_ERR ( rval );
        rval = remapper.ConvertTempestMesh ( moab::Remapper::TargetMesh ); MB_CHK_ERR ( rval );
        remapper.SetMeshType ( moab::Remapper::OverlapMesh, moab::TempestRemapper::OVERLAP_FILES );
        rval = remapper.ConvertTempestMesh ( moab::Remapper::OverlapMesh ); MB_CHK_ERR ( rval );
        rval = mbCore->write_mesh ( "tempest_intersection.h5m", &runCtx->meshsets[2], 1 ); MB_CHK_ERR ( rval );

        // print verbosely about the problem setting
        {
            moab::Range rintxverts, rintxelems;
            rval = mbCore->get_entities_by_dimension ( runCtx->meshsets[0], 0, rintxverts ); MB_CHK_ERR ( rval );
            rval = mbCore->get_entities_by_dimension ( runCtx->meshsets[0], 2, rintxelems ); MB_CHK_ERR ( rval );
            outputFormatter.printf ( 0,  "The source set contains %lu vertices and %lu elements \n", rintxverts.size(), rintxelems.size() );

            moab::Range bintxverts, bintxelems;
            rval = mbCore->get_entities_by_dimension ( runCtx->meshsets[1], 0, bintxverts ); MB_CHK_ERR ( rval );
            rval = mbCore->get_entities_by_dimension ( runCtx->meshsets[1], 2, bintxelems ); MB_CHK_ERR ( rval );
            outputFormatter.printf ( 0,  "The target set contains %lu vertices and %lu elements \n", bintxverts.size(), bintxelems.size() );
        }

        moab::EntityHandle intxset; // == remapper.GetMeshSet(moab::Remapper::OverlapMesh);

        // Compute intersections with MOAB
        {
            // Create the intersection on the sphere object
            runCtx->timer_push ( "setup the intersector" );

            moab::Intx2MeshOnSphere* mbintx = new moab::Intx2MeshOnSphere ( mbCore );
            mbintx->set_error_tolerance ( epsrel );
            mbintx->set_box_error ( boxeps );
            mbintx->set_radius_source_mesh ( radius_src );
            mbintx->set_radius_destination_mesh ( radius_dest );
#ifdef MOAB_HAVE_MPI
            mbintx->set_parallel_comm ( pcomm );
#endif
            rval = mbintx->FindMaxEdges ( runCtx->meshsets[0], runCtx->meshsets[1] ); MB_CHK_ERR ( rval );

#ifdef MOAB_HAVE_MPI
            moab::Range local_verts;
            rval = mbintx->build_processor_euler_boxes ( runCtx->meshsets[1], local_verts ); MB_CHK_ERR ( rval );

            runCtx->timer_pop();

            moab::EntityHandle covering_set;
            runCtx->timer_push ( "communicate the mesh" );
            rval = mbintx->construct_covering_set ( runCtx->meshsets[0], covering_set ); MB_CHK_ERR ( rval ); // lots of communication if mesh is distributed very differently
            runCtx->timer_pop();
#else
            moab::EntityHandle covering_set = runCtx->meshsets[0];
#endif
            // Now let's invoke the MOAB intersection algorithm in parallel with a
            // source and target mesh set representing two different decompositions
            runCtx->timer_push ( "compute intersections with MOAB" );
            rval = mbCore->create_meshset ( moab::MESHSET_SET, intxset ); MB_CHK_SET_ERR ( rval, "Can't create new set" );
            rval = mbintx->intersect_meshes ( covering_set, runCtx->meshsets[1], intxset ); MB_CHK_SET_ERR ( rval, "Can't compute the intersection of meshes on the sphere" );
            runCtx->timer_pop();

            // free the memory
            delete mbintx;
        }

        {
            moab::Range intxelems, intxverts;
            rval = mbCore->get_entities_by_dimension ( intxset, 2, intxelems ); MB_CHK_ERR ( rval );
            rval = mbCore->get_entities_by_dimension ( intxset, 0, intxverts, true ); MB_CHK_ERR ( rval );
            outputFormatter.printf ( 0,  "The intersection set contains %lu elements and %lu vertices \n", intxelems.size(), intxverts.size() );

            double initial_sarea = area_on_sphere_lHuiller ( mbCore, runCtx->meshsets[0], radius_src ); // use the target to compute the initial area
            double initial_tarea = area_on_sphere_lHuiller ( mbCore, runCtx->meshsets[1], radius_dest ); // use the target to compute the initial area
            double area_method1 = area_on_sphere_lHuiller ( mbCore, intxset, radius_src );
            double area_method2 = area_on_sphere ( mbCore, intxset, radius_src );

            outputFormatter.printf ( 0,  "initial areas: source = %12.10f, target = %12.10f \n", initial_sarea, initial_tarea );
            outputFormatter.printf ( 0,  " area with l'Huiller: %12.10f with Girard: %12.10f\n", area_method1, area_method2 );
            outputFormatter.printf ( 0,  " relative difference areas = %12.10e\n", fabs ( area_method1 - area_method2 ) / area_method1 );
            outputFormatter.printf ( 0,  " relative error w.r.t source = %12.10e, target = %12.10e \n", fabs ( area_method1 - initial_sarea ) / area_method1, fabs ( area_method1 - initial_tarea ) / area_method1 );
        }

        // Write out our computed intersection file
        rval = mbCore->write_mesh ( "moab_intersection.h5m", &intxset, 1 ); MB_CHK_ERR ( rval );

        if ( runCtx->computeWeights )
        {
            runCtx->timer_push ( "compute weights with the Tempest meshes" );
            // Call to generate an offline map with the tempest meshes
            OfflineMap weightMap;
            int err = GenerateOfflineMapWithMeshes (  weightMap, *runCtx->meshes[0], *runCtx->meshes[1], *runCtx->meshes[2],
                      "", "",     // std::string strInputMeta, std::string strOutputMeta,
                      runCtx->disc_methods[0], runCtx->disc_methods[1], // std::string strInputType, std::string strOutputType,
                      runCtx->disc_orders[0], runCtx->disc_orders[1],  // int nPin=4, int nPout=4,
                      runCtx->fNoBubble, true, runCtx->ensureMonotonicity // bool fNoBubble = false, bool fCorrectAreas = false, int fMonotoneTypeID = 0
                                                   );
            runCtx->timer_pop();

            std::map<std::string, std::string> mapAttributes;
            if ( err ) { rval = moab::MB_FAILURE; }
            else { weightMap.Write ( "outWeights.nc", mapAttributes ); }
        }
    }
    else if ( runCtx->meshType == moab::TempestRemapper::OVERLAP_MOAB )
    {
        // Usage: mpiexec -n 2 tools/mbtempest -t 5 -l mycs_2.h5m -l myico_2.h5m -f myoverlap_2.h5m
#ifdef MOAB_HAVE_MPI
        rval = pcomm->check_all_shared_handles(); MB_CHK_ERR ( rval );
#endif
        // print verbosely about the problem setting
        {
            moab::Range rintxverts, rintxelems;
            rval = mbCore->get_entities_by_dimension ( runCtx->meshsets[0], 0, rintxverts ); MB_CHK_ERR ( rval );
            rval = mbCore->get_entities_by_dimension ( runCtx->meshsets[0], 2, rintxelems ); MB_CHK_ERR ( rval );
            rval = fix_degenerate_quads ( mbCore, runCtx->meshsets[0] ); MB_CHK_ERR ( rval );
            rval = positive_orientation ( mbCore, runCtx->meshsets[0], radius_src ); MB_CHK_ERR ( rval );
            if ( !proc_id ) outputFormatter.printf ( 0,  "The source set contains %lu vertices and %lu elements \n", rintxverts.size(), rintxelems.size() );

            moab::Range bintxverts, bintxelems;
            rval = mbCore->get_entities_by_dimension ( runCtx->meshsets[1], 0, bintxverts ); MB_CHK_ERR ( rval );
            rval = mbCore->get_entities_by_dimension ( runCtx->meshsets[1], 2, bintxelems ); MB_CHK_ERR ( rval );
            rval = fix_degenerate_quads ( mbCore, runCtx->meshsets[1] ); MB_CHK_ERR ( rval );
            rval = positive_orientation ( mbCore, runCtx->meshsets[1], radius_dest ); MB_CHK_ERR ( rval );
            if ( !proc_id ) outputFormatter.printf ( 0,  "The target set contains %lu vertices and %lu elements \n", bintxverts.size(), bintxelems.size() );
        }

        // First compute the covering set such that the target elements are fully covered by the lcoal source grid
        runCtx->timer_push ( "construct covering set for intersection" );
        rval = remapper.ConstructCoveringSet ( epsrel, 1.0, 1.0, 0.1, runCtx->rrmGrids ); MB_CHK_ERR ( rval );
        runCtx->timer_pop();

        // Compute intersections with MOAB with either the Kd-tree or the advancing front algorithm 
        runCtx->timer_push ( "setup and compute mesh intersections" );
        rval = remapper.ComputeOverlapMesh ( runCtx->kdtreeSearch, false ); MB_CHK_ERR ( rval );
        runCtx->timer_pop();

        // print some diagnostic checks to see if the overlap grid resolved the input meshes correctly
        {
            double local_areas[4], global_areas[4]; // Array for Initial area, and through Method 1 and Method 2
            // local_areas[0] = area_on_sphere_lHuiller ( mbCore, runCtx->meshsets[1], radius_src );
            local_areas[0] = area_on_sphere_lHuiller ( mbCore, runCtx->meshsets[0], radius_src );
            local_areas[1] = area_on_sphere_lHuiller ( mbCore, runCtx->meshsets[1], radius_dest );
            local_areas[2] = area_on_sphere_lHuiller ( mbCore, runCtx->meshsets[2], radius_src );
            local_areas[3] = area_on_sphere ( mbCore, runCtx->meshsets[2], radius_src );

#ifdef MOAB_HAVE_MPI
            MPI_Allreduce ( &local_areas[0], &global_areas[0], 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
            global_areas[0] = local_areas[0];
            global_areas[1] = local_areas[1];
            global_areas[2] = local_areas[2];
            global_areas[3] = local_areas[3];
#endif
            if ( !proc_id )
            {
                outputFormatter.printf ( 0, "initial area: source = %12.14f, target = %12.14f, overlap with l'Huiller: %12.14f\n", global_areas[0], global_areas[1], global_areas[2] );
                // outputFormatter.printf ( 0, " area with l'Huiller: %12.14f with Girard: %12.14f\n", global_areas[2], global_areas[3] );
                // outputFormatter.printf ( 0, " relative difference areas = %12.10e\n", fabs ( global_areas[2] - global_areas[3] ) / global_areas[2] );
                outputFormatter.printf ( 0, " relative error w.r.t source = %12.14e, and target = %12.14e\n", fabs ( global_areas[0] - global_areas[2] ) / global_areas[0], fabs ( global_areas[1] - global_areas[2] ) / global_areas[1] );
            }
        }

        if ( runCtx->intxFilename.size() )
        {
            moab::EntityHandle writableOverlapSet;
            rval = mbCore->create_meshset ( moab::MESHSET_SET, writableOverlapSet ); MB_CHK_SET_ERR ( rval, "Can't create new set" );
            moab::EntityHandle meshOverlapSet = remapper.GetMeshSet ( moab::Remapper::OverlapMesh );
            moab::Range ovEnts;
            rval = mbCore->get_entities_by_dimension ( meshOverlapSet, 2, ovEnts ); MB_CHK_SET_ERR ( rval, "Can't create new set" );
            rval = mbCore->get_entities_by_dimension ( meshOverlapSet, 0, ovEnts ); MB_CHK_SET_ERR ( rval, "Can't create new set" );

#ifdef MOAB_HAVE_MPI
            // Do not remove ghosted entities if we still haven't computed weights
            // Remove ghosted entities from overlap set before writing the new mesh set to file
            if (nprocs > 1)
            {
                moab::Range ghostedEnts;
                rval = remapper.GetOverlapAugmentedEntities(ghostedEnts); MB_CHK_ERR ( rval );
                ovEnts = moab::subtract(ovEnts, ghostedEnts);
            }
#endif
            rval = mbCore->add_entities(writableOverlapSet, ovEnts);MB_CHK_SET_ERR(rval, "Deleting ghosted entities failed");

            size_t lastindex = runCtx->intxFilename.find_last_of(".");
            sstr.str("");
            sstr << runCtx->intxFilename.substr(0, lastindex) << ".h5m";
            if(!runCtx->proc_id) std::cout << "Writing out the MOAB intersection mesh file to " << sstr.str() << std::endl;

            // Write out our computed intersection file
            rval = mbCore->write_file ( sstr.str().c_str(), NULL, writeOptions, &writableOverlapSet, 1 ); MB_CHK_ERR ( rval );
        }

        if ( runCtx->computeWeights )
        {
            runCtx->meshes[2] = remapper.GetMesh ( moab::Remapper::OverlapMesh );
            if(!runCtx->proc_id) std::cout << std::endl;

            runCtx->timer_push ( "setup computation of weights" );
            // Call to generate the remapping weights with the tempest meshes
            moab::TempestOnlineMap* weightMap = new moab::TempestOnlineMap ( &remapper );
            runCtx->timer_pop();

            runCtx->timer_push ( "compute weights with TempestRemap" );

            rval = weightMap->GenerateRemappingWeights ( runCtx->disc_methods[0], runCtx->disc_methods[1],         // std::string strInputType, std::string strOutputType,
                                                   runCtx->disc_orders[0],  runCtx->disc_orders[1],                // int nPin=4, int nPout=4,
                                                   runCtx->fNoBubble, runCtx->ensureMonotonicity,                  // bool fNoBubble=true, int fMonotoneTypeID=0,
                                                   runCtx->fVolumetric, runCtx->fNoConservation, runCtx->fNoCheck, // bool fVolumetric=false, bool fNoConservation=false, bool fNoCheck=false,
                                                   runCtx->doftag_names[0], runCtx->doftag_names[1],               // std::string source_tag_name, std::string target_tag_name,
                                                   "", //"",                                                       // std::string strVariables="",
                                                   "", "",                                                         // std::string strInputData="", std::string strOutputData="",
                                                   "", true,                                                       // std::string strNColName="", bool fOutputDouble=true,
                                                   "", false, 0.0,                                                 // std::string strPreserveVariables="", bool fPreserveAll=false, double dFillValueOverride=0.0,
                                                   runCtx->fInputConcave, runCtx->fOutputConcave                   // bool fInputConcave = false, bool fOutputConcave = false
                                                 );MB_CHK_ERR ( rval );
            runCtx->timer_pop();

            // Invoke the CheckMap routine on the TempestRemap serial interface directly, if running on a single process
            if (nprocs == 1) {
                const double dNormalTolerance = 1.0E-8;
                const double dStrictTolerance = 1.0E-12;
                weightMap->CheckMap(!runCtx->fNoCheck, !runCtx->fNoCheck, !runCtx->fNoCheck && (runCtx->ensureMonotonicity), dNormalTolerance, dStrictTolerance);
            }

            if ( runCtx->outFilename.size() )
            {
                size_t lastindex = runCtx->outFilename.find_last_of(".");
                sstr.str("");
                sstr << runCtx->outFilename.substr(0, lastindex) << ".h5m";
                // Write the map file to disk in parallel
                rval = weightMap->WriteParallelMap(sstr.str().c_str());MB_CHK_ERR ( rval );

                // Write out the metadata information for the map file
                if (proc_id == 0) {
                    sstr.str("");
                    sstr << runCtx->outFilename.substr(0, lastindex) << ".meta";

                    std::ofstream metafile(sstr.str());
                    metafile << "Generator = MOAB-TempestRemap (mbtempest) Offline Regridding Weight Generator" << std::endl;
                    metafile << "domain_a = " << runCtx->inFilenames[0] << std::endl;
                    metafile << "domain_b = " << runCtx->inFilenames[1] << std::endl;
                    metafile << "grid_file_src = " << runCtx->inFilenames[0] << std::endl;
                    metafile << "grid_file_dst = " << runCtx->inFilenames[1] << std::endl;
                    metafile << "grid_file_ovr = " << (runCtx->intxFilename.size() ? runCtx->intxFilename : "outOverlap.h5m") << std::endl;
                    metafile << "mono_type = " << runCtx->ensureMonotonicity << std::endl;
                    metafile << "np_src = " << runCtx->disc_orders[0] << std::endl;
                    metafile << "np_dst = " << runCtx->disc_orders[1] << std::endl;
                    metafile << "type_src = " << runCtx->disc_methods[0] << std::endl;
                    metafile << "type_dst = " << runCtx->disc_methods[1] << std::endl;
                    metafile << "bubble = " << (runCtx->fNoBubble ? "false" : "true") << std::endl;
                    metafile << "concave_src = " << (runCtx->fInputConcave ? "true" : "false") << std::endl;
                    metafile << "concave_dst = " << (runCtx->fOutputConcave ? "true" : "false") << std::endl;
                    metafile << "version = " << "MOAB v5.1.0+" << std::endl;
                    metafile.close();
                }
            }

            if ( runCtx->verifyWeights )
            {
                // Let us pick a sampling test function for solution evaluation
                moab::TempestOnlineMap::sample_function testFunction = &sample_stationary_vortex; // &sample_slow_harmonic;

                runCtx->timer_push ( "describe a solution on source grid" );
                moab::Tag srcAnalyticalFunction;
                rval = weightMap->DefineAnalyticalSolution ( srcAnalyticalFunction, "AnalyticalSolnSrcExact", 
                                                             moab::Remapper::SourceMesh, 
                                                             testFunction);MB_CHK_ERR ( rval );
                runCtx->timer_pop();
                // rval = mbCore->write_file ( "srcWithSolnTag.h5m", NULL, writeOptions, &runCtx->meshsets[0], 1 ); MB_CHK_ERR ( rval );

                runCtx->timer_push ( "describe a solution on target grid" );
                moab::Tag tgtAnalyticalFunction;
                moab::Tag tgtProjectedFunction;
                rval = weightMap->DefineAnalyticalSolution ( tgtAnalyticalFunction, "AnalyticalSolnTgtExact", 
                                                             moab::Remapper::TargetMesh, 
                                                             testFunction,
                                                             &tgtProjectedFunction,
                                                             "ProjectedSolnTgt");MB_CHK_ERR ( rval );
                // rval = mbCore->write_file ( "tgtWithSolnTag.h5m", NULL, writeOptions, &runCtx->meshsets[1], 1 ); MB_CHK_ERR ( rval );
                runCtx->timer_pop();

                runCtx->timer_push ( "compute solution projection on target grid" );
                rval = weightMap->ApplyWeights(srcAnalyticalFunction, tgtProjectedFunction);MB_CHK_ERR ( rval );
                runCtx->timer_pop();
                rval = mbCore->write_file ( "tgtWithSolnTag2.h5m", NULL, writeOptions, &runCtx->meshsets[1], 1 ); MB_CHK_ERR ( rval );

                runCtx->timer_push ( "compute error metrics against analytical solution on target grid" );
                std::map<std::string, double> errMetrics;
                rval = weightMap->ComputeMetrics(moab::Remapper::TargetMesh, tgtAnalyticalFunction, tgtProjectedFunction, errMetrics, true);MB_CHK_ERR ( rval );
                runCtx->timer_pop();
            }

            delete weightMap;
        }
    }

    // Clean up
    remapper.clear();
    delete runCtx;
    delete mbCore;

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    exit ( 0 );
}


static moab::ErrorCode CreateTempestMesh ( ToolContext& ctx, moab::TempestRemapper& remapper, Mesh* tempest_mesh )
{
    moab::ErrorCode rval = moab::MB_SUCCESS;
    int err;
    moab::DebugOutput& outputFormatter = ctx.outputFormatter;

    if ( !ctx.proc_id ) { outputFormatter.printf ( 0,  "Creating TempestRemap Mesh object ...\n" ); }

    if ( ctx.meshType == moab::TempestRemapper::OVERLAP_FILES )
    {
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        err = GenerateOverlapMesh ( ctx.inFilenames[0], ctx.inFilenames[1], *tempest_mesh, ctx.outFilename, "NetCDF4", "exact", true );

        if ( err ) { rval = moab::MB_FAILURE; }
        else
        {
            ctx.meshes.push_back ( tempest_mesh );
        }
    }
    else if ( ctx.meshType == moab::TempestRemapper::OVERLAP_MEMORY )
    {
        // Load the meshes and validate
        ctx.meshsets.resize ( 3 );
        ctx.meshes.resize ( 3 );
        ctx.meshsets[0] = remapper.GetMeshSet ( moab::Remapper::SourceMesh );
        ctx.meshsets[1] = remapper.GetMeshSet ( moab::Remapper::TargetMesh );
        ctx.meshsets[2] = remapper.GetMeshSet ( moab::Remapper::OverlapMesh );

        // First the source
        rval = remapper.LoadMesh ( moab::Remapper::SourceMesh, ctx.inFilenames[0], moab::TempestRemapper::DEFAULT ); MB_CHK_ERR ( rval );
        ctx.meshes[0] = remapper.GetMesh ( moab::Remapper::SourceMesh );

        // Next the target
        rval = remapper.LoadMesh ( moab::Remapper::TargetMesh, ctx.inFilenames[1], moab::TempestRemapper::DEFAULT ); MB_CHK_ERR ( rval );
        ctx.meshes[1] = remapper.GetMesh ( moab::Remapper::TargetMesh );

        // Now let us construct the overlap mesh, by calling TempestRemap interface directly
        // For the overlap method, choose between: "fuzzy", "exact" or "mixed"
        err = GenerateOverlapWithMeshes ( *ctx.meshes[0], *ctx.meshes[1], *tempest_mesh, "" /*ctx.outFilename*/, "NetCDF4", "exact", false );

        if ( err ) { rval = moab::MB_FAILURE; }
        else
        {
            remapper.SetMesh ( moab::Remapper::OverlapMesh, tempest_mesh );
            ctx.meshes[2] = remapper.GetMesh ( moab::Remapper::OverlapMesh );
            // ctx.meshes.push_back(*tempest_mesh);
        }
    }
    else if ( ctx.meshType == moab::TempestRemapper::OVERLAP_MOAB )
    {
        ctx.meshsets.resize ( 3 );
        ctx.meshes.resize ( 3 );
        ctx.meshsets[0] = remapper.GetMeshSet ( moab::Remapper::SourceMesh );
        ctx.meshsets[1] = remapper.GetMeshSet ( moab::Remapper::TargetMesh );
        ctx.meshsets[2] = remapper.GetMeshSet ( moab::Remapper::OverlapMesh );

        const double radius_src = 1.0 /*2.0*acos(-1.0)*/;
        const double radius_dest = 1.0 /*2.0*acos(-1.0)*/;

        const char* additional_read_opts = (ctx.n_procs > 1 ? "NO_SET_CONTAINING_PARENTS;" : "");
        // Load the source mesh and validate
        rval = remapper.LoadNativeMesh ( ctx.inFilenames[0], ctx.meshsets[0], additional_read_opts ); MB_CHK_ERR ( rval );
        // Rescale the radius of both to compute the intersection
        rval = ScaleToRadius(ctx.mbcore, ctx.meshsets[0], radius_src);MB_CHK_ERR ( rval );
        rval = remapper.ConvertMeshToTempest ( moab::Remapper::SourceMesh ); MB_CHK_ERR ( rval );
        ctx.meshes[0] = remapper.GetMesh ( moab::Remapper::SourceMesh );

        // Load the target mesh and validate
        rval = remapper.LoadNativeMesh ( ctx.inFilenames[1], ctx.meshsets[1], additional_read_opts ); MB_CHK_ERR ( rval );
        rval = ScaleToRadius(ctx.mbcore, ctx.meshsets[1], radius_dest);MB_CHK_ERR ( rval );
        rval = remapper.ConvertMeshToTempest ( moab::Remapper::TargetMesh ); MB_CHK_ERR ( rval );
        ctx.meshes[1] = remapper.GetMesh ( moab::Remapper::TargetMesh );
    }
    else if ( ctx.meshType == moab::TempestRemapper::ICO )
    {
        err = GenerateICOMesh ( *tempest_mesh, ctx.blockSize, ctx.computeDual, ctx.outFilename, "NetCDF4" );

        if ( err ) { rval = moab::MB_FAILURE; }
        else
        {
            ctx.meshes.push_back ( tempest_mesh );
        }
    }
    else if ( ctx.meshType == moab::TempestRemapper::RLL )
    {
        err = GenerateRLLMesh ( *tempest_mesh,                    // Mesh& meshOut,
                                ctx.blockSize * 2, ctx.blockSize, // int nLongitudes, int nLatitudes,
                                0.0, 360.0,                       // double dLonBegin, double dLonEnd,
                                -90.0, 90.0,                      // double dLatBegin, double dLatEnd,
                                false, false, false,              // bool fGlobalCap, bool fFlipLatLon, bool fForceGlobal,
                                "" /*ctx.inFilename*/, "", "",    // std::string strInputFile, std::string strInputFileLonName, std::string strInputFileLatName,
                                ctx.outFilename, "NetCDF4",       // std::string strOutputFile, std::string strOutputFormat
                                true                              // bool fVerbose
                              );

        if ( err ) { rval = moab::MB_FAILURE; }
        else
        {
            ctx.meshes.push_back ( tempest_mesh );
        }
    }
    else   // default
    {
        err = GenerateCSMesh ( *tempest_mesh, ctx.blockSize, ctx.outFilename, "NetCDF4" );

        if ( err ) { rval = moab::MB_FAILURE; }
        else
        {
            ctx.meshes.push_back ( tempest_mesh );
        }
    }

    if ( ctx.meshType != moab::TempestRemapper::OVERLAP_MOAB && !tempest_mesh )
    {
        std::cout << "Tempest Mesh is not a complete object; Quitting...";
        exit ( -1 );
    }

    return rval;
}


///////////////////////////////////////////////
//         Test functions

double sample_slow_harmonic ( double dLon, double dLat )
{
  return (2.0 + cos(dLat) * cos(dLat) * cos(2.0 * dLon));
}

double sample_fast_harmonic ( double dLon, double dLat )
{
  return (2.0 + pow(sin(2.0 * dLat), 16.0) * cos(16.0 * dLon));
	//return (2.0 + pow(cos(2.0 * dLat), 16.0) * cos(16.0 * dLon));
}

double sample_constant ( double /*dLon*/, double /*dLat*/ )
{
  return 1.0;
}

double sample_stationary_vortex ( double dLon, double dLat )
{
  const double dLon0 = 0.0;
  const double dLat0 = 0.6;
  const double dR0 = 3.0;
  const double dD = 5.0;
  const double dT = 6.0;

  ///		Find the rotated longitude and latitude of a point on a sphere
  ///		with pole at (dLonC, dLatC).
  {
    double dSinC = sin(dLat0);
    double dCosC = cos(dLat0);
    double dCosT = cos(dLat);
    double dSinT = sin(dLat);

    double dTrm  = dCosT * cos(dLon - dLon0);
    double dX = dSinC * dTrm - dCosC * dSinT;
    double dY = dCosT * sin(dLon - dLon0);
    double dZ = dSinC * dSinT + dCosC * dTrm;

    dLon = atan2(dY, dX);
    if (dLon < 0.0) {
        dLon += 2.0 * M_PI;
    }
    dLat = asin(dZ);
  }

  double dRho = dR0 * cos(dLat);
  double dVt = 3.0 * sqrt(3.0) / 2.0 / cosh(dRho) / cosh(dRho) * tanh(dRho);

  double dOmega;
  if (dRho == 0.0) {
    dOmega = 0.0;
  } else {
    dOmega = dVt / dRho;
  }

  return (1.0 - tanh(dRho / dD * sin(dLon - dOmega * dT)));
}

///////////////////////////////////////////////

