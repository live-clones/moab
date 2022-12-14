#include "moab/MOABConfig.h"

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

#include "moab/iMOAB.h"
#include <string.h>

#define CHECKRC( rc, message )             \
    if( 0 != ( rc ) )                      \
    {                                      \
        printf( "Error: %s.\n", message ); \
        return 1;                          \
    }

int main( int argc, char* argv[] )
{
    /* Declarations */
    ErrCode rc;
    int nprocs = 1, rank = 0;
    char filen[1024] = "";
    int i, j, k;

    int irank = 0, num_tags_to_sync = 1;
    int* ranks;

    int appID, appDupID;
    iMOAB_AppID pid;
    iMOAB_AppID pidDup;
    int compid = 11, compdupid = 12;

    int num_global_vertices = 0, num_global_elements = 0, num_dimension = 0, num_parts = 0;
    int nverts[3], nelem[3], nblocks[3], nsbc[3], ndbc[3];

    iMOAB_GlobalID* vGlobalID;
    int* vranks;
    double* coords;
    int size_coords;

#ifdef MOAB_HAVE_MPI
    const char* read_opts   = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";
    int num_ghost_layers[1] = { 1 };
#else
    const char* read_opts   = "";
    int num_ghost_layers[1] = { 0 };
#endif

    iMOAB_GlobalID* gbIDs;
    int tagIndex[2];
    int entTypes[2]    = { 0, 1 }; /* first is on vertex, second is on elements; */
    int tagTypes[2]    = { DENSE_INTEGER, DENSE_DOUBLE };
    int num_components = 1;

    iMOAB_GlobalID *element_global_IDs, *block_IDs;
    int vertices_per_element, num_elements_in_block;
    int conn[27], nv, eindex;

    int* int_tag_vals;
    double* double_tag_vals;

    int local_index           = 0;  /* test element with local index 0 */
    int num_adjacent_elements = 10; /* we can have maximum 6 actually */
    iMOAB_LocalID adjacent_element_IDs[10], *element_connectivity, *local_element_ID;
    int size_conn, *element_ownership;
    iMOAB_GlobalID* global_element_ID;

    iMOAB_LocalID *vertBC_ID, *surfBC_ID;
    int *vertBC_value, *ref_surf, *bc_value;

    char *outputFile, *writeOptions;

#ifdef MOAB_HAVE_MPI
    MPI_Comm comm;
    MPI_Init( &argc, &argv );

    comm = MPI_COMM_WORLD;
    MPI_Comm_size( comm, &nprocs );
    MPI_Comm_rank( comm, &rank );
#endif

#ifdef MOAB_HAVE_HDF5
    strcpy( filen, MOAB_MESH_DIR );
    strcat( filen, "unittest/io/p8ex1.h5m" );
#endif

    if( argc > 1 ) strcpy( filen, argv[1] );

    if( !strlen( filen ) )
    {
        printf( "Invalid filename specified at command line.\n" );
        return 1;
    }
    /*
     * MOAB needs to be initialized; A MOAB instance will be created, and will be used by each
     * application in this framework. There is no parallel context yet.
     */
    rc = iMOAB_Initialize( argc, argv );
    CHECKRC( rc, "failed to initialize MOAB" );

    /*
     * Header information is cheap to retrieve from an hdf5 file; it can be called by each process
     * or can be called by one master task and broadcasted. The hdf5 file has to be in MOAB native
     * format, and has to be partitioned in advance. There should be at least as many partitions in
     * the file as number of tasks in the communicator.
     */
    rc = iMOAB_ReadHeaderInfo( filen, &num_global_vertices, &num_global_elements, &num_dimension, &num_parts );
    CHECKRC( rc, "failed to read header info" );

    if( 0 == rank )
    {
        printf( "file %s has %d vertices, %d elements, %d parts in partition\n", filen, num_global_vertices,
                num_global_elements, num_parts );
    }

    pid    = &appID;
    pidDup = &appDupID;

    /*
     * Each application has to be registered once. A mesh set and a parallel communicator will be
     * associated with each application. A unique application id will be returned, and will be used
     * for all future mesh operations/queries.
     */
    rc = iMOAB_RegisterApplication( "MBAPP",
#ifdef MOAB_HAVE_MPI
                                    &comm,
#endif
                                    &compid, pid );
    CHECKRC( rc, "failed to register application" );

    rc = iMOAB_RegisterApplication( "MBDUPAPP",
#ifdef MOAB_HAVE_MPI
                                    &comm,
#endif
                                    &compdupid, pidDup );
    CHECKRC( rc, "failed to register duplicate application" );

    /*
     * Loading the mesh is a parallel IO operation. Ghost layers can be exchanged too, and default
     * MOAB sets are augmented with ghost elements. By convention, blocks correspond to MATERIAL_SET
     * sets, side sets with NEUMANN_SET sets, node sets with DIRICHLET_SET sets. Each element and
     * vertex entity should have a GLOBAL ID tag in the file, which will be available for visible
     * entities
     */
    rc = iMOAB_LoadMesh( pid, filen, read_opts, num_ghost_layers );
    CHECKRC( rc, "failed to load mesh" );

    rc = iMOAB_LoadMesh( pidDup, filen, read_opts, num_ghost_layers );
    CHECKRC( rc, "failed to load mesh" );

    rc = iMOAB_SetGlobalInfo( pid, &num_global_vertices, &num_global_elements );
    CHECKRC( rc, "failed to set global info" );

    /*
     * Each process in the communicator will have access to a local mesh instance, which will
     * contain the original cells in the local partition and ghost entities. Number of vertices,
     * primary cells, visible blocks, number of sidesets and nodesets boundary conditions will be
     * returned in size 3 arrays, for local, ghost and total numbers.
     */
    rc = iMOAB_GetMeshInfo( pid, nverts, nelem, nblocks, nsbc, ndbc );
    CHECKRC( rc, "failed to get mesh info" );

    vGlobalID = (iMOAB_GlobalID*)malloc( nverts[2] * sizeof( iMOAB_GlobalID ) );
    /*
     * All visible vertices have a unique global ID. All arrays are allocated by client.
     * Currently, global ID is a 32 bit integer, total number of visible vertices was returned by
     * iMOAB_GetMeshInfo
     */
    rc = iMOAB_GetVertexID( pid, &nverts[2], vGlobalID );
    CHECKRC( rc, "failed to get vertex id info" );

    vranks = (int*)malloc( nverts[2] * sizeof( int ) );
    /*
     * In MOAB terminology, a vertex can be owned, shared or ghost. An owned vertex will be
     * owned by the current process. A shared vertex appears at the interface between partitions, it
     * can be owned by another process, while a ghost vertex is owned by another process always. All
     * vertex arrays have a natural local order, in which the ghost vertices are at the end of the
     * arrays. The local index of the vertex will vary from 0 to number of visible vertices (in this
     * example nverts[2]-1)
     */
    rc = iMOAB_GetVertexOwnership( pid, &nverts[2], vranks );
    CHECKRC( rc, "failed to get vertex ranks" );

    coords      = (double*)malloc( 3 * nverts[2] * sizeof( double ) );
    size_coords = 3 * nverts[2];
    /*
     * Coordinates are returned interleaved, for all visible vertices. Size is dimension times
     * number of visible vertices
     */
    rc = iMOAB_GetVisibleVerticesCoordinates( pid, &size_coords, coords );
    CHECKRC( rc, "failed to get coordinates" );

    /*
     * Block IDs correspond to MATERIAL_SET tag values in hdf5 file
     * They also correspond to block IDs in a Cubit model.
     * Local visible blocks might contain ghost cells, owned by other
     * processes.
     * 
     * Inside MOAB, a block will correspond to a meshset with a MATERIAL_SET tag,
     * with value the block ID.
     * A MATERIAL_SET tag in moab if SPARSE and type INTEGER.
     */
    gbIDs = (iMOAB_GlobalID*)malloc( nblocks[2] * sizeof( iMOAB_GlobalID ) );
    rc    = iMOAB_GetBlockID( pid, &nblocks[2], gbIDs );
    CHECKRC( rc, "failed to get block info" );

    /*
     * The 2 tags used in this example exist in the file, already.
     * If a user needs a new tag, it can be defined with the same call as this one
     * this method, iMOAB_DefineTagStorage, will return a local index for the tag.
     * The name of the tag is case sensitive.
     * This method is collective.
     */
    rc = iMOAB_DefineTagStorage( pid, "INTFIELD", &tagTypes[0], &num_components, &tagIndex[0] );
    CHECKRC( rc, "failed to get tag INTFIELD " );
    rc = iMOAB_DefineTagStorage( pid, "DFIELD", &tagTypes[1], &num_components, &tagIndex[1] );
    CHECKRC( rc, "failed to get tag DFIELD " );

    /*
     * Synchronize one of the tags only, just to see what happens;
     * Tags are not exchanged for ghost entities by default. The user has to explicitly
     * synchronize them with a call like this one. The tags to synchronize are referenced by the
     * index returned at their "definition" with iMOAB_DefineTagStorage method.
     */
    rc = iMOAB_SynchronizeTags( pid, &num_tags_to_sync, &tagIndex[0], &tagTypes[0] );
    CHECKRC( rc, "failed to sync tag INTFIELD " );

    for( irank = 0; irank < nprocs; irank++ )
    {
        if( irank == rank )
        {
            /* printf some of the block info */
            printf( "on rank %d, there are \n"
                    "  %3d visible vertices of which  %3d local  %3d ghost \n"
                    "  %3d visible elements of which  %3d owned  %3d ghost \n"
                    "  %3d visible blocks\n"
                    "  %3d visible neumann BCs\n"
                    "  %3d visible dirichlet BCs\n",
                    rank, nverts[2], nverts[0], nverts[1], nelem[2], nelem[0], nelem[1], nblocks[2], nsbc[2], ndbc[2] );

            /* print some of the vertex id infos */
            printf( "on rank %d vertex info:\n", rank );
            for( i = 0; i < nverts[2]; i++ )
                printf( " vertex local id: %3d, rank ID:%d  global ID: %3d  coords: %g, %g, %g\n", i, vranks[i],
                        vGlobalID[i], coords[3 * i], coords[3 * i + 1], coords[3 * i + 2] );

            element_global_IDs = (iMOAB_GlobalID*)malloc( nelem[2] * sizeof( iMOAB_GlobalID ) );
            block_IDs          = (iMOAB_GlobalID*)malloc( nelem[2] * sizeof( iMOAB_GlobalID ) );
            ranks              = (int*)malloc( nelem[2] * sizeof( int ) );
            /*
             * Visible elements info is available for all visible elements, with this call. Similar
             * information can be retrieved by looping over visible blocks
             */
            rc = iMOAB_GetVisibleElementsInfo( pid, &nelem[2], element_global_IDs, ranks, block_IDs );
            CHECKRC( rc, "failed to get all elem info" );
            for( i = 0; i < nelem[2]; i++ )
                printf( " element local id: %3d,  global ID: %3d  rank:%d  block ID: %2d \n", i, element_global_IDs[i],
                        ranks[i], block_IDs[i] );
            free( element_global_IDs );
            free( ranks );
            free( block_IDs );

            /* Get first element connectivity */
            nv     = 27;
            eindex = 0;
            /*
             * Element connectivity can be retrieved for a single element. Local vertex indices are returned.
             */
            rc = iMOAB_GetElementConnectivity( pid, &eindex, &nv, conn );
            CHECKRC( rc, "failed to get first element connectivity" );
            printf( " conn for first element: \n" );
            for( i = 0; i < nv; i++ )
                printf( " %3d", conn[i] );
            printf( "\n" );

            /*
             * Test neighbors:
             * Neighbor elements that share a face with current element can be determined with this
             * call Elements are identified by their local index
             */
            rc = iMOAB_GetNeighborElements( pid, &local_index, &num_adjacent_elements, adjacent_element_IDs );
            CHECKRC( rc, "failed to get first element neighbors" );
            printf( "  neighbors for first element:\n" );
            for( i = 0; i < num_adjacent_elements; i++ )
            {
                printf( "  %4d", adjacent_element_IDs[i] );
            }
            printf( "\n" );

            for( i = 0; i < nblocks[2]; i++ )
            {
                printf( " block index: %3d, block ID: %3d \n", i, gbIDs[i] );
                /*
                 * Blocks should have the same type of primary elements. Number of elements in block
                 * and number of vertices per element are found with this call. The method refers
                 * only to the visible elements on the process. The block is identified with its
                 * block ID. Should it be identified with its local index in the list of visible
                 * blocks?
                 */
                rc = iMOAB_GetBlockInfo( pid, &gbIDs[i], &vertices_per_element, &num_elements_in_block );
                CHECKRC( rc, "failed to elem block info" );
                printf( "    has %4d elements with %d vertices per element\n", num_elements_in_block,
                        vertices_per_element );

                size_conn            = num_elements_in_block * vertices_per_element;
                element_connectivity = (iMOAB_LocalID*)malloc( sizeof( iMOAB_LocalID ) * size_conn );
                /*
                 * Connectivity for all elements in the block can be determined with this method.
                 * Vertices are identified by their local index (from 0 to number of visible
                 * vertices -1)
                 */
                rc = iMOAB_GetBlockElementConnectivities( pid, &gbIDs[i], &size_conn, element_connectivity );
                CHECKRC( rc, "failed to get block elem connectivity" );
                element_ownership = (int*)malloc( sizeof( int ) * num_elements_in_block );

                /*
                 * Each element is owned by a particular process. Ownership can be returned for all
                 * elements in a block with this call.
                 */
                rc = iMOAB_GetElementOwnership( pid, &gbIDs[i], &num_elements_in_block, element_ownership );
                CHECKRC( rc, "failed to get block elem ownership" );
                global_element_ID = (iMOAB_GlobalID*)malloc( sizeof( iMOAB_GlobalID ) * num_elements_in_block );
                local_element_ID  = (iMOAB_LocalID*)malloc( sizeof( iMOAB_LocalID ) * num_elements_in_block );

                /*
                 * Global element IDs are determined with this call. Local indices within the
                 * process are returned too.
                 */
                rc = iMOAB_GetElementID( pid, &gbIDs[i], &num_elements_in_block, global_element_ID, local_element_ID );
                CHECKRC( rc, "failed to get block elem IDs" );
                for( j = 0; j < num_elements_in_block; j++ )
                {
                    printf( "  elem %3d owned by %d gid: %4d lid: %4d  -- ", j, element_ownership[j],
                            global_element_ID[j], local_element_ID[j] );
                    for( k = 0; k < vertices_per_element; k++ )
                        printf( " %5d", element_connectivity[j * vertices_per_element + k] );
                    printf( "\n" );
                }
                free( global_element_ID );
                free( local_element_ID );
                free( element_connectivity );
                free( element_ownership );
            }
            /*
             * Query integer tag values on vertices.
             * this tag was synchronized (see iMOAB_SynchronizeTags above)
             * Being a dense tag, all visible vertices have a value.
             * As usual, the client allocates the returned array, which is filled with data with
             * this call (deep copy)
             */
            int_tag_vals = (int*)malloc( sizeof( int ) * nverts[2] ); /* for all visible vertices on the rank */
            rc           = iMOAB_GetIntTagStorage( pid, "INTFIELD", &nverts[2], &entTypes[0], int_tag_vals );
            CHECKRC( rc, "failed to get INTFIELD tag" );
            printf( "INTFIELD tag values:\n" );
            for( i = 0; i < nverts[2]; i++ )
            {
                printf( " %4d", int_tag_vals[i] );
                if( i % 20 == 19 ) printf( "\n" );
            }
            printf( "\n" );
            free( int_tag_vals );

            /*
             * Query double tag values on elements.
             * This tag was not synchronized, so ghost elements have a default value of 0.
             */
            double_tag_vals = (double*)malloc( sizeof( double ) * nelem[2] ); /* for all visible elements on the rank */
            rc              = iMOAB_GetDoubleTagStorage( pid, "DFIELD", &nelem[2], &entTypes[1], double_tag_vals );
            CHECKRC( rc, "failed to get DFIELD tag" );
            printf( "DFIELD tag values: (not exchanged) \n" );
            for( i = 0; i < nelem[2]; i++ )
            {
                printf( " %f", double_tag_vals[i] );
                if( i % 8 == 7 ) printf( "\n" );
            }
            printf( "\n" );
            free( double_tag_vals );

            /* query surface BCs */
            surfBC_ID = (iMOAB_LocalID*)malloc( sizeof( iMOAB_LocalID ) * nsbc[2] );
            ref_surf  = (int*)malloc( sizeof( int ) * nsbc[2] );
            bc_value  = (int*)malloc( sizeof( int ) * nsbc[2] );
            /*
             * Surface boundary condition information is returned for all visible Neumann conditions
             * the reference surfaces are indexed from 1 to 6 for hexahedrons, 1 to 4 to
             * tetrahedrons, in the MOAB canonical ordering convention
             */
            rc = iMOAB_GetPointerToSurfaceBC( pid, &nsbc[2], surfBC_ID, ref_surf, bc_value );
            CHECKRC( rc, "failed to get surf boundary conditions" );
            printf( " Surface boundary conditions:\n" );
            for( i = 0; i < nsbc[2]; i++ )
            {
                printf( "  elem_localID %4d  side:%d  BC:%2d\n", surfBC_ID[i], ref_surf[i], bc_value[i] );
            }
            free( surfBC_ID );
            free( ref_surf );
            free( bc_value );

            /* Query vertex BCs */
            vertBC_ID    = (iMOAB_LocalID*)malloc( sizeof( iMOAB_LocalID ) * ndbc[2] );
            vertBC_value = (int*)malloc( sizeof( int ) * ndbc[2] );
            rc           = iMOAB_GetPointerToVertexBC( pid, &ndbc[2], vertBC_ID, vertBC_value );
            CHECKRC( rc, "failed to get vertex boundary conditions" );
            printf( "  Vertex boundary conditions:\n" );
            for( i = 0; i < ndbc[2]; i++ )
            {
                printf( "   vertex %4d   BC:%2d\n", vertBC_ID[i], vertBC_value[i] );
            }
            free( vertBC_ID );
            free( vertBC_value );
        }
#ifdef MOAB_HAVE_MPI
        MPI_Barrier( comm ); /* to avoid printing problems */
#endif
    }

    /* free allocated data */
    free( coords );
    free( vGlobalID );
    free( vranks );
    outputFile = "fnew.h5m";
#ifdef MOAB_HAVE_MPI
    writeOptions = "PARALLEL=WRITE_PART";
#else
    writeOptions            = "";
#endif
    /*
     * The file can be written in parallel, and it will contain additional tags defined by the user
     * we may extend the method to write only desired tags to the file
     */
    rc = iMOAB_WriteMesh( pid, outputFile, writeOptions );
    CHECKRC( rc, "failed to write output file" );

    /*
     * Deregistering application will delete all mesh entities associated with the application and
     * will free allocated tag storage.
     */
    rc = iMOAB_DeregisterApplication( pid );
    CHECKRC( rc, "failed to de-register application" );

    rc = iMOAB_DeregisterApplication( pidDup );
    CHECKRC( rc, "failed to de-register application" );

    /*
     * This method will delete MOAB instance
     */
    rc = iMOAB_Finalize();
    CHECKRC( rc, "failed to finalize MOAB" );
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
