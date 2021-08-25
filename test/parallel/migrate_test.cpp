/*
 * migrate_test will contain tests for migrating meshes in parallel environments, with iMOAB api
 * these  methods are tested also in the example MigrateMesh.F90, with variable
 * numbers of processes; migrate_test is launched usually on 2 processes, and it tests
 * various cases
 * a mesh is read on senders tasks, sent to receivers tasks, and then written out for verification
 * It depends on hdf5 parallel for reading and writing in parallel
 */

#include "moab/ParallelComm.hpp"
#include "moab/Core.hpp"
#include "moab_mpi.h"
#include "moab/iMOAB.h"
#include "TestUtil.hpp"

#define RUN_TEST_ARG2( A, B ) run_test( &( A ), #A, B )

using namespace moab;

#define CHECKRC( rc, message )     \
    if( 0 != ( rc ) )              \
    {                              \
        printf( "%s\n", message ); \
        return MB_FAILURE;         \
    }

int is_any_proc_error( int is_my_error )
{
    int result = 0;
    int err    = MPI_Allreduce( &is_my_error, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
    return err || result;
}

int run_test( ErrorCode ( *func )( const char* ), const char* func_name, const char* file_name )
{
    ErrorCode result = ( *func )( file_name );
    int is_err       = is_any_proc_error( ( MB_SUCCESS != result ) );
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 )
    {
        if( is_err )
            std::cout << func_name << " : FAILED!!" << std::endl;
        else
            std::cout << func_name << " : success" << std::endl;
    }

    return is_err;
}

ErrorCode migrate_1_1( const char* filename );
ErrorCode migrate_1_2( const char* filename );
ErrorCode migrate_2_1( const char* filename );
ErrorCode migrate_2_2( const char* filename );
ErrorCode migrate_4_2( const char* filename );
ErrorCode migrate_2_4( const char* filename );
ErrorCode migrate_4_3( const char* filename );
ErrorCode migrate_overlap( const char* filename );

// some global variables, used by all tests
int rank, size, ierr;

int compid1, compid2;  // component ids are unique over all pes, and established in advance;
int nghlay;            // number of ghost layers for loading the file
int groupTasks[4];     // at most 4 tasks
int startG1, startG2, endG1, endG2;

MPI_Comm jcomm;  // will be a copy of the global
MPI_Group jgroup;

int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    MPI_Comm_dup( MPI_COMM_WORLD, &jcomm );
    MPI_Comm_group( jcomm, &jgroup );

    std::string filename;
    filename = TestDir + "/field1.h5m";
    if( argc > 1 )
    {
        filename = argv[1];
    }
    int num_errors = 0;
    num_errors += RUN_TEST_ARG2( migrate_1_1, filename.c_str() );
    num_errors += RUN_TEST_ARG2( migrate_1_2, filename.c_str() );
    num_errors += RUN_TEST_ARG2( migrate_2_1, filename.c_str() );
    num_errors += RUN_TEST_ARG2( migrate_2_2, filename.c_str() );
    if( size >= 4 )
    {
        num_errors += RUN_TEST_ARG2( migrate_4_2, filename.c_str() );
        num_errors += RUN_TEST_ARG2( migrate_2_4, filename.c_str() );
        num_errors += RUN_TEST_ARG2( migrate_4_3, filename.c_str() );
        num_errors += RUN_TEST_ARG2( migrate_overlap, filename.c_str() );
    }
    if( rank == 0 )
    {
        if( !num_errors )
            std::cout << "All tests passed" << std::endl;
        else
            std::cout << num_errors << " TESTS FAILED!" << std::endl;
    }

    MPI_Group_free( &jgroup );
    MPI_Comm_free( &jcomm );
    MPI_Finalize();
    return num_errors;
}

ErrorCode migrate( const char* filename, const char* outfile )
{
    // first create MPI groups

    std::string filen( filename );
    MPI_Group group1, group2;
    for( int i = startG1; i <= endG1; i++ )
        groupTasks[i - startG1] = i;

    ierr = MPI_Group_incl( jgroup, endG1 - startG1 + 1, groupTasks, &group1 );
    CHECKRC( ierr, "can't create group1" )

    for( int i = startG2; i <= endG2; i++ )
        groupTasks[i - startG2] = i;

    ierr = MPI_Group_incl( jgroup, endG2 - startG2 + 1, groupTasks, &group2 );
    CHECKRC( ierr, "can't create group2" )

    // create 2 communicators, one for each group
    int tagcomm1 = 1, tagcomm2 = 2;
    MPI_Comm comm1, comm2;
    ierr = MPI_Comm_create_group( jcomm, group1, tagcomm1, &comm1 );
    CHECKRC( ierr, "can't create comm1" )

    ierr = MPI_Comm_create_group( jcomm, group2, tagcomm2, &comm2 );
    CHECKRC( ierr, "can't create comm2" )

    ierr = iMOAB_Initialize( 0, 0 );  // not really needed anything from argc, argv, yet; maybe we should
    CHECKRC( ierr, "can't initialize iMOAB" )

    // give some dummy values to component ids, just to differentiate between them
    // the par comm graph is unique between components
    compid1        = 4;
    compid2        = 7;
    int context_id = -1;  // default context; will be now set to compid1 or compid2

    int appID1;
    iMOAB_AppID pid1 = &appID1;
    int appID2;
    iMOAB_AppID pid2 = &appID2;

    if( comm1 != MPI_COMM_NULL )
    {
        ierr = iMOAB_RegisterApplication( "APP1", &comm1, &compid1, pid1 );
        CHECKRC( ierr, "can't register app1 " )
    }
    if( comm2 != MPI_COMM_NULL )
    {
        ierr = iMOAB_RegisterApplication( "APP2", &comm2, &compid2, pid2 );
        CHECKRC( ierr, "can't register app2 " )
    }

    int method = 0;  // trivial partition for sending
    if( comm1 != MPI_COMM_NULL )
    {

        std::string readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" );

        nghlay = 0;

        ierr = iMOAB_LoadMesh( pid1, filen.c_str(), readopts.c_str(), &nghlay, filen.length(),
                               strlen( readopts.c_str() ) );
        CHECKRC( ierr, "can't load mesh " )
        ierr = iMOAB_SendMesh( pid1, &jcomm, &group2, &compid2, &method );  // send to component 2
        CHECKRC( ierr, "cannot send elements" )
    }

    if( comm2 != MPI_COMM_NULL )
    {
        ierr = iMOAB_ReceiveMesh( pid2, &jcomm, &group1, &compid1 );  // receive from component 1
        CHECKRC( ierr, "cannot receive elements" )
        std::string wopts;
        wopts = "PARALLEL=WRITE_PART;";
        ierr =
            iMOAB_WriteMesh( pid2, (char*)outfile, (char*)wopts.c_str(), strlen( outfile ), strlen( wopts.c_str() ) );
        CHECKRC( ierr, "cannot write received mesh" )
    }

    MPI_Barrier( jcomm );

    // we can now free the sender buffers
    if( comm1 != MPI_COMM_NULL ) ierr = iMOAB_FreeSenderBuffers( pid1, &context_id );

    // exchange tag, from component to component
    // one is receiving, one is sending the tag; the one that is sending needs to have communicator
    // not null
    int size_tag  = 1;  // a double dense tag, on elements
    int tagType   = DENSE_DOUBLE;
    int tagIndex2 = 0, tagIndex1 = 0;  // these will be tag indices on each app pid

    std::string fileAfterTagMigr( outfile );  // has h5m
    int sizen = fileAfterTagMigr.length();
    fileAfterTagMigr.erase( sizen - 4, 4 );  // erase extension .h5m
    fileAfterTagMigr = fileAfterTagMigr + "_tag.h5m";

    // now send a tag from component 2, towards component 1
    if( comm2 != MPI_COMM_NULL )
    {
        ierr =
            iMOAB_DefineTagStorage( pid2, "element_field", &tagType, &size_tag, &tagIndex2, strlen( "element_field" ) );
        CHECKRC( ierr, "failed to get tag element_field " );
        // this tag is already existing in the file

        // first, send from compid2 to compid1, from comm2, using common joint comm
        // as always, use nonblocking sends
        // contex_id should be now compid1
        context_id = compid1;
        ierr       = iMOAB_SendElementTag( pid2, "element_field", &jcomm, &context_id, strlen( "element_field" ) );
        CHECKRC( ierr, "cannot send tag values" )
    }
    // receive on component 1
    if( comm1 != MPI_COMM_NULL )
    {
        ierr =
            iMOAB_DefineTagStorage( pid1, "element_field", &tagType, &size_tag, &tagIndex1, strlen( "element_field" ) );
        CHECKRC( ierr, "failed to get tag DFIELD " );
        context_id = compid2;
        ierr       = iMOAB_ReceiveElementTag( pid1, "element_field", &jcomm, &context_id, strlen( "element_field" ) );
        CHECKRC( ierr, "cannot send tag values" )
        std::string wopts;
        wopts = "PARALLEL=WRITE_PART;";
        ierr  = iMOAB_WriteMesh( pid1, (char*)fileAfterTagMigr.c_str(), (char*)wopts.c_str(), fileAfterTagMigr.length(),
                                strlen( wopts.c_str() ) );
        CHECKRC( ierr, "cannot write received mesh" )
    }

    MPI_Barrier( jcomm );

    // we can now free the sender buffers
    if( comm2 != MPI_COMM_NULL ) ierr = iMOAB_FreeSenderBuffers( pid2, &context_id );

    if( comm2 != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( pid2 );
        CHECKRC( ierr, "cannot deregister app 2 receiver" )
    }

    if( comm1 != MPI_COMM_NULL )
    {
        ierr = iMOAB_DeregisterApplication( pid1 );
        CHECKRC( ierr, "cannot deregister app 1 sender" )
    }

    ierr = iMOAB_Finalize();
    CHECKRC( ierr, "did not finalize iMOAB" )

    if( MPI_COMM_NULL != comm1 ) MPI_Comm_free( &comm1 );
    if( MPI_COMM_NULL != comm2 ) MPI_Comm_free( &comm2 );

    MPI_Group_free( &group1 );
    MPI_Group_free( &group2 );
    return MB_SUCCESS;
}
// migrate from task 0 to task 1, non overlapping
ErrorCode migrate_1_1( const char* filename )
{
    startG1 = endG1 = 0;
    startG2 = endG2 = 1;
    return migrate( filename, "migrate11.h5m" );
}
// migrate from task 0 to 2 tasks (0 and 1)
ErrorCode migrate_1_2( const char* filename )
{
    startG1 = endG1 = startG2 = 0;
    endG2                     = 1;
    return migrate( filename, "migrate12.h5m" );
}

// migrate from 2 tasks (0, 1) to 1 task (0)
ErrorCode migrate_2_1( const char* filename )
{
    startG1 = endG2 = startG2 = 0;
    endG1                     = 1;
    return migrate( filename, "migrate21.h5m" );
}

// migrate from 2 tasks to 2 tasks (overkill)
ErrorCode migrate_2_2( const char* filename )
{
    startG1 = startG2 = 0;
    endG1 = endG2 = 1;
    return migrate( filename, "migrate22.h5m" );
}
// migrate from 4 tasks to 2 tasks
ErrorCode migrate_4_2( const char* filename )
{
    startG1 = startG2 = 0;
    endG2             = 1;
    endG1             = 3;
    return migrate( filename, "migrate42.h5m" );
}

ErrorCode migrate_2_4( const char* filename )
{
    startG1 = startG2 = 0;
    endG2             = 3;
    endG1             = 1;
    return migrate( filename, "migrate24.h5m" );
}

ErrorCode migrate_4_3( const char* filename )
{
    startG1 = startG2 = 0;
    endG2             = 2;
    endG1             = 3;
    return migrate( filename, "migrate43.h5m" );
}

ErrorCode migrate_overlap( const char* filename )
{
    startG1 = 0;
    startG2 = 1;
    endG1   = 1;
    endG2   = 2;
    return migrate( filename, "migrate_over.h5m" );
}
