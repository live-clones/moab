/*
 * imoab_coupler_utils.hpp
 *
 *  Created on: Aug. 22, 2020
 *  \brief will contain utility methods for refactoring imoab*coupler tests, to avoid repetitive tasks
 *  \ even migrate tests can use some of these utilities
 *  1) create_comm_group(int start, int end, int tag, MPI_Group& group, MPI_Comm& comm)
 *
 */

#ifndef TEST_PARALLEL_IMOAB_COUPLER_UTILS_HPP_
#define TEST_PARALLEL_IMOAB_COUPLER_UTILS_HPP_

#define CHECKIERR( rc, message )                       \
    if( 0 != ( rc ) )                                  \
    {                                                  \
        printf( "%s. ErrorCode = %d\n", message, rc ); \
        return 1;                                      \
    }
#define PUSH_TIMER( localcomm, operation )               \
    {                                                    \
        MPI_Barrier(localcomm);                          \
        timer_ops = timer.time_since_birth();            \
        opName    = operation;                           \
    }
#define POP_TIMER( localcomm, localrank )                                                         \
    {                                                                                             \
        double locElapsed = timer.time_since_birth() - timer_ops, sum = 0, maxElapsed = 0;        \
        MPI_Reduce( &locElapsed, &maxElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, localcomm );             \
        MPI_Reduce( &locElapsed, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, localcomm );                    \
        int procsize;                                                                             \
        MPI_Comm_size(localcomm, &procsize);                                                      \
        if( !( localrank ) )                                                                      \
            std::cout << "[LOG] Time taken to " << opName.c_str() << ": max = " << maxElapsed     \
                      << ", avg = " << sum/procsize << "\n";                   \
        opName.clear();                                                                           \
    }

/*
 *  \brief create an MPI group and an MPI communicator for the group, in the global communicator
 */
int create_group_and_comm( int start, int end, MPI_Group worldGroup, MPI_Group* group, MPI_Comm* comm )
{
    std::vector< int > groupTasks;
    groupTasks.resize( end - start + 1, 0 );
    for( int i = start; i <= end; i++ )
        groupTasks[i - start] = i;

    int ierr = MPI_Group_incl( worldGroup, end - start + 1, &groupTasks[0], group );
    CHECKIERR( ierr, "Cannot create group" )

    ierr = MPI_Comm_create( MPI_COMM_WORLD, *group, comm );
    CHECKIERR( ierr, "Cannot create comm" )

    return 0;
}

int create_joint_comm_group( MPI_Group agroup, MPI_Group bgroup, MPI_Group* abgroup, MPI_Comm* abcomm )
{
    int ierr = MPI_Group_union( agroup, bgroup, abgroup );
    CHECKIERR( ierr, "Cannot create joint union group" )

    ierr = MPI_Comm_create( MPI_COMM_WORLD, *abgroup, abcomm );
    CHECKIERR( ierr, "Cannot create joint communicator from union group" )

    return 0;
}

int setup_component_coupler_meshes( iMOAB_AppID cmpId,
                                    int cmpTag,
                                    iMOAB_AppID cplCmpId,
                                    int cmpcouTag,
                                    MPI_Comm* cmpcomm,
                                    MPI_Group* cmpPEGroup,
                                    MPI_Comm* coucomm,
                                    MPI_Group* cplPEGroup,
                                    MPI_Comm* cmpcoucomm,
                                    std::string& filename,
                                    std::string& readopts,
                                    int nghlay,
                                    int repartitioner_scheme )
{
    int ierr = 0;
    if( *cmpcomm != MPI_COMM_NULL )
    {
        // load mesh on component side
        ierr =
            iMOAB_LoadMesh( cmpId, filename.c_str(), readopts.c_str(), &nghlay );
        CHECKIERR( ierr, "Cannot load component mesh" )
    }
    if ( *cmpcoucomm != MPI_COMM_NULL &&  5 == repartitioner_scheme )
    {
        // need to send the zoltan buffer from coupler root towards the component root
        // it assumes there is a buffer saved somewhere on the root of the coupler, because scheme 4 was called in advance
        int is_fort = 0;
        ierr = iMOAB_RetrieveZBuffer( cmpPEGroup, cplPEGroup, cmpcoucomm, &is_fort);
        CHECKIERR( ierr, "Cannot retrieve Zoltan buffer" )
    }

    if( *cmpcomm != MPI_COMM_NULL )
    {
        // then send mesh to coupler pes
        ierr = iMOAB_SendMesh( cmpId, cmpcoucomm, cplPEGroup, &cmpcouTag,
                               &repartitioner_scheme );  // send to  coupler pes
        CHECKIERR( ierr, "cannot send elements" )
    }
    // now, receive mesh, on coupler communicator; first mesh 1, atm
    if( *coucomm != MPI_COMM_NULL )
    {

        ierr = iMOAB_ReceiveMesh( cplCmpId, cmpcoucomm, cmpPEGroup,
                                  &cmpTag );  // receive from component
        CHECKIERR( ierr, "cannot receive elements on coupler app" )
    }

    // we can now free the sender buffers
    if( *cmpcomm != MPI_COMM_NULL )
    {
        int context_id = cmpcouTag;
        ierr           = iMOAB_FreeSenderBuffers( cmpId, &context_id );
        CHECKIERR( ierr, "cannot free buffers used to send atm mesh" )
    }
    return 0;
}

#endif /* TEST_PARALLEL_IMOAB_COUPLER_UTILS_HPP_ */
