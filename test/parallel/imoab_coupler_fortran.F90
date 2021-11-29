!      This program shows how to do an imoab_coupler type test in fortran 90
!      the atm and ocean files are read from repo, meshes are migrated to coupler pes,
!      the intx is caried on, weight generation,; tag migration and projection,
!      migrate back to ocean pes, and compare against baseline1.txt
!      cannot change source/target files, or the layout, just test the imoab
!       calls, that need to give the same results as C=++ equivalent test (imoab_coupler)
!          between atm SE and ocn FV

SUBROUTINE errorout(ierr, message)
      integer ierr
      character*(*) message
      if (ierr.ne.0) then
        print *, message
        call exit (1)
      end if
     return
     end
!
#include "moab/MOABConfig.h"

#ifndef MOAB_MESH_DIR
#error Specify MOAB_MESH_DIR path
#endif

program imoab_coupler_fortran

    use iso_c_binding
    use iMOAB

#include "mpif.h"
#include "moab/MOABConfig.h"
    !implicit none
    integer :: m ! for number of arguments ; if less than 1, exit
    integer :: global_comm
    integer :: rankInGlobalComm, numProcesses
    integer :: ierr
    integer :: my_id, num_procs
    integer :: jgroup   ! group for global comm
    character(:), allocatable :: atmFileName
    character(:), allocatable :: ocnFileName
    character(:), allocatable :: baselineFileName
    character(:), allocatable :: readopts, fileWriteOptions
    character :: appname*128
    character(:), allocatable :: weights_identifier1
    character(:), allocatable :: disc_methods1, disc_methods2, dof_tag_names1, dof_tag_names2
    integer :: disc_orders1, disc_orders2
    character(:), allocatable :: atmocn_map_file_name, intx_from_file_identifier

    ! all groups and comms are on the same number of processes, for simplicity
    integer :: atmGroup, atmComm, ocnGroup, ocnComm, cplGroup, cplComm
    integer :: atmCouComm, ocnCouComm

    integer :: cmpatm, cplatm, cmpocn, cplocn, atmocnid ! global IDs
    integer :: cmpAtmPID, cplAtmPID ! iMOAB app ids
    integer :: cmpOcnPID, cplOcnPID ! iMOAB app ids
    integer :: cplAtmOcnPID ! intx pid
    integer :: nghlay, partScheme, context_id
    integer :: fNoBubble, fMonotoneTypeID, fVolumetric, fNoConserve, fValidate

    integer, dimension(2) ::  tagIndex
    integer, dimension (2) :: tagTypes!  { DENSE_DOUBLE, DENSE_DOUBLE }
    integer :: atmCompNDoFs ! = disc_orders[0] * disc_orders[0],
    integer :: ocnCompNDoFs !  = 1 /*FV*/
    character(:), allocatable :: bottomTempField , bottomTempProjectedField, bottomUVelField
    character(:), allocatable :: bottomUVelProjectedField, bottomVVelField, bottomVVelProjectedField
    integer, dimension(3) ::  nverts, nelem, nblocks, nsbc, ndbc
    DOUBLE PRECISION, allocatable :: vals(:) ! to set the double values to 0
    integer :: i ! for loops
    integer :: storLeng, eetype ! for tags defs
    character(:), allocatable :: concat_fieldname, concat_fieldnameT, outputFileOcn
    integer :: tagIndexIn2 ! not really needed

    cmpatm = 5
    cplatm = 6
    cmpocn = 17
    cplocn = 18
    atmocnid = 618

    call mpi_init(ierr)
    call errorout(ierr,'mpi_init')
    call mpi_comm_rank (MPI_COMM_WORLD, my_id, ierr)
    call errorout(ierr, 'fail to get MPI rank')


    call mpi_comm_size (MPI_COMM_WORLD, num_procs, ierr)
    call errorout(ierr, 'fail to get MPI size')
    call mpi_comm_dup(MPI_COMM_WORLD, global_comm, ierr)
    call errorout(ierr, 'fail to get global comm duplicate')

    call MPI_Comm_group( global_comm, jgroup, ierr );  !  all processes in jgroup
    call errorout(ierr, 'fail to get joint group')
    ! readopts( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS" )
    ! readoptsLnd( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION" )
    atmFileName = MOAB_MESH_DIR//'unittest/wholeATM_T.h5m'//C_NULL_CHAR
    ocnFileName = MOAB_MESH_DIR//'unittest/recMeshOcn.h5m'//C_NULL_CHAR
    baselineFileName = MOAB_MESH_DIR//'unittest/baseline1.txt'//C_NULL_CHAR

    ! all comms span the whole world, for simplicity
    atmComm = MPI_COMM_NULL
    ocnComm = MPI_COMM_NULL
    cplComm = MPI_COMM_NULL
    atmCouComm = MPI_COMM_NULL
    ocnCouComm = MPI_COMM_NULL
    call mpi_comm_dup(global_comm, atmComm, ierr)
    call mpi_comm_dup(global_comm, ocnComm, ierr)
    call mpi_comm_dup(global_comm, cplComm, ierr)
    call mpi_comm_dup(global_comm, atmCouComm, ierr)
    call mpi_comm_dup(global_comm, ocnCouComm, ierr)

    ! all groups of interest are easy breezy
    call MPI_Comm_group( atmComm, atmGroup, ierr )
    call MPI_Comm_group( ocnComm, ocnGroup, ierr )
    call MPI_Comm_group( cplComm, cplGroup, ierr )

    readopts ='PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS'//C_NULL_CHAR
    nghlay = 0 ! no ghost layers
#ifdef MOAB_HAVE_ZOLTAN
    partScheme = 2  ! RCB with zoltan
#else
    partScheme = 0  ! Trivial partitioner
#endif

    if (my_id.eq.0) then
        print *,' number of tasks: ', num_procs
        print *,' Atm file: ', atmFileName
        print *,' Ocn file: ', ocnFileName
        print *,' baseline file: ', baselineFileName
        print *,' using partitioner: ', partScheme
    endif

    ierr = iMOAB_Initialize()
    appname = 'ATM'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, atmComm, cmpatm, cmpatmPid)
    appname = 'ATMX'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, cplComm, cplatm, cplatmPid)
    appname = 'OCN'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, ocnComm, cmpocn, cmpocnPid)
    appname = 'OCNX'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, cplComm, cplocn, cplocnPid)

    appname = 'ATMOCN'//C_NULL_CHAR
    ierr = iMOAB_RegisterApplication(appname, cplComm, atmocnid, cplAtmOcnPID)

    ! read atm and migrate
    if (atmComm .NE. MPI_COMM_NULL) then
       ierr = iMOAB_LoadMesh( cmpatmPid, atmFileName, readopts, nghlay )
       call errorout(ierr, 'fail to load atm')
       ierr = iMOAB_SendMesh( cmpatmPid, atmCouComm, cplGroup, cplatm, partScheme )
       call errorout(ierr, 'fail to send atm')
    endif
    if( cplComm .NE. MPI_COMM_NULL ) then
        ierr = iMOAB_ReceiveMesh( cplatmPid, atmCouComm, atmGroup, cmpatm ) !
        call errorout(ierr, 'fail to receive atm')
    endif

    ! we can now free the sender buffers
    if (atmComm .NE. MPI_COMM_NULL) then
        context_id = cplatm
        ierr           = iMOAB_FreeSenderBuffers( cmpatmPid, context_id )
        call errorout(ierr, 'fail to free atm buffers')
    endif

    ! read ocn and migrate
    if (ocnComm .NE. MPI_COMM_NULL) then
       ierr = iMOAB_LoadMesh( cmpocnPid, ocnFileName, readopts, nghlay )
       call errorout(ierr, 'fail to load ocn')
       ierr = iMOAB_SendMesh( cmpocnPid, ocnCouComm, cplGroup, cplocn, partScheme )
       call errorout(ierr, 'fail to send ocn')
    endif
    if( cplComm .NE. MPI_COMM_NULL ) then
        ierr = iMOAB_ReceiveMesh( cplocnPid, ocnCouComm, ocnGroup, cmpocn ) !
        call errorout(ierr, 'fail to receive ocn')
    endif

    ! we can now free the sender buffers
    if (ocnComm .NE. MPI_COMM_NULL) then
        context_id = cplocn
        ierr           = iMOAB_FreeSenderBuffers( cmpocnPid, context_id )
        call errorout(ierr, 'fail to free ocn buffers')
    endif

    if( cplComm .NE. MPI_COMM_NULL ) then
        ierr = iMOAB_ComputeMeshIntersectionOnSphere( cplAtmPID, cplOcnPID,  cplAtmOcnPID )
        ! coverage mesh was computed here, for cplAtmPID, atm on coupler pes
        ! basically, atm was redistributed according to target (ocean) partition, to "cover" the
        !ocean partitions check if intx valid, write some h5m intx file
        call errorout(ierr, 'cannot compute intersection')
    endif

    if( atmCouComm .NE. MPI_COMM_NULL ) then
        ! the new graph will be for sending data from atm comp to coverage mesh.
        ! it involves initial atm app; cmpAtmPID; also migrate atm mesh on coupler pes, cplAtmPID
        ! results are in cplAtmOcnPID, intx mesh; remapper also has some info about coverage mesh
        ! after this, the sending of tags from atm pes to coupler pes will use the new par comm
        ! graph, that has more precise info about what to send for ocean cover ; every time, we
        ! will use the element global id, which should uniquely identify the element
        ierr = iMOAB_CoverageGraph( atmCouComm, cmpAtmPID, cplAtmPID, cplAtmOcnPID, cmpatm, cplatm, cplocn ) ! it happens over joint communicator
        call errorout(ierr, 'cannot recompute direct coverage graph for ocean')
    endif

    weights_identifier1='scalar'//C_NULL_CHAR
    disc_methods1='cgll'//C_NULL_CHAR
    disc_methods2='fv'//C_NULL_CHAR
    disc_orders1=4
    disc_orders2=1
    dof_tag_names1='GLOBAL_DOFS'//C_NULL_CHAR
    dof_tag_names2='GLOBAL_ID'//C_NULL_CHAR
    ! fMonotoneTypeID = 0, fVolumetric = 0, fValidate = 1, fNoConserve = 0, fNoBubble = 1
    fNoBubble = 1
    fMonotoneTypeID = 0
    fVolumetric = 0
    fNoConserve = 0
    fValidate = 1

    if( cplComm .NE. MPI_COMM_NULL ) then

        ierr = iMOAB_ComputeScalarProjectionWeights(  &
            cplAtmOcnPID, weights_identifier1, disc_methods1, disc_orders1,  &
            disc_methods2,  disc_orders2, fNoBubble, fMonotoneTypeID, fVolumetric, fNoConserve, &
            fValidate, dof_tag_names1, dof_tag_names2 )
        call errorout(ierr, 'cannot compute scalar projection weights')

#ifdef MOAB_HAVE_NETCDF
        atmocn_map_file_name = 'atm_ocn_map_f.nc'//C_NULL_CHAR
        ierr = iMOAB_WriteMappingWeightsToFile( cplAtmOcnPID, weights_identifier1, atmocn_map_file_name)
        call errorout(ierr, 'failed to write map file to disk')
        intx_from_file_identifier = 'map-from-file'//C_NULL_CHAR
        ierr = iMOAB_LoadMappingWeightsFromFile_Old( cplAtmOcnPID, intx_from_file_identifier,  &
                                                     atmocn_map_file_name)
        call errorout(ierr, 'failed to load map file from disk')
#endif
    endif

    ! start copy
    tagTypes(1) =  1 ! somehow, DENSE_DOUBLE give 0, while it should be 1; maybe moab::DENSE_DOUBLE ?
    tagTypes(2) =  1 ! ! DENSE_DOUBLE
    atmCompNDoFs = disc_orders1 * disc_orders1
    ocnCompNDoFs = 1 ! /*FV*/

    bottomTempField          = 'a2oTbot'//C_NULL_CHAR
    bottomTempProjectedField = 'a2oTbot_proj'//C_NULL_CHAR
    bottomUVelField          = 'a2oUbot'//C_NULL_CHAR
    bottomUVelProjectedField = 'a2oUbot_proj'//C_NULL_CHAR
    bottomVVelField          = 'a2oVbot'//C_NULL_CHAR
    bottomVVelProjectedField = 'a2oVbot_proj'//C_NULL_CHAR

    if( cplComm .NE. MPI_COMM_NULL ) then
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomTempField, tagTypes(1), atmCompNDoFs, tagIndex(1))
        call errorout(ierr, 'failed to define the field tag a2oTbot')
        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomTempProjectedField, tagTypes(2), ocnCompNDoFs, tagIndex(2))
        call errorout(ierr, 'failed to define the field tag a2oTbot_proj')
        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomUVelField,  tagTypes(1), atmCompNDoFs, tagIndex(1))
        call errorout(ierr, 'failed to define the field tag a2oUbot')
        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomUVelProjectedField, tagTypes(2), ocnCompNDoFs, tagIndex(2))
        call errorout(ierr, 'failed to define the field tag a2oUbot_proj')

        ierr = iMOAB_DefineTagStorage( cplAtmPID, bottomVVelField,  tagTypes(1), atmCompNDoFs, tagIndex(1))
        call errorout(ierr, 'failed to define the field tag a2oUbot')

        ierr = iMOAB_DefineTagStorage( cplOcnPID, bottomVVelProjectedField,tagTypes(2), ocnCompNDoFs, tagIndex(2))
        call errorout(ierr, 'failed to define the field tag a2oUbot_proj')
    endif

    ! make the tag 0, to check we are actually sending needed data
    if( cplAtmPID .ge. 0 ) then

        !  Each process in the communicator will have access to a local mesh instance, which
        !  will contain the original cells in the local partition and ghost entities. Number of
        !  vertices, primary cells, visible blocks, number of sidesets and nodesets boundary
        !  conditions will be returned in numProcesses 3 arrays, for local, ghost and total
        !  numbers.

        ierr = iMOAB_GetMeshInfo( cplAtmPID, nverts, nelem, nblocks, nsbc, ndbc )
        call errorout(ierr, 'failed to get num primary elems')
        storLeng = nelem(3)*atmCompNDoFs
        allocate (vals(storLeng ) )
        eetype = 1 ! double type

        do i = 1, storLeng
            vals(:) = 0.
        enddo

        ! set the tag values to 0.0
        ierr = iMOAB_SetDoubleTagStorage( cplAtmPID, bottomTempField, storLeng, eetype, vals)
        call errorout(ierr, 'cannot make tag nul')
        ierr = iMOAB_SetDoubleTagStorage( cplAtmPID, bottomUVelField, storLeng, eetype, vals)
        call errorout(ierr, 'cannot make tag nul')
        ierr = iMOAB_SetDoubleTagStorage( cplAtmPID, bottomVVelField, storLeng, eetype, vals)
        call errorout(ierr, 'cannot make tag nul')

    endif

    ! Define the field variables to project
    concat_fieldname  = 'a2oTbot;a2oUbot;a2oVbot;'//C_NULL_CHAR
    concat_fieldnameT = 'a2oTbot_proj;a2oUbot_proj;a2oVbot_proj;'//C_NULL_CHAR

    if( atmComm .NE. MPI_COMM_NULL ) then

        ! As always, use nonblocking sends
        ! this is for projection to ocean:
        ierr = iMOAB_SendElementTag( cmpAtmPID, concat_fieldname, atmCouComm, cplocn)
        call errorout(ierr, 'cannot send tag values')

    endif

    if( cplComm .NE. MPI_COMM_NULL ) then
        !// receive on atm on coupler pes, that was redistributed according to coverage
        ierr = iMOAB_ReceiveElementTag( cplAtmPID, concat_fieldname, atmCouComm, cplocn)
        call errorout(ierr, 'cannot receive tag values')
    endif

        ! we can now free the sender buffers
    if( atmComm .NE. MPI_COMM_NULL ) then

        ierr = iMOAB_FreeSenderBuffers( cmpAtmPID, cplocn ) !context is for ocean
        call errorout(ierr, 'cannot free buffers used to resend atm tag towards the coverage mesh')

    endif
    if( cplComm .ne. MPI_COMM_NULL )  then

        outputFileOcn= "AtmOnCplF.h5m"//C_NULL_CHAR
        fileWriteOptions = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
        ierr                 = iMOAB_WriteMesh( cplAtmPID, outputFileOcn, fileWriteOptions)
        call errorout(ierr, 'could not write AtmOnCpl.h5m to disk')

    endif
    if( cplComm .ne. MPI_COMM_NULL ) then

        ! We have the remapping weights now. Let us apply the weights onto the tag we defined
        ! on the source mesh and get the projection on the target mesh
        ierr = iMOAB_ApplyScalarProjectionWeights( cplAtmOcnPID, weights_identifier1, concat_fieldname, &
                        concat_fieldnameT)
        call errorout(ierr, 'failed to compute projection weight application')

        outputFileOcn= "OcnOnCplF.h5m"//C_NULL_CHAR
        fileWriteOptions = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
        ierr                 = iMOAB_WriteMesh( cplOcnPID, outputFileOcn, fileWriteOptions)
        call errorout(ierr, 'could not write OcnOnCpl.h5m to disk')

    endif
    ! send the projected tag back to ocean pes, with send/receive tag
    ! first makje sure the tags are defined, otherwise they cannot be received
    if( ocnComm .ne. MPI_COMM_NULL ) then

        ierr = iMOAB_DefineTagStorage( cmpOcnPID, bottomTempProjectedField, tagTypes(2), ocnCompNDoFs, tagIndexIn2 )
        call errorout(ierr, 'failed to define the field tag for receiving back the tag a2oTbot_proj on ocn pes')

        ierr = iMOAB_DefineTagStorage( cmpOcnPID, bottomUVelProjectedField, tagTypes(2), ocnCompNDoFs, tagIndexIn2 )
        call errorout(ierr, 'failed to define the field tag for receiving back the tag a2oUbot_proj on ocn pes')

        ierr = iMOAB_DefineTagStorage( cmpOcnPID, bottomVVelProjectedField, tagTypes(2), ocnCompNDoFs, tagIndexIn2 )
        call errorout(ierr, 'failed to define the field tag for receiving back the tag a2oVbot_proj on ocn pes')

    endif

    !  send the tag to ocean pes, from ocean mesh on coupler pes
    !  from couComm, using common joint comm ocn_coupler
    !  as always, use nonblocking sends
    !  original graph (context is -1_
    if( cplComm .ne. MPI_COMM_NULL ) then
        context_id = cmpocn
        ierr = iMOAB_SendElementTag( cplOcnPID, concat_fieldnameT, ocnCouComm, context_id)
        call errorout(ierr, 'cannot send tag values back to ocean pes')
    endif

    if( ocnComm .ne. MPI_COMM_NULL ) then
       context_id = cplocn
       ierr = iMOAB_ReceiveElementTag( cmpOcnPID, concat_fieldnameT, ocnCouComm, context_id)
       call errorout(ierr, 'cannot receive tag values from ocean mesh on coupler pes')
    endif

    if( cplComm .ne. MPI_COMM_NULL ) then
        context_id = cmpocn
        ierr = iMOAB_FreeSenderBuffers( cplOcnPID, context_id )
        call errorout(ierr, 'cannot free sender buffers on coupler')
    endif

    if( ocnComm .ne. MPI_COMM_NULL )  then

        outputFileOcn= "OcnWithProjF.h5m"//C_NULL_CHAR
        fileWriteOptions = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
        if (my_id.eq.0) then
            print *,' Writing ocean mesh file with projected solution to disk: ', outputFileOcn
        endif
        ierr = iMOAB_WriteMesh( cmpOcnPID, outputFileOcn, fileWriteOptions)
        call errorout(ierr, 'could not write OcnWithProjF.h5m to disk')

    endif

    ! free up resources
    ierr = iMOAB_DeregisterApplication(cplAtmOcnPID)
    call errorout(ierr, 'could not de-register OCN component')

    ierr = iMOAB_DeregisterApplication(cplOcnPID)
    call errorout(ierr, 'could not de-register OCN component')
    ierr = iMOAB_DeregisterApplication(cmpOcnPID)
    call errorout(ierr, 'could not de-register OCN component')
    ierr = iMOAB_DeregisterApplication(cplAtmPID)
    call errorout(ierr, 'could not de-register OCN component')
    ierr = iMOAB_DeregisterApplication(cmpAtmPID)
    call errorout(ierr, 'could not de-register OCN component')

    ! Free all MPI datastructures
    call MPI_Comm_free(atmComm, ierr)
    call MPI_Group_free(atmGroup, ierr)
    call MPI_Comm_free(ocnComm, ierr)
    call MPI_Group_free(ocnGroup, ierr)
    call MPI_Comm_free(cplComm, ierr)
    call MPI_Group_free(cplGroup, ierr)

    call MPI_Comm_free(atmCouComm, ierr)
    call MPI_Comm_free(ocnCouComm, ierr)
    call MPI_Comm_free(global_comm, ierr)
    call MPI_Group_free(jgroup, ierr)
    call MPI_Finalize(ierr)

end program imoab_coupler_fortran
