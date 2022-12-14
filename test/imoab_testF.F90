
SUBROUTINE errorout(ierr, message)
integer ierr
character*(*) message
if (ierr.ne.0) then
      print *, message
      call exit (1)
end if
return
end

program fdriver

use iso_c_binding
use iMOAB

#include "moab/MOABConfig.h"
#ifdef MOAB_HAVE_MPI
include 'mpif.h'
#endif

#ifndef MOAB_MESH_DIR
#error Specify MOAB_MESH_DIR path
#endif

      integer ierr, num_procs, my_id
      integer pid
      integer compid
      character :: appname*32
      character :: filename*1024
      character :: fname*1024
      character :: readopts*1024
      integer ngv, nge, ndim, nparts
      integer nghlay
      integer nverts(3), nelem(3), nblocks(3), nsbc(3), ndbc(3)
      !      large enough work arrays
      integer iwork(100000)
      double precision  dwork(100000)
      !      indices in work arrays for vertex ids, ranks, coordinates
      integer vID, vRA, vCO
      !
      !      indices for free memory (index in integer or double work arrays)
      integer ifree, dfree
      !      size of coordinates array
      integer  nCO
      integer  bID
      !     for some tags in the file
      character :: tagname1*128, tagname2*128
      integer      tagtype(2), enttype(2), num_co
      integer      tagindex(2)
      integer      ntsync

      integer      eRA, beID, eID
      !     vertices per element, number of elements in block
      integer  vpere, nebl, blockID
      !      iWORK(eCO) start for connectivity
      integer    sizeconn, eCO
      !      IWORK(egID) , IWORK(elID) starts for global el ID, local elem ID
      integer    egID, elID, eOWN
      integer    iTAG, dTAG

      !       indices for surface BC element, reference surf BC, value
      integer    isBC, irBC, ivBC

      character outfile*1024, wopts*1024
      my_id = 0
      fname = 'unittest/io/p8ex1.h5m'//C_NULL_CHAR

#ifdef MOAB_HAVE_MPI
      call MPI_INIT ( ierr )
      call errorout(ierr, 'fail to initialize MPI')

      !     find out MY process ID, and how many processes were started.

      call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
      call errorout(ierr, 'fail to get MPI rank')

      call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
      call errorout(ierr, 'fail to get MPI size')

      if (my_id .eq. 0) then
            print *, " I'm process ", my_id, " out of ", num_procs, " processes."
      end if
#endif

      ierr = iMOAB_Initialize()
      call errorout(ierr, 'fail to initialize iMOAB')

      appname = 'PROTEUS'//C_NULL_CHAR
      ! give a unique external id; here we just use 8?
      compid = 8
      !  number of ghost layers needed by application
      nghlay = 0


#ifdef MOAB_HAVE_MPI
      ierr = iMOAB_RegisterApplication(appname, MPI_COMM_WORLD, compid, pid)
#else
      ierr = iMOAB_RegisterApplication(appname, compid, pid)
#endif
      call errorout(ierr, 'fail to register application')

      filename = MOAB_MESH_DIR//fname
      ierr = iMOAB_ReadHeaderInfo ( filename, ngv, nge, ndim, nparts)
      call errorout(ierr, 'fail to read header info')

      if (0.eq.my_id) then
            print *, filename, ' has ', nparts, ' parts in partition', ngv, ' vertices ', nge, ' elements of dimension ', ndim
      endif
#ifdef MOAB_HAVE_MPI
      nghlay = 1
      readopts='PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;'// &
       'PARALLEL_RESOLVE_SHARED_ENTS'//C_NULL_CHAR
#else
      readopts=C_NULL_CHAR
#endif

      !  now let us load the mesh in parallel
      ierr = iMOAB_LoadMesh(pid, filename, readopts, nghlay)
      call errorout(ierr, 'fail to read file in parallel')

      !  number of vertices/elements/blocks/sidesets in the mesh
      ierr = iMOAB_GetMeshInfo(pid, nverts, nelem, nblocks, nsbc, ndbc)
      call errorout(ierr, 'fail to get mesh info')

      vID = 1
      ierr = iMOAB_GetVertexID(pid, nverts(3), IWORK(vID) )
      call errorout(ierr, 'failed to get vertex id info')

      vRA = vID + nverts(3)
      ierr = iMOAB_GetVertexOwnership(pid, nverts(3), IWORK(vRA) )
      call errorout(ierr, 'failed to get vertex owner ranks')

      ifree = vRA + nverts(3)

      !  double * coords = (double*) malloc(3*nverts[2]*sizeof(double));
      vCO = 1
      nCO = 3 * nverts(3)

      ierr = iMOAB_GetVisibleVerticesCoordinates(pid, nCO, DWORK(vCO))
      print *,  'nCO = ', nCO
      print *, 'ierr = ', ierr
      call errorout(ierr, 'failed to get coordinates')
      dfree = vCO + 3 * nverts(3)

      bID = ifree
      ierr = iMOAB_GetBlockID(pid, nblocks(3), IWORK(bID))
      call errorout(ierr, 'failed to get block info')
      ifree = ifree + nblocks(3)

      !      the 2 tags used in this example exist in the file, already
      ! first tag, INTFIELD is on vertices, integer
      ! second tag DFIELD is on elements, double
      tagtype(1)=0                                  !dense, int
      tagtype(2)=1                                  !dense, double
      enttype(1)=0                                  ! on verts
      enttype(2)=1                                  ! on elem
      num_co = 1
      tagname1 ='INTFIELD'//C_NULL_CHAR
      ierr = iMOAB_DefineTagStorage(pid, tagname1, tagtype(1), num_co, tagindex(1) )
      call errorout(ierr, 'failed to get tag INTFIELD')

      tagname2 ='DFIELD'//C_NULL_CHAR
      ierr = iMOAB_DefineTagStorage(pid, tagname2, tagtype(2), num_co, tagindex(2) )
      call errorout(ierr, 'failed to get tag DFIELD')

      ! synchronize one of the tags only, just to see what happens
      ! put in the sync array just first tag index (INTFIELD)
      ntsync =2

      ierr = iMOAB_SynchronizeTags(pid, ntsync, tagIndex, tagType )
      call errorout(ierr, 'failed to sync tag INTFIELD')

      !  start printing some information, retrieved from each task
      do irk=0, num_procs-1
            if (irk .eq. my_id) then

            ! printf some of the block info */
            print *, 'on rank ', my_id, ' there are '
            print *, nverts(3), ' visible vertices of which ',nverts(1), ' local ', nverts(2), ' ghost'
            print *,  nblocks(3), ' visible blocks'
            print *,  nsbc(3), ' visible Neumann BCs'
            print *,   ndbc(3), ' visible dirichlet BCs'

            print *, 'on rank ', my_id, ' vertex info:'
            do i=1,nverts(3)
            write(*, 100)  i, IWORK(vRA+i-1), IWORK( vID+i-1), DWORK(vCO+3*i-3), DWORK(vCO+3*i-2), DWORK(vCO+3*i-1)
100         FORMAT(' vertex local id ', I3, ' rank ID', I3, ' global ID:', I3, '  coords:',  3F11.3)
      enddo

      eID = ifree
      beID = eID + nelem(3)
      eRA = beID + nelem(3)
      ierr = iMOAB_GetVisibleElementsInfo(pid, nelem(3),IWORK(eID), IWORK(eRA), IWORK(beID) )
      call errorout(ierr, 'failed to get all elem info')
      ifree = eRA + nelem(3)
      do i=1, nelem(3)
          write(*, 101) IWORK(eID+i-1), IWORK(eRA+i-1), IWORK(beID+i-1)
101       FORMAT( ' global ID ', I5, ' rank: ', I3, ' block ID: ', I4)
      enddo

      do  i=1,nblocks(3)

            print *,' block index:', i, ' block ID ', IWORK(bID+i-1)
            blockID = IWORK(bID+i-1)
            ierr = iMOAB_GetBlockInfo(pid, blockID , vpere, nebl)
            call errorout(ierr, 'failed to elem block info')
            print *, '  has' , nebl, ' elements with ', vpere, 'verts'

            sizeconn = nebl * vpere
            eCO = ifree
            ierr = iMOAB_GetBlockElementConnectivities(pid, blockID, sizeconn, IWORK(eCO) )
            print *, ierr
            call errorout(ierr, 'failed to get block elem connectivity')

            ifree = ifree + sizeconn
            eOWN = ifree
            ierr = iMOAB_GetElementOwnership(pid, blockID, nebl, IWORK(eOWN))
            call errorout(ierr, 'failed to get block elem ownership')
            ifree = ifree+nebl

            egID = ifree
            elID = ifree + nebl
            ierr = iMOAB_GetElementID(pid, blockID, nebl, IWORK(egID), IWORK(elID)  )
            call errorout(ierr, 'failed to get block elem IDs')
            ifree = elID + nebl

            do j=1, nebl
                  write (*, 102) j,  IWORK(eOWN+j-1),IWORK(egID+j-1), IWORK(elID+j-1), (IWORK(eCO-1+(j-1)*vpere+k), k=1,vpere)
102               FORMAT(' elem ', I3, ' owned by', I3, ' gid:', I3, ' lid:', I3, ' : ', 10I5)
            enddo
      enddo

      ! query int tag values on vertices
      iTAG= ifree

      ierr = iMOAB_GetIntTagStorage(pid, tagname1, nverts(3), enttype(1), IWORK(iTAG) )
      call errorout(ierr, 'failed to get INTFIELD tag')
      ifree = iTAG + nverts(3)
      print * , 'INTFIELD tag values'
      write(*, 103) (IWORK(iTAG+k-1), k=1,nverts(3) )
103     FORMAT (10I8)

      dTAG = dfree
      ! query double tag values on elements

      ierr = iMOAB_GetDoubleTagStorage(pid, tagname2, nelem(3), entType(2), DWORK(dTAG) )
      call errorout(ierr, 'failed to get DFIELD tag')
      dfree = dTAG + nelem(3)
      print *, 'DFIELD tag values: (not exchanged) '
      write (*, 104) (DWORK(dTAG+k-1), k=1,nelem(3) )
104     FORMAT ( 10F7.2 )

      !  query surface BCs
      isBC = ifree
      irBC = isBC + nsbc(3)
      ivBC = irBC + nsbc(3)

      ierr = iMOAB_GetPointerToSurfaceBC(pid, nsbc(3), IWORK(isBC),IWORK(irBC), IWORK(ivBC))
      call errorout(ierr, 'failed to get surf boundary conditions')
      ifree = ivBC + nsbc(3)
      print * , 'Surface boundary conditions '
      write (*, 105) (IWORK(isBC+k-1),IWORK(irBC+k-1), IWORK(ivBC+k-1), k=1, nsbc(3))
105     FORMAT (' elem localID: ', I3, ' side:', I1, ' val:', I4)

      ! query vertex BCs
      iveBC = ifree
      ivaBC = iveBC + ndbc(3)

      ierr = iMOAB_GetPointerToVertexBC(pid, ndbc(3), IWORK(iveBC),IWORK(ivaBC))
      call errorout(ierr, 'failed to get vertex boundary conditions')
      ifree = ivaBC + ndbc(3)

      print *, '  Vertex boundary conditions:'
      write (*, 106) (IWORK(iveBC+k-1),IWORK(ivaBC+k-1), k=1, ndbc(3))
106     FORMAT (' vertex: ', I3, ' BC:', I6 )

      endif

#ifdef MOAB_HAVE_MPI
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
      call errorout(ierr, 'fail at barrier')
enddo

      outfile = 'fnew2.h5m'//C_NULL_CHAR
#ifdef MOAB_HAVE_MPI
      wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
#else
      wopts=C_NULL_CHAR
#endif

      !     write out the mesh file to disk
      ierr = iMOAB_WriteMesh(pid, outfile, wopts)
      call errorout(ierr, 'fail to write the mesh file')

      !     all done. de-register and finalize
      ierr = iMOAB_DeregisterApplication(pid)
      call errorout(ierr, 'fail to deregister application')

      ierr = iMOAB_Finalize()
      call errorout(ierr, 'fail to finalize iMOAB')

#ifdef MOAB_HAVE_MPI
      call MPI_FINALIZE ( ierr )
#endif
      call errorout(ierr, 'fail to finalize MPI')

      stop
end   
