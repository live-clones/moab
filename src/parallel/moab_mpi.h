#ifndef MOAB_MPI_H
#define MOAB_MPI_H
#include "moab/MOABConfig.h"

#ifndef __cplusplus
#  include <mpi.h>
#elif !defined(MOAB_MPI_CXX_CONFLICT)
#  ifndef MPICH_IGNORE_CXX_SEEK
#    define MPICH_IGNORE_CXX_SEEK
#  endif
#  include <mpi.h>
#else
#  include <stdio.h>
#  ifdef SEEK_SET
#    undef SEEK_SET
#    ifdef MB_SEEK_SET
#      define MOAB_RESTORE_SEEK_SET
#    endif
#  endif
#  ifdef SEEK_CUR
#    undef SEEK_CUR
#    ifdef MOAB_SEEK_CUR
#      define MOAB_RESTORE_SEEK_CUR
#    endif
#  endif
#  ifdef SEEK_END
#    undef SEEK_END
#    ifdef MOAB_SEEK_END
#      define MOAB_RESTORE_SEEK_END
#    endif
#  endif
#  include <mpi.h>
#  ifdef MOAB_RESTORE_SEEK_SET
#    undef MOAB_RESTORE_SEEK_SET
#    define SEEK_SET MOAB_SEEK_SET
#  endif
#  ifdef MOAB_RESTORE_SEEK_CUR
#    undef MOAB_RESTORE_SEEK_CUR
#    define SEEK_CUR MOAB_SEEK_CUR
#  endif
#  ifdef MOAB_RESTORE_SEEK_END
#    undef MOAB_RESTORE_SEEK_END
#    define SEEK_END MOAB_SEEK_END
#  endif
#endif


#endif
