#ifndef MOAB_MPI_H
#define MOAB_MPI_H
#include "moab_mpi_config.h"
// following is to disable inclusion of an openmpi header that causes a compile warning
// since we don't use c++ mpi bindings, we don't need it, and doing this allows us to mark warnings as errors
#define OMPI_SKIP_MPICXX 1

#ifndef __cplusplus
#include <mpi.h>
#elif !defined( MB_MPI_CXX_CONFLICT )
#ifndef MPICH_IGNORE_CXX_SEEK
#define MPICH_IGNORE_CXX_SEEK
#endif
#include <mpi.h>
#else
#include <stdio.h>
#ifdef SEEK_SET
#undef SEEK_SET
#ifdef MB_SEEK_SET
#define MB_RESTORE_SEEK_SET
#endif
#endif
#ifdef SEEK_CUR
#undef SEEK_CUR
#ifdef MB_SEEK_CUR
#define MB_RESTORE_SEEK_CUR
#endif
#endif
#ifdef SEEK_END
#undef SEEK_END
#ifdef MB_SEEK_END
#define MB_RESTORE_SEEK_END
#endif
#endif
#include <mpi.h>
#ifdef MB_RESTORE_SEEK_SET
#undef MB_RESTORE_SEEK_SET
#define SEEK_SET MB_SEEK_SET
#endif
#ifdef MB_RESTORE_SEEK_CUR
#undef MB_RESTORE_SEEK_CUR
#define SEEK_CUR MB_SEEK_CUR
#endif
#ifdef MB_RESTORE_SEEK_END
#undef MB_RESTORE_SEEK_END
#define SEEK_END MB_SEEK_END
#endif
#endif

#endif
