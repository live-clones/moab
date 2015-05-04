#ifndef MBProcConfig_HEADER
#define MBProcConfig_HEADER

#include "MBTypes.h"
#include "MBRange.hpp"
#ifdef MOAB_HAVE_MPI
#  include "MBmpi.h"
#endif

#include "moab/ProcConfig.hpp"
typedef moab::ProcConfig MBProcConfig;

#endif
