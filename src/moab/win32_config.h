#ifndef MOAB_WIN32_CONFIG_H
#define MOAB_WIN32_CONFIG_H

#ifdef WIN32
#include "MOAB_export.h"

#define _USE_MATH_DEFINES // for C++  
#include <cmath> 

#else
#define MOAB_EXPORT
#endif

#endif
