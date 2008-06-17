/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#include "MBProcConfig.hpp"

//! Constructor
MBProcConfig::MBProcConfig(MPI_Comm proc_comm) 
    : procComm(proc_comm),
      crystalInit(false)
{
#ifdef USE_MPI
  int rank, size;
  MPI_Comm_rank(procComm, &rank); 
  procRank = (unsigned int) rank;
  MPI_Comm_size(procComm, &size); 
  procSize = (unsigned int) size;
#else
  procRank = 0;
  procSize = 1;
#endif
}

crystal_data *MBProcConfig::crystal_router(bool construct_if_missing)
{
  if (!crystalInit && construct_if_missing)
#ifdef USE_MPI
    crystal_init(&crystalData, procComm), crystalInit = true;
#else
  ;
#endif

  return &crystalData;
}

MBProcConfig::~MBProcConfig() 
{
#ifdef USE_MPI
  if (crystalInit) 
    crystal_free(&crystalData), crystalInit = false;
#endif
}

