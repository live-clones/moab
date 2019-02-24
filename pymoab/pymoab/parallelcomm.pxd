"""Implements ParallelComm functionality."""

from pymoab cimport moab
cimport mpi4py.MPI as MPI
from mpi4py.libmpi cimport *

cdef class ParallelComm:

    cdef moab.ParallelComm * inst
    cdef moab.Interface * interface
