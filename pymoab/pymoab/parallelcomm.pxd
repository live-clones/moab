"""Implements ParallelComm functionality."""

from pymoab cimport moab
from .core cimport Core

cdef class ParallelComm:

    cdef moab.ParallelComm * inst
    cdef moab.Interface * interface
    cdef Core mbCore


from libcpp.string cimport string as std_string

# cdef extern from "moab/ParallelComm.hpp":

#     cdef std_string PARALLEL_PARTITION_SET_TAG_NAME

