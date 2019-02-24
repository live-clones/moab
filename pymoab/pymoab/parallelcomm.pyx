# distutils: language = c++
"""
Implements parallelcomm functionality.
"""

cimport mpi4py.MPI as MPI
from mpi4py.libmpi cimport *

# from mpi4py import MPI
cdef extern from "mpi-compat.h": pass

from pymoab cimport moab
from cython.operator cimport dereference as deref

from .core cimport Core
from .tag cimport Tag, _tagArray
from .rng cimport Range
from .types import check_error, np_tag_type, validate_type, _convert_array, _eh_array, _eh_py_type
#from . import types
from libcpp.vector cimport vector
from libcpp.string cimport string as std_string
from libc.stdlib cimport malloc

#from collections import Iterable

cdef void *null = NULL

cdef class ParallelComm(object):

    def __cinit__(self, Core c, MPI.Comm comm, pid=None):
        """
        Constructor.

        Requires a moab core, c, to operate on.
        Requires a MPI Comm instance to build the ParallelComm.

        If no argument is provided, an empty ParallelComm will be created and returned.
        """
        self.interface  = <moab.Interface*> c.inst
        # self.comm = comm

        # Get the C interface copy of the comm instance
        cdef MPI_Comm c_comm = comm.ob_mpi
        self.inst = new moab.ParallelComm(self.interface, c_comm)


    def __del__(self):
        """
        Destructor.
        """
        del self.inst

    def get_id(self):
        """
        The identifier for this ParallelComm.
        """
        return self.inst.get_id()

    def size(self):
        """
        Get the number of MPI processes in the parallel comm instance
        """
        return self.inst.size()

    def rank(self):
        """
        Get the current MPI rank in the parallel comm instance
        """
        return self.inst.rank()

    # def comm(self):
    #     """
    #     Get the current MPI comm instance
    #     """
    #     cdef MPI_Comm c_comm = self.inst.comm()
    #     return self.comm


    def assign_global_ids(self, setid, int dimension, int startid=1, largestdimonly=True, isparallel=True, ownedonly=False, exceptions = ()):
        """
        Invoke the algorithm to assign global IDs in parallel based on dimension of entities and 
        whether it is filtered over owned/shared entities.
        """
        cdef start_id = startid
        cdef moab.ErrorCode err = self.inst.assign_global_ids(setid, dimension)#, startid, largestdimonly, isparallel, ownedonly)
        check_error(err, exceptions)
        return err

    def create_part(self, exceptions = ()):
        """
        Create a new part and return EntityHandle
        """
        cdef moab.EntityHandle partHandle = 0
        cdef moab.ErrorCode err = self.inst.create_part(partHandle)
        check_error(err, exceptions)
        return _eh_py_type(partHandle)

    def __str__(self):
        """
        ParallelComm as a string
        """
        rank = self.inst.rank()
        size = self.inst.size()
        sid = self.inst.get_id()
        return "[%d] ParallelComm with ID=%d, has total size = %d" % (rank,sid,size)

    def __repr__(self):
        """
        Representation of class as a string
        """
        return self.__str__()

