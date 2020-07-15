"""Implements parallelcomm functionality."""

cdef extern from "mpi-compat.h": pass

cimport mpi4py.MPI as MPI
import mpi4py.MPI as MPI
from mpi4py.libmpi cimport *


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
from libcpp cimport bool


#from collections import Iterable

cdef void *null = NULL

# from pymoab cimport parallelcomm
# PARALLEL_PARTITION_SET_TAG_NAME    = str(parallelcomm.PARALLEL_PARTITION_SET_TAG_NAME.decode())

cdef class ParallelComm(object):

    def __cinit__(self, Core c, MPI.Comm comm, pid=None):
        """
        Constructor.

        Requires a moab core, c, to operate on.
        Requires a MPI Comm instance to build the ParallelComm.

        If no argument is provided, an empty ParallelComm will be created and returned.
        """
        self.mbCore = c
        self.interface  = <moab.Interface*> c.inst
        # self.comm = comm

        # Get the C interface copy of the comm instance
        cdef MPI_Comm c_comm = comm.ob_mpi
        self.inst = new moab.ParallelComm(self.interface, c_comm)
        # self.partTag = self.inst.partition_tag()


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

    def moab_instance(self):
        """
        Get the current MPI rank in the parallel comm instance
        """
        return self.mbCore

    # def comm(self):
    #     """
    #     Get the current MPI comm instance
    #     """
    #     cdef MPI_Comm c_comm = self.inst.comm()
    #     return self.comm


    def assign_global_ids(self, setid, int dimension=3, int startid=1, bint largestdimonly=True, bint isparallel=True, bint ownedonly=False, exceptions = ()):
        """
        Invoke the algorithm to assign global IDs in parallel based on dimension of entities and 
        whether it is filtered over owned/shared entities.
        """
        # cdef moab.EntityHandle setHandle = setid
        cdef moab.ErrorCode err
        err = self.inst.assign_global_ids(<unsigned long> setid, dimension, startid, largestdimonly, isparallel, ownedonly)
        check_error(err, exceptions)
        # return err

    def resolve_shared_ents (self, setid, proc_ents, int resolve_dim=3, int shared_dim=-1, Tag idTag=None, exceptions = ()):
        """
        Resolve shared entities in parallel
        """
        cdef Range r
        cdef moab.Tag* gidTag = NULL
        cdef moab.ErrorCode err
        if isinstance(proc_ents, _eh_py_type):
            proc_ents = Range(proc_ents)
        r = proc_ents
        if idTag is not None:
            gidTag = &(idTag.inst)
        err = self.inst.resolve_shared_ents(<unsigned long> setid, deref(r.inst), resolve_dim, shared_dim, NULL, gidTag)
        check_error(err, exceptions)

    def resolve_shared_entities (self, setid, int resolve_dim=3, int shared_dim=-1, Tag idTag=None, exceptions = ()):
        """
        Resolve shared entities in parallel
        """
        cdef moab.ErrorCode err
        cdef moab.Tag* gidTag = NULL
        if idTag is not None:
            gidTag = &(idTag.inst)
        err = self.inst.resolve_shared_ents(<unsigned long> setid, resolve_dim, shared_dim, gidTag)
        check_error(err, exceptions)

    def create_part(self, exceptions = ()):
        """
        Create a new part and return EntityHandle
        """
        cdef moab.EntityHandle partHandle = 0
        cdef moab.ErrorCode err = self.inst.create_part(partHandle)
        check_error(err, exceptions)
        # coreobj = Core(self.inst.get_moab())
        # rank = self.inst.rank()
        # err = coreobj.tag_set_data(self.partition_tag(), partHandle, rank)
        # check_error(err, exceptions)
        return _eh_py_type(partHandle)

    def partition_tag(self):
        """
        Return the partition tag for PComm object
        """
        ret_tag = self.mbCore.tag_get_handle("PARALLEL_PARTITION")
        return ret_tag

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

