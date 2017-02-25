"""Implements core functionality."""
from cython.operator cimport dereference as deref

cimport numpy as np
import numpy as np
import ctypes

from pymoab cimport moab
from .tag cimport Tag, TagArray
from .range cimport Range
from .types import check_error, np_tag_type, validate_type
from . import types
from libcpp.vector cimport vector
from libc.stdlib cimport malloc

cdef void* null = NULL


cdef class Core(object):

#    cdef moab.Core *inst

    def __cinit__(self):
        self.inst = new moab.Core()

    def __del__(self):
        del self.inst

    def impl_version(self):
        """MOAB implementation number as a float."""
        return self.inst.impl_version()

    def load_file(self, str fname, exceptions = ()):
        cfname = fname.decode()
        cdef const char * file_name = cfname
        cdef moab.ErrorCode err = self.inst.load_file(fname)
        check_error(err, exceptions)

    def write_file(self, str fname, exceptions = ()):
        """Writes the MOAB data to a file."""
        cfname = fname.decode()
        cdef const char * file_name = cfname
        cdef moab.ErrorCode err = self.inst.write_file(fname)
        check_error(err, exceptions)

    def create_meshset(self, unsigned int options = 0x02, exceptions = ()):
        cdef moab.EntityHandle ms_handle = 0
        cdef moab.EntitySetProperty es_property = <moab.EntitySetProperty> options
        cdef moab.ErrorCode err = self.inst.create_meshset(es_property, ms_handle)
        check_error(err, exceptions)
        return ms_handle

    def add_entities(self, moab.EntityHandle ms_handle, entities, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        if isinstance(entities, Range):
           r = entities
           err = self.inst.add_entities(ms_handle, deref(r.inst))
        else:
           arr = entities
           err = self.inst.add_entities(ms_handle, <unsigned long*> arr.data, len(entities))
        check_error(err, exceptions)

    def create_vertices(self, np.ndarray[np.float64_t, ndim=1] coordinates, exceptions = ()):
        cdef Range rng = Range()
        cdef moab.ErrorCode err = self.inst.create_vertices(<double *> coordinates.data,
                                                            len(coordinates)//3,
                                                            deref(rng.inst))
        check_error(err, exceptions)
        return rng

    def create_element(self, int t, np.ndarray[np.uint64_t, ndim=1] connectivity, exceptions = ()):
        cdef moab.EntityType typ = <moab.EntityType> t
        cdef moab.EntityHandle handle = 0
        cdef int nnodes = len(connectivity)
        cdef moab.ErrorCode err = self.inst.create_element(typ,
            <unsigned long*> connectivity.data, nnodes, handle)
        check_error(err, exceptions)
        return handle

    def create_elements(self, int t, np.ndarray[np.uint64_t, ndim=2] connectivity, exceptions = ()):
        cdef int i
        cdef moab.ErrorCode err
        cdef moab.EntityType typ = <moab.EntityType> t
        #cdef moab.EntityHandle handle = 0
        cdef int nelems = connectivity.shape[0]
        cdef int nnodes = connectivity.shape[1]
        cdef np.ndarray[np.uint64_t, ndim=1] connectivity_i
        cdef np.ndarray[np.uint64_t, ndim=1] handles = np.empty(nelems, 'uint64')
        for i in range(nelems):
            connectivity_i = connectivity[i]
            err = self.inst.create_element(typ, <unsigned long*> connectivity_i.data,
                                           nnodes, deref((<unsigned long*> handles.data)+i))
            check_error(err, exceptions)
        return handles

    def tag_get_handle(self,
                       const char* name,
                       size = None,
                       tag_type = None,
                       create_if_missing = False,
                       storage_type = types.MB_TAG_DENSE,
                       exceptions = ()):
        cdef Tag tag = Tag()
        cdef moab.ErrorCode err
        cdef moab.DataType tt
        cdef int s
        err = self.inst.tag_get_handle(name, tag.inst)
        if err == types.MB_TAG_NOT_FOUND and create_if_missing:
            if tag_type is None or size is None:
                print "ERROR: Not enough information provided to create tag."
                raise ValueError
            else:
                tt = tag_type
                s = size
                err = self.inst.tag_get_handle(name, s, tt, tag.inst, storage_type|types.MB_TAG_CREAT)
        check_error(err, exceptions)
        return tag

    def tag_set_data(self, Tag tag, entity_handles, np.ndarray data, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        cdef moab.DataType tag_type = moab.MB_MAX_DATA_TYPE
        err = self.inst.tag_get_data_type(tag.inst, tag_type);
        check_error(err, ())
        cdef int length = 0
        err = self.inst.tag_get_length(tag.inst,length);
        check_error(err,())
        data = validate_type(tag_type,length,data)
        if isinstance(entity_handles,Range):
            r = entity_handles
            err = self.inst.tag_set_data(tag.inst, deref(r.inst), <const void*> data.data)
            check_error(err, exceptions)
        elif isinstance(entity_handles,np.ndarray):
            assert entity_handles.dtype == 'uint64'
            arr = entity_handles
            err = self.inst.tag_set_data(tag.inst, <unsigned long*> arr.data, len(entity_handles), <const void*> data.data)
            check_error(err, exceptions)
        else:
            check_error(types.MB_FAILURE)

    def tag_get_data(self, Tag tag, entity_handles, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        cdef moab.DataType tag_type = moab.MB_MAX_DATA_TYPE
        err = self.inst.tag_get_data_type(tag.inst, tag_type);
        check_error(err,())
        cdef int length = 0
        err = self.inst.tag_get_length(tag.inst,length);
        check_error(err,())
        cdef np.ndarray data
        if tag_type is types.MB_TYPE_OPAQUE:
            data = np.empty((len(entity_handles),),dtype='S'+str(length))
        else:
            data = np.empty((length*len(entity_handles),),dtype=np.dtype(np_tag_type(tag_type)))
        if isinstance(entity_handles,Range):
            r = entity_handles
            err = self.inst.tag_get_data(tag.inst, deref(r.inst), <void*> data.data)
            check_error(err, exceptions)
        elif isinstance(entity_handles,np.ndarray):
            assert entity_handles.dtype == 'uint64'
            arr = entity_handles
            err = self.inst.tag_get_data(tag.inst, <unsigned long*> arr.data, len(entity_handles), <void*> data.data)
            check_error(err,exceptions)
        else:
            check_error(types.MB_FAILURE)
        return data

    def get_adjacencies(self, entity_handles, int to_dim, bint create_if_missing = False, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        cdef Range adj = Range()
        if isinstance(entity_handles, Range):
            r = entity_handles
            err = self.inst.get_adjacencies(deref(r.inst), to_dim, create_if_missing, deref(adj.inst))
        else:
            arr = entity_handles
            err = self.inst.get_adjacencies(<unsigned long*> arr.data, len(entity_handles), to_dim, create_if_missing, deref(adj.inst))
        check_error(err, exceptions)
        return adj

    def type_from_handle(self, entity_handle):
        cdef moab.EntityType t
        t = self.inst.type_from_handle(<unsigned long> entity_handle)
        return t

    def get_child_meshsets(self, meshset_handle, num_hops = 1, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r = Range()
        err = self.inst.get_child_meshsets(<unsigned long> meshset_handle, deref(r.inst), num_hops)
        check_error(err, exceptions)
        return r

    def add_child_meshset(self, parent_meshset, child_meshset, exceptions = ()):
        cdef moab.ErrorCode err
        err = self.inst.add_child_meshset(<unsigned long> parent_meshset, <unsigned long> child_meshset)
        check_error(err, exceptions)

    def get_root_set(self):
        return <unsigned long> 0

    def get_coords(self, entities, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range r
        cdef np.ndarray[np.uint64_t, ndim=1] arr
        cdef np.ndarray coords
        if isinstance(entities, Range):
            r = entities
            coords = np.empty((3*r.size(),),dtype='float64')
            err = self.inst.get_coords(deref(r.inst), <double*> coords.data)
        else:
            arr = entities
            coords = np.empty((3*len(arr),),dtype='float64')
            err = self.inst.get_coords(<unsigned long*> arr.data, len(entities), <double*> coords.data)
        check_error(err, exceptions)
        return coords

    def get_entities_by_type(self, meshset, t, bint recur = False, exceptions = ()):
        cdef moab.ErrorCode err
        cdef vector[moab.EntityHandle] entities
        cdef moab.EntityType typ = t
        err = self.inst.get_entities_by_type(<unsigned long> meshset,
                                             typ,
                                             entities,
                                             recur)
        check_error(err, exceptions)
        return entities

    def get_entities_by_type_and_tag(self, meshset, t, tags, np.ndarray vals, int cond = 0, bint recur = False, exceptions = ()):
        cdef moab.ErrorCode err
        assert len(tags) == len(vals)
        cdef int num_tags = len(tags)
        cdef moab.EntityType typ = t
        cdef TagArray ta = TagArray(tags)
        cdef moab.DataType tag_type = moab.MB_MAX_DATA_TYPE
        cdef Range ents = Range()
        cdef np.ndarray tmp
        #create temporary array
        cdef int num_none = 0
        for val in vals:
            if val == None:
                num_none +=1
        cdef char** vals_arr = <char**> malloc(num_tags*sizeof(char*))
        for i in range(num_tags):
            if vals[i] == None:
                vals_arr[i] =  NULL
            else:
                tmp = vals[i:i+1]
                vals_arr[i] = tmp.data
        #setup vectors to hold data
        cdef vector[int] int_vec
        int_vec.resize(num_tags)
        cdef vector[double] double_vec
        double_vec.resize(num_tags)
        cdef void** arr = <void**> malloc(num_tags*sizeof(void*))
        #get the tag type
        cdef moab.DataType this_tag_type = moab.MB_MAX_DATA_TYPE
        cdef Tag this_tag
        for i in range(num_tags):
            # if None is passed, set pointer to NULL and continue
            if vals[i] == None:
                arr[i] = NULL
                continue
            # otherwise get the tag type
            this_tag = tags[i]
            err = self.inst.tag_get_data_type(this_tag.inst, this_tag_type)
            check_error(err)
            # cast tag value as type
            if this_tag_type == types.MB_TYPE_INTEGER:
                int_vec[i] = <int> vals[i]
                arr[i] = &(int_vec[i])
            if this_tag_type == types.MB_TYPE_DOUBLE:
                double_vec[i] = <double> vals[i]
                arr[i] = &(double_vec[i])
        #here goes nothing
        err = self.inst.get_entities_by_type_and_tag(<unsigned long> meshset, typ, ta.ptr, <const void**> arr, len(tags), deref(ents.inst), cond, recur)
        check_error(err, exceptions)
        return ents

    def get_entities_by_handle(self, meshset, bint recur = False, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range ents = Range()
        err = self.inst.get_entities_by_handle(<unsigned long> meshset, deref(ents.inst), recur)
        check_error(err, exceptions)
        return ents

    def get_entities_by_dimension(self, meshset, int dimension, bint recur = False, exceptions = ()):
        cdef moab.ErrorCode err
        cdef Range ents = Range()
        err = self.inst.get_entities_by_dimension(<unsigned long> meshset, dimension, deref(ents.inst), recur)
        check_error(err, exceptions)
        return ents
