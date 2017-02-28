"""Implements range functionality."""
from cython.operator cimport dereference as deref

from pymoab cimport moab

cdef class Range(object):

    #def __cinit__(self, moab.EntityHandle val1=None, moab.EntityHandle val2=None):
    #    if val1 is None or val2 is None:
    #        self.inst = new moab.Range()
    #    else:
    #        self.inst = new moab.Range(val1, val2)

    def __cinit__(self):
        self.inst = new moab.Range()

    def __del__(self):
        del self.inst

    def size(self):
        """The number of values this Ranges represents."""
        return len(self)

    def __len__(self):
        return self.inst.size()

    def psize(self):
        """The number of range pairs in the list."""
        return self.inst.psize()

    def empty(self):
        """Is the range empty?"""
        return self.inst.empty()

    def clear(self):
        """clears the contents of the list."""
        self.inst.clear()

    def __iter__(self):
        cdef int i = 0
        for i in range(0, self.inst.size()):
            yield self[i]

    def __getitem__(self, int index):
        cdef moab.EntityHandle rtn = deref(self.inst)[index]
        if index < self.size():
            return rtn
        else:
            raise StopIteration

    def __str__(self):
        prefix= "["
        delim = ", "
        suffix= "]"
        out_str=prefix
        for entity in self:
            out_str+=str(entity)
            out_str+=delim
        #replace final delimiter
        out_str=out_str[:-2]
        out_str+=suffix
        return out_str
            
