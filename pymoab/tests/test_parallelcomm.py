from pymoab import core
from pymoab import parallelcomm
from pymoab import types
from pymoab.rng import Range
from pymoab.scd import ScdInterface
from pymoab.hcoord import HomCoord
from subprocess import call
from driver import test_driver, CHECK, CHECK_EQ, CHECK_NOT_EQ, CHECK_ITER_EQ

from mpi4py import MPI
import numpy as np
import os

bytes_per_char_ = np.array(["a"]).nbytes

def test_parallel_load_mesh():
    mb = core.Core()
    pc = parallelcomm.ParallelComm(mb, comm=MPI.COMM_WORLD)
    try:
        mb.load_file("parallel_file.h5m", file_set = None, readopts = 'PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS')
    except:
        try:
            print( """
            WARNING: .h5m file load failed. If hdf5 support is enabled in this
            build there could be a problem.
            """)
            mb.load_file("parallel_file.vtk", file_set = 0, readopts = '')
        except:
            raise(IOError, "Failed to load MOAB file.")

    #load into file_set
    mb1 = core.Core()
    pc1 = parallelcomm.ParallelComm(mb1, comm=MPI.COMM_WORLD)
    file_set = mb1.create_meshset()
    try:
        mb1.load_file("parallel_file.h5m", file_set, readopts = 'PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS')
    except:
        try:
            print("""
            WARNING: .h5m parallel file load failed. If HDF5 support is enabled in this
            build there could be a problem.
            """)
            mb1.load_file("parallel_file.vtk",file_set, readopts = '')
        except:
            raise(IOError, "Failed to parallel load MOAB file.")

    ents = mb1.get_entities_by_type(file_set,types.MBMAXTYPE)
    CHECK_NOT_EQ(len(ents),0)

def test_parallel_write_mesh():
    mb = core.Core()
    pc = parallelcomm.ParallelComm(mb, comm=MPI.COMM_WORLD)
    mb.create_vertices(np.ones(3))

    try:
        mb.write_file("poutfile.h5m", output_sets = None, writeopts = 'PARALLEL=WRITE_PART')
        assert os.path.isfile("poutfile.h5m")
    except:
        try:
            print("""
            WARNING: .h5m file write failed. If HDF5 support is enabled in this
            build there could be a problem.
            """)
            mb.write_file("poutfile.vtk", output_sets = None, writeopts = '')
            assert os.path.isfile("outfile.vtk")
        except:
            raise(IOError, "Failed to parallel write MOAB file.")


def test_parallel_write_tags():
    """
    Test write tag functionality
    """

    # test values
    outfile = "parallel_write_tag_test.h5m"

    mb = core.Core()
    pc = parallelcomm.ParallelComm(mb, comm=MPI.COMM_WORLD)
    rank = pc.rank()

    coords = np.array((rank+0,rank+0,rank+0,rank+1,rank+0,rank+0,rank+1,rank+1,rank+1),dtype='float64')
    verts = mb.create_vertices(coords)
    CHECK_EQ(len(verts),3)
    #create elements
    verts = np.array(((verts[0],verts[1],verts[2]),),dtype='uint64')
    tris = mb.create_elements(types.MBTRI,verts)
    CHECK_EQ(len(tris),1)

    #check that the element is there via GLOBAL_ID tag
    global_id_tag = mb.tag_get_handle(types.GLOBAL_ID_TAG_NAME,1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
    gids = np.array(rank+0, rank+1, rank+2, dtype='uint64')
    mb.tag_set_data(global_id_tag, verts, gids)

    pc.resolve_shared_ents(mb.get_root_set(), tris, 2, 1)

    pset = pc.create_part()
    mb.add_entities(pset, vs)
    pc.assign_global_ids(setid=pset, dimension=0, startid=5, largestdimonly=False)
    rank = pc.rank()
    mb.tag_set_data(pc.partition_tag(), pset, rank)

    # create writing tag
    write_tag = mb.tag_get_handle("WRITE",
                                  3,
                                  types.MB_TYPE_DOUBLE,
                                  types.MB_TAG_DENSE,
                                  create_if_missing=True)
    # set some data on that tag
    data = [0.7071, 0.7071, 0.0]
    mb.tag_set_data(write_tag, vs, data)

    mb.write_file(outfile, output_tags = [write_tag], writeopts = 'PARALLEL=WRITE_PART')

    mb2 = core.Core()
    mb2.load_file(outfile, readopts = 'PARALLEL=READ_PART')

    vs = mb2.get_entities_by_type(0, types.MBVERTEX)

    # get the write tag
    new_write_tag = mb2.tag_get_handle("WRITE")

    # make sure we can still get data for the write tag
    d = mb2.tag_get_data(new_write_tag, vs)

    # make sure the second tag is not there
    try:
        no_write_tag = mb2.tag_get_handle("NO_WRITE")
        raise AssertionError("Tag get handle succeeded when it should not.")
    except(RuntimeError):
        pass

    # write multiple tags
    mb.write_file(outfile, output_tags = [write_tag, no_write_tag], writeopts = 'PARALLEL=WRITE_PART')

    mb2 = core.Core()
    pc2 = parallelcomm.ParallelComm(mb2, comm=MPI.COMM_WORLD)
    mb2.load_file(outfile, readopts = 'PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS')

    vs = mb2.get_entities_by_type(0, types.MBVERTEX)

    # get the write tag
    new_write_tag = mb2.tag_get_handle("WRITE")

    # make sure we can still get data for the write tag
    d = mb2.tag_get_data(new_write_tag, vs)

    # make sure the tag is not there
    new_no_write_tag = mb2.tag_get_handle("NO_WRITE")

    # make sure we can now get data for the "no-write" tag
    d = mb2.tag_get_data(new_no_write_tag, vs)


def test_parallel_delete_mesh():
    mb = core.Core()
    mb.create_vertices(np.ones(9))
    rs = mb.get_root_set()
    ents = mb.get_entities_by_handle(rs)
    CHECK_EQ(len(ents),3)
    # now delete all mesh entities
    mb.delete_mesh()
    ents = mb.get_entities_by_handle(rs)
    CHECK_EQ(len(ents),0)


if __name__ == "__main__":
    tests = [test_parallel_load_mesh,
             test_parallel_write_mesh,
             test_parallel_write_tags,
             test_parallel_delete_mesh
             ]
    test_driver(tests)
