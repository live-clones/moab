#!/usr/bin/env python
# add to the grid.h5m data some data from land files D2d.clm2.r.0001-01-02-00000.nc
# 

import sys
import os
import netCDF4
import numpy as np
import math

from optparse import OptionParser

from pymoab import core
from pymoab import types
from pymoab.rng import Range


print "== Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser()
parser.description = "This script takes a grid.h5m MOAB file and adds some data from land output"

parser.add_option("-g", "--grid", dest="grid", help="grid.h5m grid file to input.", default="grid.h5m", metavar="FILENAME")
parser.add_option("-f", "--landf", dest="land", help="land nc file to process.", default="/home/iulian/acme/scratch/D2d/run/D2d.clm2.r.0001-01-02-00000.nc", metavar="FILENAME")
parser.add_option("-m", "--moab", dest="moabFile", help="MOAB grid file name used as output.", default="outland.h5m", metavar="FILENAME")

for option in parser.option_list:
        if option.default != ("NO", "DEFAULT"):
                option.help += (" " if option.help else "") + "[default: %default]"
                
options, args = parser.parse_args()

if not options.grid:
        sys.exit('Error: input grid file is required.  Specify with -g command line argument.')
if not options.land:
        sys.exit('Error: input land result file is required.  Specify with -f command line argument.')
if not options.moabFile:
        sys.exit('Error: MOAB output file is required.  Specify with -m command line argument.')
        

fin = netCDF4.Dataset(options.land, 'r')
#fout = netCDF4.Dataset(options.moabFile, 'w')  # This will clobber existing files
mb = core.Core()

# read moab file
mb.load_file(options.grid)
# Dimensions
land_size=len(fin.dimensions["gridcell"])

grid1d_ixy=fin.variables['grid1d_ixy'][:]

mask_tag = mb.tag_get_handle("MASK",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)

polys = mb.get_entities_by_type(0, types.MBPOLYGON)
print 'number of polygons:', len(polys)

npolys = len(polys)
data = np.zeros(npolys)

for e in xrange(land_size):
  data[grid1d_ixy[e]-1]= 1
  
mb.tag_set_data(mask_tag,polys,data)



try:
        #mb.write_file(options.moabFile)
        mb.write_file(options.moabFile)
        assert os.path.isfile(options.moabFile)
except:
        try:
            print("""
            WARNING: .h5m file write failed. If hdf5 support is enabled in this
            build there could be a problem.
            """)
            mb.write_file("outfile.vtk")
            assert os.path.isfile("outfile.vtk")
        except:
            raise(IOError, "Failed to write MOAB file.")
