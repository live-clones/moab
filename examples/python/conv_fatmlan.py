#!/usr/bin/env python
# Create a moab h5m file from a fatmlndfrc file
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
parser.description = "This script takes a fatmlndfrc file and generates a MOAB h5m file."

parser.add_option("-f", "--falf", dest="falf", help="fatmlndfrc grid file to input.", default="domain.lnd.ne11np4_oQU240.160614.nc", metavar="FILENAME")
parser.add_option("-m", "--moab", dest="moabFile", help="MOAB grid file name used as output.", default="grid.h5m", metavar="FILENAME")

for option in parser.option_list:
        if option.default != ("NO", "DEFAULT"):
                option.help += (" " if option.help else "") + "[default: %default]"
                
options, args = parser.parse_args()

if not options.falf:
        sys.exit('Error: atm land fraction input grid file is required.  Specify with -f command line argument.')
if not options.moabFile:
        sys.exit('Error: MOAB output grid file is required.  Specify with -m command line argument.')
        

fin = netCDF4.Dataset(options.falf, 'r')
#fout = netCDF4.Dataset(options.moabFile, 'w')  # This will clobber existing files
mb = core.Core()
# Write to output file
# Dimensions
grid_size=len(fin.dimensions["n"])
grid_corners=len(fin.dimensions["nv"])
ni=len(fin.dimensions["ni"])
nj=len(fin.dimensions["nj"])

grid_corner_lat = fin.variables['yv'][:]
grid_corner_lon = fin.variables['xv'][:]

nverts=grid_size*grid_corners
# create vertices
coords=np.zeros(nverts*3)

connect=np.zeros((grid_size,grid_corners), dtype='uint64')
for e in xrange(grid_size):
   for j in xrange(grid_corners):
      latd = grid_corner_lat[0,e, j] /180. * math.pi
      lond = grid_corner_lon[0,e, j] /180. * math.pi
      #  cart%x=sphere%r*COS(sphere%lat)*COS(sphere%lon)
      #  cart%y=sphere%r*COS(sphere%lat)*SIN(sphere%lon)
      #  cart%z=sphere%r*SIN(sphere%lat)
      coords[3*(e*grid_corners+j)]   = math.cos(latd) * math.cos(lond)
      coords[3*(e*grid_corners+j)+1] = math.cos(latd) * math.sin(lond)
      coords[3*(e*grid_corners+j)+2] = math.sin(latd)
      connect[e,j] = e*grid_corners+j+1
      
mb.create_vertices(coords)
mb.create_elements(types.MBPOLYGON, connect)

# set global id for polygons
polys = mb.get_entities_by_type(0,types.MBPOLYGON, False)
data = range(1,grid_size+1)
global_id_tag = mb.tag_get_handle("GLOBAL_ID",1,types.MB_TYPE_INTEGER,types.MB_TAG_DENSE,True)
gids = np.array(data)
mb.tag_set_data(global_id_tag,polys,gids)
# set global id for vertices
data = range(1,nverts+1)
verts = mb.get_entities_by_type(0,types.MBVERTEX, False)
gids = np.array(data)
mb.tag_set_data(global_id_tag,verts,gids)

# set mask, frac,  tag on elements, to help a little in debugging
mask = fin.variables['mask'][:]
frac = fin.variables['frac'][:]
area = fin.variables['area'][:]
mask_array   = np.zeros(grid_size)
frac_array   = np.zeros(grid_size)
area_array   = np.zeros(grid_size)

for e in xrange(grid_size):
  mask_array[e]=mask[0,e]
  area_array[e]=area[0,e]
  frac_array[e]=frac[0,e]


mask_tag = mb.tag_get_handle("mask",1,types.MB_TYPE_DOUBLE ,types.MB_TAG_DENSE,True)
mb.tag_set_data(mask_tag, polys, mask_array)

frac_tag = mb.tag_get_handle("frac",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
mb.tag_set_data(frac_tag, polys, frac_array)
area_tag = mb.tag_get_handle("area",1,types.MB_TYPE_DOUBLE,types.MB_TAG_DENSE,True)
mb.tag_set_data(area_tag, polys, area_array)


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
