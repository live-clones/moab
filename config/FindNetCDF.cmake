#
# Find NetCDF include directories and libraries
#
# NetCDF_INCLUDES            - list of include paths to find netcdf.h
# NetCDF_LIBRARIES           - list of libraries to link against when using NetCDF
# NetCDF_FOUND               - Do not attempt to use NetCDF if "no", "0", or undefined.

set (NetCDF_DIR "" CACHE PATH "Path to search for NetCDF header and library files" )
set (NetCDF_FOUND NO CACHE INTERNAL "Found NetCDF components successfully." )

find_path( NetCDF_INCLUDE_DIR netcdf.h
  ${NetCDF_DIR}
  ${NetCDF_DIR}/include
  /usr/local/include
  /usr/include
)

find_library( NetCDF_C_LIBRARY
  NAMES netcdf
  HINTS ${NetCDF_DIR}
  ${NetCDF_DIR}/lib64
  ${NetCDF_DIR}/lib
  /usr/local/lib64
  /usr/lib64
  /usr/lib64/netcdf
  /usr/local/lib
  /usr/lib
  /usr/lib/netcdf
)

find_library( NetCDF_CXX_LIBRARY
  NAMES netcdf_c++
  HINTS ${NetCDF_DIR}
  ${NetCDF_DIR}/lib64
  ${NetCDF_DIR}/lib
  /usr/local/lib64
  /usr/lib64
  /usr/lib64/netcdf
  /usr/local/lib
  /usr/lib
  /usr/lib/netcdf
)

find_library( NetCDF_FORTRAN_LIBRARY
  NAMES netcdf_g77 netcdf_ifc netcdf_x86_64
  HINTS ${NetCDF_DIR}
  ${NetCDF_DIR}/lib64
  ${NetCDF_DIR}/lib
  /usr/local/lib64
  /usr/lib64
  /usr/lib64/netcdf
  /usr/local/lib
  /usr/lib
  /usr/lib/netcdf
)

IF (NOT NetCDF_FOUND)
  if ( NetCDF_INCLUDE_DIR AND NetCDF_C_LIBRARY )
    set( NetCDF_FOUND YES )
    set(NetCDF_INCLUDES "-I${NetCDF_INCLUDE_DIR}")
    set(NetCDF_LIBRARIES ${NetCDF_C_LIBRARY})
    if ( NetCDF_CXX_LIBRARY )
      set(NetCDF_LIBRARIES ${NetCDF_LIBRARIES} ${NetCDF_CXX_LIBRARY})
    endif ( NetCDF_CXX_LIBRARY )
    if ( NetCDF_FORTRAN_LIBRARY )
      set(NetCDF_LIBRARIES ${NetCDF_LIBRARIES} ${NetCDF_FORTRAN_LIBRARY})
    endif ( NetCDF_FORTRAN_LIBRARY )
    message (STATUS "---   NetCDF Configuration ::")
    message (STATUS "        INCLUDES  : ${NetCDF_INCLUDES}")
    message (STATUS "        LIBRARIES : ${NetCDF_LIBRARIES}")
  else ( NetCDF_INCLUDE_DIR AND NetCDF_C_LIBRARY )
    set( NetCDF_FOUND NO )
    message("finding NetCDF failed, please try to set the var NetCDF_DIR")
  endif ( NetCDF_INCLUDE_DIR AND NetCDF_C_LIBRARY )
ENDIF (NOT NetCDF_FOUND)

mark_as_advanced(
  NetCDF_DIR
  NetCDF_INCLUDES
  NetCDF_LIBRARIES
)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF "NetCDF not found, check environment variables NetCDF_DIR"
  NetCDF_DIR NetCDF_INCLUDES NetCDF_LIBRARIES)
