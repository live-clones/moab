set( moab_PARALLEL_SRCS
  ParallelComm.cpp
  ProcConfig.cpp
  ParallelData.cpp
  ParallelMergeMesh.cpp
  ReadParallel.hpp      ReadParallel.cpp
  SharedSetData.hpp     SharedSetData.cpp
  ParCommGraph.cpp
  gs.cpp
)

set(moab_PARALLEL_INSTALL_HDRS moab/gs.hpp
                               moab_mpi.h
                               moab/ParallelComm.hpp
                               moab/ParallelMergeMesh.hpp
                               moab/ProcConfig.hpp
                               moab/ParallelData.hpp
                               moab/ParCommGraph.hpp
                               MBParallelConventions.h )

include_directories(
    ${MOAB_BINARY_DIR}/src
    ${MOAB_SOURCE_DIR}/src
    ${MOAB_SOURCE_DIR}/src/parallel
    ${CMAKE_CURRENT_BINARY_DIR}
  )

if ( MOAB_HAVE_HDF5 AND MOAB_HAVE_HDF5_PARALLEL )
  set( moab_PARALLEL_SRCS ${moab_PARALLEL_SRCS}
                          WriteHDF5Parallel.cpp
                          WriteHDF5Parallel.hpp )
endif()

if ( MOAB_HAVE_HDF5 )
  include_directories(
    ${HDF5_INCLUDES}
    ${MOAB_SOURCE_DIR}/mhdf/include
  )
endif ( MOAB_HAVE_HDF5 )

if ( MOAB_HAVE_ZOLTAN )
  include_directories( ${ZOLTAN_INCLUDE_DIR} )
endif ( MOAB_HAVE_ZOLTAN )

set_source_files_properties( ${moab_PARALLEL_SRCS}
  COMPILE_FLAGS "-DIS_BUILDING_MB -DSRCDIR=\"${MOAB_SOURCE_DIR}\" ${MOAB_DEFINES}"
)

add_library(
  moabpar
  OBJECT ${moab_PARALLEL_INSTALL_HDRS} ${moab_PARALLEL_SRCS}
  )

include(AutoconfHeader)
moab_install_headers(${moab_PARALLEL_INSTALL_HDRS})

