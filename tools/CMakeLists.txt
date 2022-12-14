# Documentation
# Optional tools
option ( MOAB_BUILD_MBCONVERT    "Build the MOAB mesh converter tool?"           ON )
option ( MOAB_BUILD_HEXMODOPS    "Build the hexahedral modops tool?"             ON )
option ( MOAB_BUILD_MBSIZE       "Build the mbsize tool?"                        ON )
option ( MOAB_BUILD_MBMERGE      "Build the mbmerge tool?"                       ON )
option ( MOAB_BUILD_MBMEM        "Build the mbmem tool?"                         ON )
option ( MOAB_BUILD_MBSKIN       "Build the mbskin tool?"                        ON )
option ( MOAB_BUILD_MBDEPTH      "Build the mbdepth tool?"                       ON )
option ( MOAB_BUILD_MBTAGPROP    "Build the mbtagprop tool?"                     ON )
option ( MOAB_BUILD_MBGSETS      "Build the mbgsets tool?"                       ON )
option ( MOAB_BUILD_SPHEREDECOMP "Build the sphere decomposition tool?"          ON )
option ( MOAB_BUILD_MBSURFPLOT   "Build the mbsurfplot application?"             ON )
option ( MOAB_BUILD_SPHEREDECOMP "Build the spheredecomp application?"           ON )
option ( MOAB_BUILD_MBPART       "Build the mbpart (Zoltan-based) partitioner?" OFF )
option ( MOAB_BUILD_MBSLAVEPART  "Build the slave partitioner tool?"             ON )
option ( MOAB_BUILD_MBCOUPLER    "Build the mesh coupler tool?"                  ON )
option ( MOAB_BUILD_MBHONODES    "Build the hex8 to hex27 converter tool?"       ON )
option ( MOAB_BUILD_MBUMR        "Build the uniform mesh refinement tool?"       ON )
option ( MOAB_BUILD_MBQUALITY    "Build the quality metrics tool?"               ON )
option ( MOAB_BUILD_MBTEMPEST    "Build the TempestRemap coupler tool?"         OFF )

#if ( MOAB_HAVE_MPI AND MPI_FOUND )
#  add_subdirectory( refiner )
#endif ( )

#if (ENABLE_MESQUITE AND MOAB_HAVE_CGM)
#  add_subdirectory(mesquite)
#endif(ENABLE_MESQUITE AND MOAB_HAVE_CGM)

include_directories(
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src
  ${CMAKE_BINARY_DIR}/src/parallel
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/parallel
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if (MOAB_HAVE_TEMPESTREMAP)
  include_directories( ${TEMPESTREMAP_INCLUDES} ${EIGEN3_DIR} )
endif (MOAB_HAVE_TEMPESTREMAP)

set( TOOLS_LINK_LIBS MOAB ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
if ( MOAB_HAVE_CGM )
  set( TOOLS_LINK_LIBS ${TOOLS_LINK_LIBS} ${CGM_LIBRARIES} )
  add_subdirectory( geometry )
endif ( MOAB_HAVE_CGM )

# COMPAREFILES
if ( MOAB_HAVE_HDF5 )
  add_executable( compareFiles compareFiles.cpp)
  set_source_files_properties( compareFiles.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( compareFiles ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_COMPAREFILES 1 )
  install(
    TARGETS compareFiles
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( MOAB_HAVE_HDF5 )

# MBSIZE
if ( MOAB_BUILD_MBSIZE )
  add_executable( mbsize size.cpp measure.cpp)
  set_source_files_properties( size.cpp measure.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbsize ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBSIZE 1 )
  install(
    TARGETS mbsize
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  SET(MOAB_DOCS ${MOAB_DOCS} README.mbsize)
endif ( )

# MBMERGE
if ( MOAB_BUILD_MBMERGE )
  add_executable( mbmerge merge.cpp)
  set_source_files_properties( merge.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbmerge ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBMERGE 1 )
  install(
    TARGETS mbmerge
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( )



# MBMEM
if ( MOAB_BUILD_MBMEM )
  add_executable( mbmem mbmem.cpp)
  set_source_files_properties( mbmem.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbmem ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBMEM 1 )
  install(
    TARGETS mbmem
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( )

# MBCONVERT
if ( MOAB_BUILD_MBCONVERT )
  add_executable( mbconvert convert.cpp)
  set_source_files_properties( convert.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbconvert ${TOOLS_LINK_LIBS} ${TEMPESTREMAP_LIBRARIES} )
  set( MOAB_HAVE_MBCONVERT 1 )
  install(
    TARGETS mbconvert
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  SET(MOAB_DOCS ${MOAB_DOCS} README.mbconvert)
  INSTALL(FILES mbconvert.man DESTINATION share/man/man1 RENAME mbconvert.1)
endif ( )

# MBHONODES
if ( MOAB_BUILD_MBHONODES )
  add_executable( mbhonodes honodes.cpp)
  set_source_files_properties( honodes.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbhonodes ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBHONODES 1 )
  install(
    TARGETS mbhonodes
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( )

# MBDEPTH
if ( MOAB_BUILD_MBDEPTH )
  add_executable( mbdepth depth.cpp)
  set_source_files_properties( depth.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbdepth ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBDEPTH 1 )
  install(
    TARGETS mbdepth
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( )

# MBSKIN
if ( MOAB_BUILD_MBSKIN )
  add_executable( mbskin skin.cpp)
  set_source_files_properties( skin.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbskin ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBSKIN 1 )
  install(
    TARGETS mbskin
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  SET(MOAB_DOCS ${MOAB_DOCS} README.mbskin)
endif ( )

# MBSURFPLOT
if ( MOAB_BUILD_MBSURFPLOT )
  add_executable( mbsurfplot surfplot.cpp )
  set_source_files_properties( surfplot.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbsurfplot ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBSURFPLOT 1 )
  install(
    TARGETS mbsurfplot
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  SET(MOAB_DOCS ${MOAB_DOCS} README.mbsurfplot)
endif ( )

# MBTAGPROP
if ( MOAB_BUILD_MBTAGPROP )
  add_executable( mbtagprop parse.cpp propagate_tags.cpp  )
  set_source_files_properties( parse.cpp propagate_tags.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbtagprop ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBTAGPROP 1 )
  install(
    TARGETS mbtagprop
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  SET(MOAB_DOCS ${MOAB_DOCS} README.mbtagprop)
endif ( )

# MBGSETS
if ( MOAB_BUILD_MBGSETS )
  add_executable( mbgsets gsets.cpp )
  set_source_files_properties( gsets.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbgsets ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBGSETS 1 )
  install(
    TARGETS mbgsets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif( )

# MBUMR
if ( MOAB_BUILD_MBUMR )
  add_executable( mbumr umr.cpp )
  set_source_files_properties( umr.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbumr ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBUMR 1 )
  install(
    TARGETS mbumr
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif( )

# MBSLAVEPART
if ( MOAB_BUILD_MBSLAVEPART )
  add_executable( mbslavepart mbslavepart.cpp)
  set_source_files_properties( mbslavepart.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbslavepart ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBSLAVEPART 1 )
  install(
    TARGETS mbslavepart
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif( )

# MBTEMPEST
if ( MOAB_BUILD_MBTEMPEST )
  add_executable( mbtempest mbtempest.cpp)
  set_source_files_properties( mbtempest.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbtempest ${TOOLS_LINK_LIBS} ${TEMPESTREMAP_LIBRARIES} )
  set( MOAB_HAVE_MBTEMPEST 1 )
  install(
    TARGETS mbtempest
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

endif( )

# MBQUALITY
if ( MOAB_BUILD_MBQUALITY )
  add_executable( mbquality quality.cpp )
  set_source_files_properties( quality.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbquality ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_MBQUALITY 1 )
  install(
    TARGETS mbquality
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif( )

# Zoltan
if ( MOAB_BUILD_MBPART OR ( MOAB_HAVE_ZOLTAN OR MOAB_HAVE_METIS ) )
  set(MOAB_BUILD_MBPART ON)
  # Build the mbpart executable
  add_executable( mbpart mbpart.cpp )
  include_directories( ${ZOLTAN_INCLUDES} ${PARMETIS_INCLUDES} ${METIS_INCLUDES} )
  if (MOAB_HAVE_CGM)
    include_directories( ${CGM_INCLUDES} )
  endif (MOAB_HAVE_CGM)
  set_source_files_properties( mbpart.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( mbpart ${TOOLS_LINK_LIBS} ${ZOLTAN_LIBRARIES} ${PARMETIS_LIBRARIES} ${METIS_LIBRARIES} )
  set( MOAB_HAVE_MBPART 1 )
  install(
    TARGETS mbpart
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( MOAB_BUILD_MBPART OR ( MOAB_HAVE_ZOLTAN OR MOAB_HAVE_METIS ) )

# MBCoupler
if ( ENABLE_IMESH AND MOAB_HAVE_MPI AND MPI_FOUND AND MOAB_BUILD_MBCOUPLER )
  add_subdirectory( mbcoupler )
endif ( )

# hexmodops
if( MOAB_BUILD_HEXMODOPS )
  add_executable( hexmodops makeops.cpp )
  target_link_libraries( hexmodops ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_HEXMODOPS 1 )
  install(
    TARGETS hexmodops
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  SET(MOAB_DOCS ${MOAB_DOCS} README.hexmodops)
endif()

# MBSPHEREDECOMP
if ( MOAB_BUILD_SPHEREDECOMP )
  add_executable( spheredecomp SphereDecomp.cpp SphereDecomp.hpp sphere_decomp.cpp )
  set_source_files_properties( SphereDecomp.cpp SphereDecomp.hpp sphere_decomp.cpp COMPILE_FLAGS "${MOAB_DEFINES}" )
  target_link_libraries( spheredecomp ${TOOLS_LINK_LIBS} )
  set( MOAB_HAVE_SPHEREDECOMP 1 )
  install(
    TARGETS spheredecomp
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  SET(MOAB_DOCS ${MOAB_DOCS} README.spheredecomp)
endif ( )

install(FILES ${MOAB_DOCS} DESTINATION share/doc/moab)
