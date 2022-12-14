include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/io
    ${CMAKE_SOURCE_DIR}/itaps
    ${CMAKE_BINARY_DIR}/itaps
    ${CMAKE_SOURCE_DIR}/itaps/imesh
    ${CMAKE_BINARY_DIR}/itaps/imesh
    ${CMAKE_SOURCE_DIR}/src/parallel
    ${CMAKE_BINARY_DIR}/src/parallel )

add_subdirectory(io)
add_subdirectory(dual)
if( MOAB_HAVE_HDF5 )
  add_subdirectory(h5file)
endif()
add_subdirectory(obb)
add_subdirectory(perf)

set(TEST_COMP_FLAGS "-DMESHDIR=${MOAB_ABSSRC_DIR}/MeshFiles/unittest")

set( TESTS range_test.cpp
           verdict_test.cpp
           scdseq_test.cpp
           scd_test_partn.cpp
           test_adj.cpp
           GeomUtilTests.cpp
           oriented_box_test.cpp
           adaptive_kd_tree_tests.cpp
           kd_tree_test.cpp
           bsp_tree_test.cpp
           reorder_test.cpp
           elem_eval_test.cpp
           VarLenTagTest.cpp
           TagTest.cpp
           spatial_locator_test.cpp
           bsp_tree_poly_test.cpp
           test_prog_opt.cpp
           coords_connect_iterate.cpp
           urefine_mesh_test.cpp
           adj_moab_test.cpp
           test_boundbox.cpp
           spherical_area_test.cpp
           arc_intx_tests.cpp
           intx_in_plane_test.cpp
           intx_on_sphere_test.cpp )

if(MOAB_HAVE_HDF5)
  set( TESTS ${TESTS}
             lloyd_smoother_test.cpp
             mhdf_public_test.cpp
             mergemesh_test.cpp
             mbfacet_test.cpp
             mbground_test.cpp
             gttool_test.cpp
             restore_topo_geom_incl_test.cpp
             crop_vol_test.cpp
	     gqt_simple_test.cpp
	     gqt_rayfire_test.cpp
	     gqt_pointinvol_test.cpp
            gqt_findvol_test.cpp
	     test_geom_gqt.cpp)
  set( HDF5_C_TESTS  imoab_test.c )
endif()

if(MOAB_HAVE_NETCDF)
  set ( TESTS ${TESTS}
              intx_rll_cs_sphere_test.cpp )
endif()

foreach( fname ${HDF5_C_TESTS} )
  string( REPLACE ".c" "" base ${fname} )
  add_executable( ${base} ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
  target_link_libraries( ${base} MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
endforeach()

foreach( fname ${TESTS} )
  string( REPLACE ".cpp" "" base ${fname} )
  add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${fname})
  set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
  target_link_libraries( ${base} MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
  add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
endforeach()

if(ENABLE_PYMOAB)
  SET(PYMOAB_TESTS test_core test_hom test_scd test_range test_tag test_topo_util test_skinner)
  FOREACH(pymoab_test ${PYMOAB_TESTS})
    ADD_TEST(NAME pymoab_${pymoab_test}
      COMMAND "env" "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib/:$ENV{LD_LIBRARY_PATH}" "PYTHONPATH=${PROJECT_BINARY_DIR}/pymoab:$ENV{PYTHONPATH}" "${PYTHON_EXECUTABLE}" "${pymoab_test}.py"
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/pymoab/tests/")
  ENDFOREACH(pymoab_test)
endif()

add_executable( TestTypeSequenceManager ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp TestTypeSequenceManager.cpp)
set_target_properties( TestTypeSequenceManager PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST -DIS_BUILDING_MB" )

set(TestTypeSequenceManager_LIBS MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
if(${Threads_FOUND})
  list(APPEND TestTypeSequenceManager_LIBS Threads::Threads)
endif()
target_link_libraries(TestTypeSequenceManager ${TestTypeSequenceManager_LIBS})

add_test( TestTypeSequenceManager ${EXECUTABLE_OUTPUT_PATH}/TestTypeSequenceManager )

#Test_MBMeshSet.cpp
add_executable( Test_MBMeshSet ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp Test_MBMeshSet.cpp)
set_target_properties( Test_MBMeshSet PROPERTIES COMPILE_FLAGS "${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST -DIS_BUILDING_MB" )
target_link_libraries( Test_MBMeshSet MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
add_test( Test_MBMeshSet ${EXECUTABLE_OUTPUT_PATH}/Test_MBMeshSet )

if(ZCAT_EXE AND NOT WIN32)
  find_program( ZCAT_EXE NAMES zcap ZCAT PATHS /bin /usr/bin )
  add_custom_command( OUTPUT mb_big_test.g
                      COMMAND ${ZCAT_EXE} < ${MOAB_SOURCE_DIR}/MeshFiles/unittest/mb_big_test.g.gz > mb_big_test.g
                      DEPENDS ${MOAB_SOURCE_DIR}/MeshFiles/unittest/mb_big_test.g.gz )
  add_custom_command( OUTPUT cell1.gen
                      COMMAND ${ZCAT_EXE} < ${MOAB_SOURCE_DIR}/MeshFiles/unittest/cell1.gen.gz > cell1.gen
                      DEPENDS ${MOAB_SOURCE_DIR}/MeshFiles/unittest/cell1.gen.gz )
  add_custom_command( OUTPUT cell2.gen
                      COMMAND ${ZCAT_EXE} < ${MOAB_SOURCE_DIR}/MeshFiles/unittest/cell2.gen.gz > cell2.gen
                      DEPENDS ${MOAB_SOURCE_DIR}/MeshFiles/unittest/cell2.gen.gz)

  add_custom_target(moab_test_files DEPENDS mb_big_test.g cell1.gen cell2.gen)

  #MBTest.cpp:
  set_source_files_properties( MBTest.cpp
                               COMPILE_FLAGS "-DSRCDIR=${CMAKE_CURRENT_SOURCE_DIR} ${TEST_COMP_FLAGS}" )
  add_executable( moab_test MBTest.cpp )
  target_link_libraries( moab_test MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
  add_dependencies(moab_test moab_test_files)
  add_test( TestMOAB ${EXECUTABLE_OUTPUT_PATH}/moab_test )
endif()

add_executable( var_len_test_no_template VarLenTagTest.cpp )
set_source_files_properties( VarLenTagTest.cpp
  COMPILE_FLAGS "-DMOAB_TEMPLATE_SPECIALIZATION ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
target_link_libraries( var_len_test_no_template MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
add_test( var_len_test_no_template ${EXECUTABLE_OUTPUT_PATH}/var_len_test_no_template )

add_executable( mbcn_test mbcn_test.cpp )
set_target_properties( mbcn_test PROPERTIES
  COMPILE_FLAGS "-DMOAB_TEMPLATE_SPECIALIZATION ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
target_link_libraries( mbcn_test MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
add_test( mbcn_test ${EXECUTABLE_OUTPUT_PATH}/mbcn_test )

add_executable( file_options_test ${CMAKE_SOURCE_DIR}/src/FileOptions.cpp )
set_source_files_properties( ${CMAKE_SOURCE_DIR}/src/FileOptions.cpp
                             COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
target_link_libraries( file_options_test MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
add_test( file_options_test ${EXECUTABLE_OUTPUT_PATH}/file_options_test )

add_executable( homxform_test homxform_test.cpp )
set_source_files_properties( homxform_test.cpp
                             COMPILE_FLAGS "${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
target_link_libraries( homxform_test MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
add_test( TestHomXform ${EXECUTABLE_OUTPUT_PATH}/homxform_test )

add_executable( xform_test affinexform_test.cpp )
set_source_files_properties( affinexform_test.cpp
                             COMPILE_FLAGS "${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
target_link_libraries( xform_test MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
add_test( xform_test ${EXECUTABLE_OUTPUT_PATH}/xform_test )

#add imesh tests
if(ENABLE_IMESH)
  set( TESTS MOAB_iMesh_unit_tests.cpp
             MOAB_iMesh_extensions_tests.cpp )
  if(MOAB_HAVE_MPI)
    set(TESTS ${TESTS} MOAB_iMeshP_unit_tests.cpp )
  endif()
  foreach( fname ${TESTS} )
    string( REPLACE ".cpp" "" base ${fname} )
    add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${CMAKE_SOURCE_DIR}/itaps/imesh/${fname})
    set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "-DSRCDIR=${CMAKE_SOURCE_DIR}/itaps/imesh/ ${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
    target_link_libraries( ${base} iMesh MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
    add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
  endforeach()
endif()

if (ENABLE_IREL AND MOAB_HAVE_CGM)
  include_directories( ${CGM_INCLUDES} ${CMAKE_SOURCE_DIR}/itaps/irel ${CMAKE_BINARY_DIR}/itaps/irel )
  set( TESTS test_entset.cpp
             test_both.cpp
             test_status.cpp
             normals.cpp)
  set(irel_deplibs iRel)
  if(MOAB_HAVE_FBIGEOM)
    include_directories( ${CMAKE_SOURCE_DIR}/itaps/fbigeom ${CMAKE_BINARY_DIR}/itaps/fbigeom )
    set(TESTS ${TESTS} test_entsetFB.cpp )
    set(irel_deplibs ${irel_deplibs} FBiGeomMOAB)
  endif()
  set(irel_deplibs ${irel_deplibs} iMesh MOAB)
  foreach( fname ${TESTS} )
    string( REPLACE ".cpp" "" base ${fname} )
    add_executable( ${base} ${CMAKE_SOURCE_DIR}/test/TestUtil.hpp ${CMAKE_SOURCE_DIR}/itaps/irel/${fname})
    set_target_properties( ${base} PROPERTIES COMPILE_FLAGS "-DSRCDIR=${CMAKE_SOURCE_DIR}/itaps/irel/ ${TEST_COMP_FLAGS} ${MOAB_DEFINES} -DTEST" )
    target_link_libraries( ${base} ${irel_deplibs} ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
    add_test( ${base} ${EXECUTABLE_OUTPUT_PATH}/${base} )
  endforeach()
endif()

if (ENABLE_MESQUITE AND MOAB_HAVE_CGM)
  # add_subdirectory(mesquite)
endif(ENABLE_MESQUITE AND MOAB_HAVE_CGM)

if ( MOAB_HAVE_MPI AND MPI_FOUND )

  add_executable ( mbparallelcomm_test parallel/mbparallelcomm_test.cpp )
  target_link_libraries( mbparallelcomm_test MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
  set_source_files_properties( parallel/mbparallelcomm_test.cpp
    COMPILE_FLAGS "-DIS_BUILDING_MB ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
  add_test( TestParallelComm-BcastDelete
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} 0 ${CMAKE_SOURCE_DIR}/parallel/ptest.cub )
  add_test( TestParallelComm-ReadDelete
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} -1 ${CMAKE_SOURCE_DIR}/parallel/ptest.cub )
  add_test( TestParallelComm-ReadParallel
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} -2 ${CMAKE_SOURCE_DIR}/parallel/ptest.cub )
  add_test( TestParallelComm-Broadcast
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/mbparallelcomm_test ${MPIEXEC_POSTFLAGS} -3 ${CMAKE_SOURCE_DIR}/parallel/ptest.cub )

  if ( MOAB_HAVE_HDF5 )
    include_directories(
      ${HDF5_INCLUDES}
      ${CMAKE_SOURCE_DIR}/src/io/mhdf/include
    )
    set_source_files_properties( parallel/mhdf_parallel.c
      COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
    add_executable( mhdf_parallel parallel/mhdf_parallel.c )
    target_link_libraries( mhdf_parallel MOAB ${HDF5_LIBRARIES} ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
    add_test( TestMHDFParallel ${EXECUTABLE_OUTPUT_PATH}/mhdf_parallel )

    set_source_files_properties( parallel/parmerge_test.cpp
      COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
    add_executable( parmerge_test parallel/parmerge_test.cpp )
    target_link_libraries( parmerge_test MOAB ${HDF5_LIBRARIES} ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
    add_test( Testparmerge_test ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${EXECUTABLE_OUTPUT_PATH}/parmerge_test ${MPIEXEC_POSTFLAGS}  )

  endif ( MOAB_HAVE_HDF5 )

  set_source_files_properties( parallel/parallel_unit_tests.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
  add_executable ( parallel_unit_tests parallel/parallel_unit_tests.cpp )
  target_link_libraries( parallel_unit_tests MOAB ${CGM_LIBRARIES} ${MPI_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
  add_test( TestParallel
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} ${EXECUTABLE_OUTPUT_PATH}/parallel_unit_tests ${MPIEXEC_POSTFLAGS} )

  set_source_files_properties( parallel/pcomm_serial.cpp
    COMPILE_FLAGS "-DTEST ${MOAB_DEFINES} ${TEST_COMP_FLAGS}" )
  add_executable ( pcomm_serial parallel/pcomm_serial.cpp )
  target_link_libraries( pcomm_serial MOAB ${CGM_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )
  add_test( TestPCommSerial
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
    ${EXECUTABLE_OUTPUT_PATH}/pcomm_serial ${MPIEXEC_POSTFLAGS} 1 ${MOAB_SOURCE_DIR}/MeshFiles/unittest/ptest.cub )

endif ( MOAB_HAVE_MPI AND MPI_FOUND )
