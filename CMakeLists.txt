cmake_minimum_required( VERSION 3.5 )
SET(PACKAGE_NAME "MOAB")
SET(PACKAGE_VERSION "5.4.0")

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT WIN32)
  #This all breaks on windows.
  SET(CMAKE_Fortran_COMPILER_INIT ${CMAKE_GENERATOR_FC})
  SET(CMAKE_CXX_FLAGS_INIT "-fPIC -DPIC")
  SET(CMAKE_CXX_FLAGS_DEBUG_INIT "-O0 -g")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
  SET(CMAKE_CXX_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
  SET(CMAKE_C_FLAGS_INIT "-fPIC -DPIC")
  SET(CMAKE_C_FLAGS_DEBUG_INIT "-O0 -g")
  SET(CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
  SET(CMAKE_C_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
  SET(CMAKE_Fortran_FLAGS_INIT "-fPIC")
  SET(CMAKE_Fortran_FLAGS_DEBUG_INIT "-O0 -g")
  SET(CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")
  SET(CMAKE_Fortran_FLAGS_RELEASE_INIT "-O2")
  SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS_INIT "")
  SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS_INIT "")
endif()

project( MOAB )

include(GNUInstallDirs)
include(FeatureSummary)

#Add our Cmake directory to the module search path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/config ${CMAKE_MODULE_PATH})
################################################################################
# Set up version info
################################################################################
include (config/GetAcInitVersion.cmake)
get_ac_init_version()
set ( MOAB_VERSION_STRING "${PACKAGE_NAME} ${VERSION_STRING}" )
set ( MOAB_VERSION        "${VERSION_STRING}" )
set ( MOAB_VERSION_MAJOR  ${MAJOR_VERSION}  )
set ( MOAB_VERSION_MINOR  ${MINOR_VERSION}  )
if ( DEFINED PATCH_VERSION )
  set ( MOAB_VERSION_PATCH "${PATCH_VERSION}" )
else ( DEFINED PATCH_VERSION )
  if ( MOAB_VERSION_MINOR EQUAL 99 )
    set ( MOAB_VERSION_STRING "${MOAB_VERSION_STRING} (alpha)" )
  else ( MOAB_VERSION_MINOR EQUAL 99 )
    set ( MOAB_VERSION_STRING "${MOAB_VERSION_STRING} (beta)" )
  endif ( MOAB_VERSION_MINOR EQUAL 99 )
endif ( DEFINED PATCH_VERSION )

IF(CMAKE_VERSION VERSION_EQUAL "3.0.0" OR CMAKE_VERSION VERSION_GREATER "3.0.0")
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0020 NEW)
  #cmake_policy(SET CMP0042 OLD)
ENDIF ()
if (NOT (CMAKE_VERSION VERSION_LESS "3.12.0"))
# Allow use of <PackageName>_ROOT as a cmake variable
  cmake_policy(SET CMP0074 NEW)
endif ()


set ( abs_srcdir   ${CMAKE_SOURCE_DIR} )
set ( abs_builddir ${PROJECT_BINARY_DIR} )

################################################################################
# Install Related Settings
################################################################################

get_filename_component( MOAB_ABSSRC_DIR2 moab.make.in REALPATH )
execute_process( COMMAND dirname ${MOAB_ABSSRC_DIR2}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE MOAB_ABSSRC_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

## Set the directory where the binaries will be stored
set( EXECUTABLE_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/bin
  CACHE PATH
  "Directory where all executable will be stored"
)

## Set the directory where the libraries will be stored
set( LIBRARY_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/lib
  CACHE PATH
  "Directory where all the libraries will be stored"
)
mark_as_advanced(
  MOAB_ABSSRC_DIR
  EXECUTABLE_OUTPUT_PATH
  LIBRARY_OUTPUT_PATH)

include ( CheckIncludeFile )
include ( CheckFunctionExists )
include ( CheckTypeSize )

# Compiler defines... this should really be in a config file.
set( MOAB_DEFINES "" )
set( MOAB_LIBS )
set( MOAB_INSTALL_TARGETS )
set( MOAB_MESH_DIR        "${CMAKE_SOURCE_DIR}/MeshFiles/" )

################################################################################
# Options that the user controls
################################################################################
option ( BUILD_SHARED_LIBS   "Should shared or static libraries be created?"   ON  )
option ( ENABLE_CGM        "Should build with CGM support?"                  OFF )
option ( ENABLE_CGNS       "Should build with CGNS support?"                 OFF )
option ( ENABLE_MPI        "Should MOAB be compiled with MPI support?"       OFF )
option ( ENABLE_HDF5       "Include HDF I/O interfaces in the build?"                   OFF )
option ( ENABLE_ZLIB       "Include Zlib compression libraries (optionally used in HDF5)?"  OFF )
option ( ENABLE_SZIP       "Should build with szip support?"                 OFF )
option ( ENABLE_NETCDF     "Include NetCDF (ExodusII) interfaces in the build?" OFF )
option ( ENABLE_PNETCDF    "Include PNetCDF interfaces in the build (Requires NetCDF) ?" OFF )
option ( ENABLE_METIS      "Include Metis support for partitioning algorithms?" OFF )
option ( ENABLE_PARMETIS   "Include ParMetis support for partitioning algorithms?" OFF )
option ( ENABLE_ZOLTAN     "Include Zoltan support for partitioning algorithms?" OFF )
option ( ENABLE_VTK        "Include VTK I/O interfaces in the build?" OFF )
option ( ENABLE_TESTING    "Enable Testing"                                  ON  )
option ( MOAB_FORCE_64_BIT_HANDLES "Force MBEntityHandle to be 64 bits (uint64_t)" OFF )
option ( MOAB_FORCE_32_BIT_HANDLES "Force MBEntityHandle to be 32 bits (uint32_t)" OFF )
option ( ENABLE_PYMOAB       "Build python interface to MOAB?" OFF )
option ( ENABLE_IMESH        "Should build IMESH?"       ON )
option ( ENABLE_FBIGEOM      "Should build FBIGEOM?"     OFF )
option ( ENABLE_IREL         "Should build IREL?"        OFF )
option ( ENABLE_MESQUITE     "Should build Mesquite mesh optimization interfaces?" OFF )
option ( ENABLE_BLASLAPACK   "Should require BLAS/LAPACK?" ON )
option ( ENABLE_FORTRAN      "Should enable FORTRAN?" ON )
option ( ENABLE_TEMPESTREMAP "Should enable TempestRemap dependency?" OFF )

add_feature_info(SharedLibs BUILD_SHARED_LIBS "Should shared or static libraries be created?")
add_feature_info(EnableFortran ENABLE_FORTRAN "Should FORTRAN support be enabled?")

mark_as_advanced(
  MOAB_FORCE_64_BIT_HANDLES
  MOAB_FORCE_32_BIT_HANDLES
)

if ( ENABLE_IMESH OR ENABLE_FBIGEOM )
  include (${CMAKE_ROOT}/Modules/CMakeDetermineFortranCompiler.cmake)
endif ()

# patch added to ensure libs have a full-path ID set on Mac
if(APPLE)
  set(CMAKE_MACOSX_RPATH FALSE)
  set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

################################################################################
# Check for system include files
################################################################################
check_include_file( inttypes.h   MOAB_HAVE_INTTYPES_H )
check_include_file( stdint.h     MOAB_HAVE_STDINT_H )
check_include_file( stddef.h     MOAB_HAVE_STDDEF_H )
check_include_file( stdlib.h     MOAB_HAVE_STDLIB_H )
check_include_file( sys/types.h  MOAB_HAVE_SYS_TYPES_H )
set( HAVE_SYS_TYPES_H ${MOAB_HAVE_SYS_TYPES_H} )
set( HAVE_STDDEF_H    ${MOAB_HAVE_STDDEF_H} )
set( HAVE_STDINT_H    ${MOAB_HAVE_STDINT_H} )
set( HAVE_INTTYPES_H    ${MOAB_HAVE_INTTYPES_H} )
set( HAVE_STDLIB_H    ${MOAB_HAVE_STDLIB_H} )
check_include_file( memory.h     MOAB_HAVE_MEMORY_H )

INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(WORDS_BIGENDIAN)
################################################################################
# Integer size Related Settings
################################################################################
if ( MOAB_FORCE_64_BIT_HANDLES AND MOAB_FORCE_32_BIT_HANDLES )
  message( FATAL_ERROR
      "You may not turn both MOAB_FORCE_64_BIT_HANDLES and MOAB_FORCE_32_BIT_HANDLES on. Turn one off to continue."
    )
endif ( MOAB_FORCE_64_BIT_HANDLES AND MOAB_FORCE_32_BIT_HANDLES )

if ( NOT MOAB_FORCE_64_BIT_HANDLES AND NOT MOAB_FORCE_32_BIT_HANDLES )
  if ( MOAB_HAVE_INTTYPES_H )
    set ( CMAKE_EXTRA_INCLUDE_FILES "${CMAKE_EXTRA_INCLUDE_FILES};inttypes.h" )
  endif ( MOAB_HAVE_INTTYPES_H )
  if ( MOAB_HAVE_STDLIB_H )
    set ( CMAKE_EXTRA_INCLUDE_FILES "${CMAKE_EXTRA_INCLUDE_FILES};stdlib.h" )
    #set ( CHECK_TYPE_SIZE_PREMAIN "${CHECK_TYPE_SIZE_PREMAIN}\n#include <stdlib.h>\n" )
  endif ( MOAB_HAVE_STDLIB_H )
  check_type_size(  size_t       MOAB_HAVE_SIZE_T )
  check_type_size(  ptrdiff_t    MOAB_HAVE_PTRDIFF_T )
  set ( HAVE_SIZE_T ${MOAB_HAVE_SIZE_T} )
  set ( HAVE_PTRDIFF_T ${MOAB_HAVE_PTRDIFF_T} )
endif ( NOT MOAB_FORCE_64_BIT_HANDLES AND NOT MOAB_FORCE_32_BIT_HANDLES )

# Check of presence of finite or isfinite macros
include (CheckFunctionExists)
CHECK_FUNCTION_EXISTS(std::isfinite MOAB_HAVE_STDISFINITE)
CHECK_FUNCTION_EXISTS(isfinite MOAB_HAVE_ISFINITE)
CHECK_FUNCTION_EXISTS(finite MOAB_HAVE_FINITE)

################################################################################
# Find packages
################################################################################

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

include(config/CMakeLibraryFuncs.cmake)
include(GenerateExportHeader)
if (BUILD_SHARED_LIBS)
  if (WIN32)
    add_definitions(/DMOAB_DLL)
    if (NOT DEFINED CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS)
# https://cmake.org/cmake/help/v3.4/prop_tgt/WINDOWS_EXPORT_ALL_SYMBOLS.html
# Export all symbols
      set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_USE_MATH_DEFINES")
    endif ()
  else ()
    # add_definitions(-DMOAB_DLL) # Not needed?
  endif ()
endif ()

set (MOAB_HAVE_BLAS OFF CACHE INTERNAL "Configured MOAB with BLAS library support." )
set (MOAB_HAVE_LAPACK OFF CACHE INTERNAL "Configured MOAB with LAPACK library support." )

if ( ENABLE_BLASLAPACK  )
  # Find BLAS
  FIND_PACKAGE(BLAS REQUIRED)
  if (BLAS_FOUND)
    set(MOAB_HAVE_BLAS ON)
  endif(BLAS_FOUND)
  set_package_properties(BLAS PROPERTIES
                       TYPE RECOMMENDED
                       PURPOSE "Native BLAS implementations")

  # Find LAPACK
  FIND_PACKAGE(LAPACK REQUIRED)
  if (LAPACK_FOUND)
    set(MOAB_HAVE_LAPACK ON)
  endif(LAPACK_FOUND)

  set_package_properties(LAPACK PROPERTIES
                       TYPE RECOMMENDED
                       PURPOSE "Native LAPACK implementations")

endif (ENABLE_BLASLAPACK)

if ( LAPACK_FOUND )

  SET (ENABLE_FORTRAN ON)

  # Eigen3 is optional if LAPACK has been found
  FIND_PACKAGE(Eigen3 OPTIONAL_COMPONENTS)

  set_package_properties(Eigen3 PROPERTIES
                       TYPE RECOMMENDED
                       PURPOSE "Eigen3 headers as supplement to BLAS/LAPACK implementations")

else ( LAPACK_FOUND )

  # Find Eigen3 is required if LAPACK is not found (or requested)
  FIND_PACKAGE(Eigen3 REQUIRED)

  set_package_properties(Eigen3 PROPERTIES
                        TYPE REQUIRED
                       PURPOSE "Eigen3 headers as replacement for BLAS/LAPACK implementations")

endif ( LAPACK_FOUND )

set (MOAB_HAVE_EIGEN3 OFF CACHE INTERNAL "Configured MOAB with Eigen3 library support." )
if ( EIGEN3_FOUND )
  INCLUDE_DIRECTORIES ( "${EIGEN3_DIR}" )
  SET(MOAB_DEFINES "${MOAB_DEFINES} ${EIGEN3_INCLUDES}")
  SET(MOAB_HAVE_EIGEN3 ON)
endif (EIGEN3_FOUND)

if (ENABLE_FORTRAN)
  enable_language(Fortran)

  # Include CMake's Fortran mangling evaluation macros
  include(FortranCInterface)

  if ("${FortranCInterface_GLOBAL_CASE}" MATCHES "UPPER")
    set(FC_CASE "NAME")
  else ("${FortranCInterface_GLOBAL_CASE}" MATCHES "UPPER")
    set(FC_CASE "name")
  endif ("${FortranCInterface_GLOBAL_CASE}" MATCHES "UPPER")
  if (NOT "${FortranCInterface_GLOBAL_SUFFIX}" STREQUAL "")
    set(FortranCInterface_GLOBAL_SUFFIX "## ${FortranCInterface_GLOBAL_SUFFIX}")
  endif (NOT "${FortranCInterface_GLOBAL_SUFFIX}" STREQUAL "")
  if (NOT "${FortranCInterface_GLOBAL__SUFFIX}" STREQUAL "")
    set(FortranCInterface_GLOBAL__SUFFIX "## ${FortranCInterface_GLOBAL__SUFFIX}")
  endif (NOT "${FortranCInterface_GLOBAL__SUFFIX}" STREQUAL "")

  set(MOAB_F77_FUNC  "${FortranCInterface_GLOBAL_PREFIX} ${FC_CASE} ${FortranCInterface_GLOBAL_SUFFIX}")
  set(MOAB_F77_FUNC_ "${FortranCInterface_GLOBAL__PREFIX} ${FC_CASE} ${FortranCInterface_GLOBAL__SUFFIX}")
  set(MOAB_FC_FUNC   "${FortranCInterface_GLOBAL_PREFIX} ${FC_CASE} ${FortranCInterface_GLOBAL_SUFFIX}")
  set(MOAB_FC_FUNC_  "${FortranCInterface_GLOBAL__PREFIX} ${FC_CASE} ${FortranCInterface_GLOBAL__SUFFIX}")
endif (ENABLE_FORTRAN)

# check for MPI package
set (MOAB_HAVE_MPI OFF CACHE INTERNAL "Found necessary MPI components. Configure MOAB with it." )
if ( ENABLE_MPI )
  set (MPI_C_COMPILER ${CC} )
  set (MPI_CXX_COMPILER ${CXX} )
  set (MPI_Fortran_COMPILER ${FC} )
  find_package( MPI REQUIRED )
  # CMake FindMPI script is sorely lacking:
  if (MPI_C_COMPILER AND (NOT MPI_INCLUDE_PATH OR NOT MPI_LIBRARY))
    get_filename_component(MPI_BINDIR ${MPI_C_COMPILER} DIRECTORY)
    get_filename_component(MPI_ROOTDIR ${MPI_BINDIR} DIRECTORY)
  endif ()
  if (MPI_ROOTDIR AND (NOT MPI_INCLUDE_PATH) AND EXISTS ${MPI_ROOTDIR}/include)
    set(MPI_INCLUDE_PATH ${MPI_ROOTDIR}/include)
  endif ()
  if (MPI_ROOTDIR AND (NOT MPI_LIBRARY))
    find_library(MPI_LIBRARY mpi PATHS ${MPI_ROOTDIR} PATH_SUFFIXES lib64 lib)
  endif ()
  if ( MPI_LIBRARY AND MPI_INCLUDE_PATH )
    set( MPI_FOUND ON )
  endif ()

  if ( MPI_FOUND )
    set ( MOAB_DEFINES "${MOAB_DEFINES} -DUSE_MPI" )
    set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
    include_directories(${MPI_INCLUDE_PATH})
    set( MOAB_HAVE_MPI ON )
    message(STATUS "MPI libraries and includes found..")
    foreach (var C_COMPILER CXX_COMPILER BINDIR ROOTDIR INCLUDE_PATH LIBRARY)
      message(STATUS "    [VAR]: MPI_${var} = ${MPI_${var}}.")
    endforeach ()

    set_package_properties(MPI PROPERTIES
                       TYPE RECOMMENDED
                       PURPOSE "MPI library to enable node-level parallelism")

  endif ( MPI_FOUND )
endif ( ENABLE_MPI )

include (config/CheckCompilerFlags.cmake)

set (MOAB_HAVE_ZLIB OFF CACHE INTERNAL "Found necessary Zlib components. Configure MOAB with it." )
if ( ENABLE_ZLIB )
  find_package( ZLIB REQUIRED )
  set (MOAB_HAVE_ZLIB ON)
endif (ENABLE_ZLIB)

set (MOAB_HAVE_SZIP OFF CACHE INTERNAL "Found necessary Zlib components. Configure MOAB with it." )
if ( ENABLE_SZIP )
  find_package( SZIP REQUIRED )
  set (MOAB_HAVE_SZIP ON)
endif (ENABLE_SZIP)

find_package( CURL )

set (MOAB_HAVE_HDF5 OFF CACHE INTERNAL "Found necessary HDF5 components. Configure MOAB with it." )
set (MOAB_HAVE_HDF5_PARALLEL OFF CACHE INTERNAL "Found necessary parallel HDF5 components. Configure MOAB with it." )
if ( ENABLE_HDF5 )
  if (MOAB_HAVE_MPI)
    set(HDF5_PREFER_PARALLEL TRUE)
  endif()

  #TODO: test on all architectures with CMake-native FindHDF5 module instead of a custom one
  IF(WIN32)
    find_package( HDF5 REQUIRED )
  ELSE(WIN32)
    find_package( HDF5_MOAB REQUIRED )
  ENDIF(WIN32)

  set (MOAB_HAVE_HDF5 ON)
  # if ( HDF5_IS_PARALLEL AND MOAB_HAVE_MPI )
  if ( MOAB_HAVE_MPI )
    set (MOAB_HAVE_HDF5_PARALLEL ON)
  endif ( )
  if ( CURL_FOUND )
    set( HDF5_LIBRARIES ${HDF5_LIBRARIES} ${CURL_LIBRARIES} )
  endif ( CURL_FOUND )
  set( MOAB_LIBS ${HDF5_LIBRARIES} ${MOAB_LIBS} )
  include_directories( ${HDF5_INCLUDE_DIR} src/io/mhdf/include )
endif ( ENABLE_HDF5 )
set_package_properties(HDF5 PROPERTIES
                    TYPE RECOMMENDED
                    PURPOSE "HDF5 library interfaces for MOAB native file format support")

# Check for NetCDF availability
set (MOAB_HAVE_NETCDF OFF CACHE INTERNAL "Found necessary NetCDF components. Configure MOAB with it." )
if ( ENABLE_NETCDF )
  find_package( NetCDF REQUIRED )
  # set( MOAB_DEFINES "-DNETCDF_FILE ${MOAB_DEFINES}" )
  set (MOAB_HAVE_NETCDF ON)
  set (MOAB_HAVE_NETCDF_H 1)
  include_directories( ${NETCDF_INCLUDES} )
  set( MOAB_LIBS ${NETCDF_LIBRARIES} ${MOAB_LIBS} )
endif ( ENABLE_NETCDF )
set_package_properties(NetCDF PROPERTIES
                    TYPE RECOMMENDED
                    PURPOSE "NetCDF library interfaces for MOAB ExodusII file format support")

# Check for PNetCDF availability
set (MOAB_HAVE_PNETCDF OFF CACHE INTERNAL "Found necessary PNetCDF components. Configure MOAB with it." )
if ( ENABLE_PNETCDF AND MOAB_HAVE_NETCDF )
	find_package( PNetCDF REQUIRED )
	set (MOAB_HAVE_PNETCDF ON)
	include_directories( ${PNETCDF_INCLUDES} )
	set( MOAB_LIBS ${PNETCDF_LIBRARIES} ${MOAB_LIBS} )
endif ( ENABLE_PNETCDF AND MOAB_HAVE_NETCDF )
set_package_properties(PNetCDF PROPERTIES
                    TYPE RECOMMENDED
                    PURPOSE "Parallel-NetCDF library interfaces for MOAB nc-file format support")

# Check for Metis partitioner availability
set (MOAB_HAVE_METIS OFF CACHE INTERNAL "Found necessary Metis components. Configure MOAB with it." )
if ( ENABLE_METIS )
  find_package( Metis REQUIRED )
  set (MOAB_HAVE_METIS ON)
endif (ENABLE_METIS )
set_package_properties(Metis PROPERTIES
                    TYPE RECOMMENDED
                    PURPOSE "Metis library interfaces for MOAB offline (serial) mesh partitioning support")

# Check for ParMetis partitioner availability
set (MOAB_HAVE_PARMETIS OFF CACHE INTERNAL "Found necessary ParMetis components. Configure MOAB with it." )
if ( ENABLE_PARMETIS AND ENABLE_METIS AND METIS_FOUND )
  find_package( ParMetis REQUIRED )
  set (MOAB_HAVE_PARMETIS ON)
endif ( ENABLE_PARMETIS AND ENABLE_METIS AND METIS_FOUND )
set_package_properties(ParMetis PROPERTIES
                    TYPE RECOMMENDED
                    PURPOSE "ParMetis library interfaces for MOAB parallel mesh partitioning support (through Zoltan)")

# Check for Zoltan partitioner availability
set (MOAB_HAVE_ZOLTAN OFF CACHE INTERNAL "Found necessary Zoltan components. Configure MOAB with it." )
if ( ENABLE_ZOLTAN )
  find_package( Zoltan REQUIRED )
  set (MOAB_HAVE_ZOLTAN ON)
endif (ENABLE_ZOLTAN )
set_package_properties(Zoltan PROPERTIES
                    TYPE RECOMMENDED
                    PURPOSE "Zoltan library interfaces for MOAB offline/online parallel mesh partitioning support")

# Check for TempestRemap availability if Eigen3 headers have been found already
if ( EIGEN3_FOUND AND ENABLE_TEMPESTREMAP )
  find_package( TempestRemap REQUIRED )
  set(MOAB_HAVE_TEMPESTREMAP ON)
  set_package_properties(TempestRemap PROPERTIES
                    TYPE OPTIONAL
                    PURPOSE "TempestRemap library interfaces for MOAB climate remapping support")

endif( EIGEN3_FOUND AND ENABLE_TEMPESTREMAP )

# Check if user has asked for configuration with CGM
set (MOAB_HAVE_CGM OFF CACHE INTERNAL "Found necessary CGM components. Configure MOAB with it." )
if ( ENABLE_CGM )
   find_package( CGM REQUIRED )
   set( MOAB_DEFINES "${CGM_DEFINES} -DCGM ${MOAB_DEFINES}" )
   set (MOAB_HAVE_CGM ON)
   set (ENABLE_IREL ON)
   set (ENABLE_FBIGEOM ON)
   set_package_properties(CGM PROPERTIES
                    TYPE OPTIONAL
                    PURPOSE "CGM library interfaces for MOAB native-geometry support")
endif (ENABLE_CGM)

set (MOAB_HAVE_CGNS OFF CACHE INTERNAL "Found necessary CGNS components. Configure MOAB with it." )
if (ENABLE_CGNS)
  set( MOABIO_LIBS ${MOABIO_LIBS} ${CGNS_LIBRARIES} )
  set (MOAB_HAVE_CGNS ON)
endif(ENABLE_CGNS)

# Strip our custom flags
string(STRIP "${CMAKE_EXE_LINKER_FLAGS} " CMAKE_EXE_LINKER_FLAGS)
string(STRIP "${MOAB_CXX_FLAGS} " MOAB_CXX_FLAGS)
string(STRIP "${MOAB_DEFINES} " MOAB_DEFINES)

################################################################################
# Add Directories
################################################################################
message(STATUS "Adding src subdirectory.")
add_subdirectory( src )
message(STATUS "Adding itaps subdirectory.")
add_subdirectory( itaps )
if (NOT WIN32 OR MSYS OR MINGW) # Need further work to prepare for windows
  add_subdirectory( tools )
endif ()
if (ENABLE_PYMOAB)
  add_subdirectory( pymoab )
endif()
################################################################################
# Testing Related Settings
################################################################################
#turn on ctest if we want testing
if ( ENABLE_TESTING )
  enable_testing()
  add_subdirectory( test )
endif()

###############################################################################
#
###############################################################################
export(
  TARGETS ${MOAB_INSTALL_TARGETS}
  FILE "${PROJECT_BINARY_DIR}/MOABTargets.cmake"
  )
install(
  EXPORT MOABTargets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOAB/
  )

# Create some custom command

# For consistency with autoconf, create a "make check" equivalent to "ctest"
ADD_CUSTOM_TARGET( check COMMAND ctest DEPENDS )

# For consistency with autoconf, create "make dist" and "make distcheck" targets
INCLUDE( CMakeDistCheck )
DEFINE_DISTCHECK ( )

################################################################################
# Generate the MOABConfig.cmake file
################################################################################
set(CXX ${CMAKE_CXX_COMPILER})
set(CC ${CMAKE_C_COMPILER})
if (ENABLE_FORTRAN)
  set(F77 ${CMAKE_Fortran_COMPILER} )
  set(FC ${CMAKE_Fortran_COMPILER})
endif (ENABLE_FORTRAN)

set(CPPFLAGS "${CMAKE_CXX_FLAGS} ${MOAB_DEFINES}")
IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  set(CXXFLAGS ${CMAKE_CXX_FLAGS_DEBUG})
  set(CFLAGS ${CMAKE_C_FLAGS_DEBUG})
  if (ENABLE_FORTRAN)
    set(FFLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
  endif (ENABLE_FORTRAN)
ELSE (CMAKE_BUILD_TYPE MATCHES "Debug")
  set(CXXFLAGS ${CMAKE_CXX_FLAGS_RELEASE})
  set(CFLAGS ${CMAKE_C_FLAGS_RELEASE})
  if (ENABLE_FORTRAN)
    set(FFLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
  endif (ENABLE_FORTRAN)
ENDIF (CMAKE_BUILD_TYPE MATCHES "Debug")

# CMake does not define uninstall target. Look at install manifest and remove
# all the files recursively.
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/config/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
    )

## Find and use clang-format if available
find_program( CLANG_FORMAT_EXE NAMES "clang-format" DOC "Path to clang-format executable" )
if(CLANG_FORMAT_EXE)
	message(STATUS "clang-format found: ${CLANG_FORMAT_EXE}")
	set(DO_CLANG_FORMAT "${CLANG_FORMAT_EXE}" "-i -style=file")
endif()

# get all project files
file(GLOB_RECURSE MOAB_PROJECT_SOURCE_FILES
			RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
			LIST_DIRECTORIES false
			src/*.cpp src/*.c src/*.hpp
			tools/*.cpp tools/*.hpp
			test/*.cpp
			itaps/*.cpp itaps/*.h itaps/*.hpp
			examples/*.cpp
		)

function(prepend var prefix)
  set(listVar "")
  foreach(f ${ARGN})
    list(APPEND listVar "${prefix}/${f}")
  endforeach()
  set(${var} "${listVar}" PARENT_SCOPE)
endfunction()

if(CLANG_FORMAT_EXE)
	prepend(FILES_TO_FORMAT ${CMAKE_CURRENT_SOURCE_DIR} ${MOAB_PROJECT_SOURCE_FILES})
	add_custom_target( clang-format-project-files COMMAND ${CLANG_FORMAT_EXE} -i -style=file ${FILES_TO_FORMAT} )
endif()

## Check for clang-tidy as well for performing pre-commit checks
find_program(CLANG_TIDY_EXE NAMES "clang-tidy" DOC "Path to clang-tidy executable")
if(CLANG_TIDY_EXE)
  message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
  SET(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_EXE}
    -format-style=file
    -checks=clan*,mpi-*,port*,bug*,llvm-namespace-comment,perf*,modernize-deprecated-headers,modernize-avoid-bind,modernize-use-bool-literals,misc-noexcept-move-constructor,-readability-braces-around-statements,-misc-unused-parameters
    -header-filter="\"moab/*\""
    -fix -fix-errors
    -p=${CMAKE_CURRENT_BINARY_DIR} )

  find_program(RUN_CLANG_TIDY_BIN /usr/local/Cellar/llvm/10.0.0_3/share/clang/run-clang-tidy.py)
  find_program(RUN_CLANG_APPLY_REPLACEMENTS_BIN clang-apply-replacements PATH /usr/local/bin /usr/local/opt/llvm/bin)
  list(APPEND RUN_CLANG_TIDY_BIN_ARGS
        -clang-tidy-binary ${CLANG_TIDY_EXE}
        -clang-apply-replacements-binary ${RUN_CLANG_APPLY_REPLACEMENTS_BIN}
        -header-filter="\"moab/*\""
        -checks=clan*,mpi-*,port*,bug*,llvm-namespace-comment,perf*,modernize-deprecated-headers,modernize-avoid-bind,modernize-use-bool-literals,misc-noexcept-move-constructor,-readability-braces-around-statements,-misc-unused-parameters
        -fix
        )

	add_custom_target(
		tidy
		COMMAND ${RUN_CLANG_TIDY_BIN} ${RUN_CLANG_TIDY_BIN_ARGS}
		COMMENT "running clang tidy"
		)
	add_custom_target( clang-tidy-project-files COMMAND ${CMAKE_CXX_CLANG_TIDY} ${FILES_TO_FORMAT} )

endif(CLANG_TIDY_EXE)

# Generate the CMake input header file and then pass it through configuration
## execute_process( COMMAND cat MOABConfig.h.in COMMAND sed -e "s/#undef /#cmakedefine MOAB_/g"
#execute_process( COMMAND cat MOABConfig.h.in COMMAND sed -r "s/#undef ([A-Z].*)/#cmakedefine MOAB_\\1 @MOAB_\\1@/"
#  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/config
#  OUTPUT_FILE MOABConfig.cmake2.h.in
#  OUTPUT_STRIP_TRAILING_WHITESPACE
#  )
CONFIGURE_FILE(config/MOABConfig.cmake.h.in "${PROJECT_BINARY_DIR}/src/moab/MOABConfig.h" @ONLY)

CONFIGURE_FILE(config/MOABConfig.cmake.cmake "${PROJECT_BINARY_DIR}/MOABConfig.cmake" @ONLY)
CONFIGURE_FILE(moab.make.cmake "${PROJECT_BINARY_DIR}/lib/moab.make" @ONLY)
if ( MOAB_HAVE_FBIGEOM )
  CONFIGURE_FILE(itaps/fbigeom/FBiGeom-Defs.inc.in "${PROJECT_BINARY_DIR}/lib/FBiGeom-Defs.inc" @ONLY)
  install(
    FILES "${PROJECT_BINARY_DIR}/lib/FBiGeom-Defs.inc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( MOAB_HAVE_FBIGEOM )

if ( MOAB_HAVE_IMESH )
  CONFIGURE_FILE(itaps/imesh/iMesh-Defs.inc.in "${PROJECT_BINARY_DIR}/lib/iMesh-Defs.inc" @ONLY)
  INSTALL(
    FILES "${PROJECT_BINARY_DIR}/lib/iMesh-Defs.inc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ( MOAB_HAVE_IMESH )
INSTALL(
  FILES "${PROJECT_BINARY_DIR}/lib/moab.make"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
INSTALL(
  CODE
  "file(APPEND  ${CMAKE_INSTALL_FULL_LIBDIR}/moab.make MOAB_LIBDIR=${CMAKE_INSTALL_FULL_LIBDIR}\\n )"
  )
INSTALL(
  CODE
  "file(APPEND  ${CMAKE_INSTALL_FULL_LIBDIR}/moab.make MOAB_INCLUDES=-I${CMAKE_INSTALL_FULL_INCLUDEDIR}\\n )"
  )
INSTALL(
  FILES "${PROJECT_BINARY_DIR}/MOABConfig.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOAB/
  )
INSTALL(
  FILES "${PROJECT_SOURCE_DIR}/config/ResolveCompilerPaths.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOAB
  )

ADD_CUSTOM_TARGET(
  uninstall
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
  )

CONFIGURE_FILE(config/MOAB.pc.cmake.in "${PROJECT_BINARY_DIR}/MOAB.pc" @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/MOAB.pc
  DESTINATION lib/pkgconfig)

###################
# INSTALL documentation and examples
####################

SET(MOAB_DOCS ${MOAB_DOCS}
  LICENSE README.md RELEASE_NOTES.md TODO.md
)
install(FILES ${MOAB_DOCS} DESTINATION share/doc/moab)

CONFIGURE_FILE(examples/makefile.config.cmake "${PROJECT_BINARY_DIR}/share/examples/makefile.config" @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/share/examples/makefile.config
  DESTINATION share/examples)

SET(MOAB_EXAMPLES
  examples/makefile
  examples/makefile.config.in
  examples/CMakeLists.txt )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/beginner/README
  examples/beginner/makefile
  examples/beginner/mbex1.cpp
  examples/beginner/mbex2.cpp
  examples/beginner/mbex3.cpp
  examples/beginner/mbex4.cpp )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/basic/makefile
  examples/basic/ErrorHandlingSimulation.cpp
  examples/basic/GetEntities.cpp
  examples/basic/HelloMOAB.cpp
  examples/basic/HelloParMOAB.cpp
  examples/basic/ReadWriteTest.cpp
  examples/basic/SetsNTags.cpp
  examples/basic/StructuredMeshSimple.cpp
  examples/basic/TestErrorHandling.cpp
  examples/basic/TestErrorHandlingPar.cpp
  examples/basic/TestExodusII.cpp
  examples/basic/UniformRefinement.cpp )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/intermediate/makefile
  examples/intermediate/DirectAccessNoHoles.cpp
  examples/intermediate/DirectAccessWithHoles.cpp
  examples/intermediate/LoadPartial.cpp
  examples/intermediate/PointInElementSearch.cpp
  examples/intermediate/ReduceExchangeTags.cpp )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/advanced/makefile
  examples/advanced/CrystalRouterExample.cpp
  examples/advanced/ExtrudePoly.cpp
  examples/advanced/GenLargeMesh.cpp
  examples/advanced/LloydRelaxation.cpp
  examples/advanced/VisTags.cpp
  examples/advanced/CoupleMGen.cpp )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/fortran/makefile
  examples/fortran/DirectAccessNoHolesF90.F90
  examples/fortran/PushParMeshIntoMoabF90.F90
  examples/fortran/MigrateMesh.F90 )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/python/makefile
  examples/python/laplaciansmoother.py )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/apps/makefile )
SET(MOAB_EXAMPLES ${MOAB_EXAMPLES}
  examples/apps/earth/makefile
  examples/apps/earth/ContinentsOnGlobe.cpp
  examples/apps/earth/BoundaryContinents.cpp
  examples/apps/earth/addPCdata.cpp )

include(AutoconfHeader)
moab_install_docs(${MOAB_EXAMPLES})


feature_summary(WHAT ALL #ENABLED_FEATURES DISABLED_FEATURES PACKAGES_FOUND
                DESCRIPTION "MOAB Configuration Summary"
                INCLUDE_QUIET_PACKAGES
                FILENAME ${CMAKE_BINARY_DIR}/summary.log
                )
MESSAGE("-- MOAB Configuration summary has been written to ${CMAKE_BINARY_DIR}/summary.log.")
# All done.
