# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign subdir-objects

# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I config

# Subdirectories to build
SUBDIRS = src MeshFiles $(itaps_dir) tools test $(pymoab_dir)

if ENABLE_fbigeom
  itaps_dir_fbigeom = itaps
else
  itaps_dir_fbigeom =
endif

if ENABLE_imesh
  itaps_dir = itaps
else
  itaps_dir = $(itaps_dir_fbigeom)
endif

if ENABLE_pymoab
  pymoab_dir = pymoab
else
  pymoab_dir =
endif

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check

doc_DATA = .gitignore AUTHORS LICENSE README.md RELEASE_NOTES.md TODO.md CODE_OF_CONDUCT.md CONTRIBUTING.md FAQ.md

pkgconfig_DATA = MOAB.pc
pkgconfigdir = $(libdir)/pkgconfig

examples_DATA = \
               examples/makefile \
               examples/makefile.config \
               examples/makefile.config.in \
               examples/CMakeLists.txt \
               examples/makefile.config.cmake
examples_cmakedir = $(datadir)/examples/cmake
examples_cmake_DATA = examples/cmake/CMakeLists.txt
examples_beginnerdir = $(datadir)/examples/beginner
examples_beginner_DATA = \
               examples/beginner/README \
               examples/beginner/makefile \
               examples/beginner/mbex1.cpp \
               examples/beginner/mbex2.cpp \
               examples/beginner/mbex3.cpp \
               examples/beginner/mbex4.cpp
examples_basicdir = $(datadir)/examples/basic
examples_basic_DATA = \
               examples/basic/makefile \
               examples/basic/ErrorHandlingSimulation.cpp \
               examples/basic/DeleteEdges.cpp \
               examples/basic/GetEntities.cpp \
               examples/basic/HelloMOAB.cpp \
               examples/basic/HelloParMOAB.cpp \
               examples/basic/ReadWriteTest.cpp \
               examples/basic/SetsNTags.cpp \
               examples/basic/StructuredMeshSimple.cpp \
               examples/basic/TestErrorHandling.cpp \
               examples/basic/TestErrorHandlingPar.cpp \
               examples/basic/TestExodusII.cpp \
               examples/basic/ReadPartFile.cpp \
               examples/basic/UniformRefinement.cpp
examples_intermediatedir = $(datadir)/examples/intermediate
examples_intermediate_DATA = \
               examples/intermediate/makefile \
               examples/intermediate/DirectAccessNoHoles.cpp \
               examples/intermediate/DirectAccessWithHoles.cpp \
               examples/intermediate/LoadPartial.cpp \
               examples/intermediate/PointInElementSearch.cpp \
               examples/intermediate/ReduceExchangeTags.cpp \
               examples/intermediate/QuadTriConv.cpp \
               examples/intermediate/ExtractClose.cpp \
               examples/intermediate/ComputeTriDual.cpp \
               examples/intermediate/FixPolys.cpp
examples_advanceddir = $(datadir)/examples/advanced
examples_advanced_DATA = \
               examples/advanced/makefile \
               examples/advanced/CrystalRouterExample.cpp \
               examples/advanced/ExtrudePoly.cpp \
               examples/advanced/GenLargeMesh.cpp \
               examples/advanced/LloydRelaxation.cpp \
               examples/advanced/CoupleMGen.cpp \
               examples/advanced/VisTags.cpp
examples_advancedsmoothdir = $(datadir)/examples/advanced/smoothing
examples_advancedsmooth_DATA = \
               examples/advanced/smoothing/makefile \
               examples/advanced/smoothing/LaplacianSmoother.cpp \
               examples/advanced/smoothing/OptimizeMeshMesquite.cpp
examples_fortrandir = $(datadir)/examples/fortran
examples_fortran_DATA = \
               examples/fortran/makefile \
               examples/fortran/DirectAccessNoHolesF90.F90 \
               examples/fortran/PushParMeshIntoMoabF90.F90 \
               examples/fortran/MigrateMesh.F90 \
               examples/fortran/spherecube.h5m
examples_pythonexdir = $(datadir)/examples/python
examples_pythonex_DATA = \
               examples/python/makefile \
               examples/python/laplaciansmoother.py \
               examples/python/H5MToScrip.py examples/python/ScripToH5M.py \
               examples/python/gr3rd.py \
               examples/python/interrogate_mesh.py \
               examples/python/scrip3d.py 
examples_appsdir = $(datadir)/examples/apps
examples_apps_DATA = \
               examples/apps/makefile
examples_apps_earthdir = $(datadir)/examples/apps/earth
examples_apps_earth_DATA = \
               examples/apps/earth/makefile \
               examples/apps/earth/ContinentsOnGlobe.cpp \
               examples/apps/earth/BoundaryContinents.cpp \
               examples/apps/earth/addPCdata.cpp \
               examples/apps/earth/AddFieldtoPG2.cpp \
               examples/apps/earth/copyPartition.cpp

CMAKE_CONFIGURATION_FILES = CMakeLists.txt \
   config/AutoconfHeader.cmake      config/FindCGM.cmake          config/FindNetCDF.cmake        config/FindZoltan.cmake     config/GetAcInitVersion.cmake \
   config/CMakeDistCheck.cmake      config/FindEigen3.cmake       config/FindPNetCDF.cmake       config/Findgraphviz.cmake   config/MOABConfig.cmake.cmake \
   config/CMakeLibraryFuncs.cmake   config/FindHDF5_MOAB.cmake    config/FindParMetis.cmake      config/Findverdict.cmake    config/ResolveCompilerPaths.cmake \
   config/CheckCompilerFlags.cmake  config/FindMetis.cmake        config/FindTempestRemap.cmake  config/ForceAddFlags.cmake  config/TestBigEndian.cmake \
   config/CMakeReplicateConfig.sh   config/TemplateSpecialization.cxx   config/MOABConfig.cmake.h.in   config/cmake_uninstall.cmake.in moab.make.cmake  \
   pymoab/CMakeLists.txt            pymoab/cmake/FindPythonDev.cmake  pymoab/setup.py.cmake.in     pymoab/install.sh.cmake

EXTRA_DIST = $(doc_DATA) \
						 debian \
						 config/MOAB.pc.in config/MOAB.pc.cmake.in \
						 $(CMAKE_CONFIGURATION_FILES) \
						 $(examples_DATA) \
						 $(examples_cmake_DATA) \
						 $(examples_beginner_DATA) \
						 $(examples_basic_DATA) \
						 $(examples_intermediate_DATA) \
						 $(examples_advanced_DATA) \
						 $(examples_advancedsmooth_DATA) \
						 $(examples_fortran_DATA) \
						 $(examples_pythonex_DATA) \
						 $(examples_apps_DATA) \
						 $(examples_apps_earth_DATA) \
						 suggest_configuration.sh

# Automake doesn't seem to have a directory defined for
# platform-dependent data (or include) files. So put
# in $(libdir).  Define a $(cfgdir) to get around automake's
# check that only libraries are going in $(libdir)
cfgdir = $(libdir)
cfg_DATA = moab.make moab.config
cmakedir = $(libdir)/cmake/MOAB
cmake_DATA = MOABConfig.cmake config/ResolveCompilerPaths.cmake

# By default, moab.make will define these to $(srcdir).  We
# want to override that during the INSTALL of the file so
# that the correct values are set (e.g. if someone does
# 'make prefix=/foo install', we don't know the correct install
# directory until we're doing the install.
install-data-hook:
	$(AM_V_at)echo "MOAB_LIBDIR=${libdir}" >> $(DESTDIR)$(cfgdir)/moab.make
	$(AM_V_at)echo "MOAB_INCLUDES=-I${includedir}" >> $(DESTDIR)$(cfgdir)/moab.make

# Generate a file to be installed in $libdir containing the configuration
# options used for this MOAB build.
# Note: If you modify this, verify that it works correctly for option
#       values containing spaces.
# First line: get unformatted data from config.status
# Second line: extract just the list of options
# Third line: put individual options on separate lines
# Fourth line: remove any empty lines and trim spaces and then write to config file.
moab.config: config.status Makefile.am
	$(AM_V_GEN)./config.status --version | \
	  sed -e 's/.*options "\(.*\)"/\1/p' -e 'd' | \
	  tr "'" "\n" | \
	  sed -e '/^ *$$/d' > $@

CLEANFILES = moab.config configs.sed
DISTCLEANFILES = a.out
