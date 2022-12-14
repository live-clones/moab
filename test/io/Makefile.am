# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign no-dependencies

LDADD = $(top_builddir)/src/libMOAB.la

MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS  = $(BGDEF)-DSRCDIR=$(srcdir) \
               $(BGDEF)-DMESHDIR=$(MESHDIR) \
               -I$(builddir) -I.. -I$(srcdir)/.. \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/io \
               -I$(top_srcdir)/src/parallel

if HAVE_NETCDF
  EXODUS_TEST = exodus_test
  NC_TEST = read_nc read_ucd_nc read_mpas_nc write_nc read_gcrm_nc
  LDADD += $(NETCDF_LIBS) $(PNETCDF_LIBS) $(HDF5_LIBS)
else
  EXODUS_TEST =
  NC_TEST =
endif

if HAVE_CCMIO
  CCMIO_TEST = ccmio_test
  LDADD += $(CCMIO_LIBS)
else
  CCMIO_TEST =
endif

if HAVE_CGM
if HAVE_OCC
  CGM_TEST = read_cgm_load_test read_cgm_basic_test read_cgm_senses_test read_cgm_connectivity_test \
             read_cgm_group_test read_cgm_facet_test
endif
else
  CGM_TEST =
endif

if HAVE_CGNS
  CGNS_TEST = cgns_test
  LDADD += $(CGNS_LIBS)
else
  CGNS_TEST =
endif

TESTS = readutil_test \
	tqdcfr \
        vtk_test \
        cub_file_test \
	$(EXODUS_TEST) \
        gmsh_test \
	nastran_test \
	$(NC_TEST) \
        ideas_test \
        stl_test \
        smf_test \
	$(CCMIO_TEST) \
	$(CGM_TEST) \
	$(CGNS_TEST) \
	readrtt_test \
	read_obj_test

check_PROGRAMS = $(TESTS)

tqdcfr_SOURCES = tqdcfr_test.cpp

#test_rms_SOURCES = test_rms.cpp
#test_exo_SOURCES = test_exo.cpp
vtk_test_SOURCES = $(srcdir)/../TestUtil.hpp VtkTest.cpp
cub_file_test_SOURCES = $(srcdir)/../TestUtil.hpp cub_file_test.cpp
exodus_test_SOURCES = $(srcdir)/../TestUtil.hpp exodus_test.cpp
gmsh_test_SOURCES = $(srcdir)/../TestUtil.hpp gmsh_test.cpp
nastran_test_SOURCES = $(srcdir)/../TestUtil.hpp nastran_test.cpp
read_nc_SOURCES = $(srcdir)/../TestUtil.hpp read_nc.cpp
read_ucd_nc_SOURCES = $(srcdir)/../TestUtil.hpp read_ucd_nc.cpp
read_mpas_nc_SOURCES = $(srcdir)/../TestUtil.hpp read_mpas_nc.cpp
read_gcrm_nc_SOURCES = $(srcdir)/../TestUtil.hpp read_gcrm_nc.cpp
write_nc_SOURCES = $(srcdir)/../TestUtil.hpp write_nc.cpp
ideas_test_SOURCES = $(srcdir)/../TestUtil.hpp ideas_test.cpp
stl_test_SOURCES = $(srcdir)/../TestUtil.hpp stl_test.cpp
smf_test_SOURCES = $(srcdir)/../TestUtil.hpp smf_test.cpp
ccmio_test_SOURCES = $(srcdir)/../TestUtil.hpp ccmio_test.cpp
read_cgm_load_test_SOURCES = $(srcdir)/../TestUtil.hpp read_cgm_load_test.cpp
read_cgm_load_test_CXXFLAGS = $(CGM_CPPFLAGS)  $(CXXFLAGS)
read_cgm_basic_test_SOURCES = $(srcdir)/../TestUtil.hpp read_cgm_basic_test.cpp
read_cgm_basic_test_CXXFLAGS = $(CGM_CPPFLAGS) $(CGM_LIBS)  $(CXXFLAGS)
read_cgm_senses_test_SOURCES = $(srcdir)/../TestUtil.hpp read_cgm_senses_test.cpp
read_cgm_senses_test_CXXFLAGS = $(CGM_CPPFLAGS) $(CGM_LIBS)  $(CXXFLAGS)
read_cgm_connectivity_test_SOURCES = $(srcdir)/../TestUtil.hpp read_cgm_connectivity_test.cpp
read_cgm_connectivity_test_CXXFLAGS = $(CGM_CPPFLAGS) $(CGM_LIBS)  $(CXXFLAGS)
read_cgm_group_test_SOURCES = $(srcdir)/../TestUtil.hpp read_cgm_group_test.cpp
read_cgm_group_test_CXXFLAGS = $(CGM_CPPFLAGS) $(CGM_LIBS)  $(CXXFLAGS)
read_cgm_facet_test_SOURCES = $(srcdir)/../TestUtil.hpp read_cgm_facet_test.cpp
read_cgm_facet_test_CXXFLAGS = $(CGM_CPPFLAGS) $(CGM_LIBS) $(CXXFLAGS)
readutil_test_SOURCES = $(srcdir)/../TestUtil.hpp readutil_test.cpp
cgns_test_SOURCES=$(srcdir)/../TestUtil.hpp cgns_test.cpp
readrtt_test_SOURCES = $(srcdir)/../TestUtil.hpp readrtt_test.cpp
read_obj_test_SOURCES = $(srcdir)/../TestUtil.hpp read_obj_test.cpp

MOSTLYCLEANFILES = dumped_acis.sat \
                   gcrm.h5m \
                   mb_write.g \
                   mb_write2.g \
                   polyg.exo \
                   polyh.exo \
                   singlecyl_tmp.ccmg \
                   test_eul_across_files.nc \
                   test_eul_append.nc \
                   test_eul_ghosting.nc \
                   test_eul_T2.nc \
                   test_eul_T.nc \
                   test_fv_T.nc \
                   test_gcrm_vars.nc \
                   test_homme_T.nc \
                   test_mpas_vars.nc \
                   rtttest.rtt
DISTCLEANFILES = a.out
