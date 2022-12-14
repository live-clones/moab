# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign subdir-objects no-dependencies

SUBDIRS = . io dual obb perf
if HAVE_HDF5
  SUBDIRS += h5file
endif
if PARALLEL
  SUBDIRS += parallel
endif
if ENABLE_mesquite
  SUBDIRS += mesquite
endif
LDADD = $(top_builddir)/src/libMOAB.la

MESHDIR = $(abs_top_srcdir)/MeshFiles/unittest
RMESHDIR = $(top_srcdir)/MeshFiles/unittest

AM_CPPFLAGS  = -DIS_BUILDING_MB \
               -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/src  \
               -I$(top_srcdir)/src/LocalDiscretization \
               -I$(top_srcdir)/src/verdict \
               -I$(top_srcdir)/src/DiscreteGeometry \
               -I$(top_srcdir)/src/IntxMesh


if PARALLEL
  AM_CPPFLAGS += -I$(top_srcdir)/src/parallel
endif

TESTS = range_test \
        moab_test \
        homxform_test \
        scdseq_test \
        scd_test_partn \
        test_adj \
        test_determinism \
        seq_man_test \
        geom_util_test \
        xform_test \
        oriented_box_test \
        adaptive_kd_tree_tests \
        bsp_tree_test \
        file_options_test \
        kd_tree_test \
        var_len_test var_len_test_no_template \
        tag_test \
        mesh_set_test \
        mbcn_test \
        bsp_tree_poly_test \
        reorder_test \
        test_prog_opt \
        coords_connect_iterate \
        elem_eval_test \
        spatial_locator_test \
        test_boundbox \
        adj_moab_test \
        uref_mesh_test \
        verdict_test \
        test_Matrix3 \
        hireconst_test \
        uref_hirec_test \
        spherical_area_test \
        arc_intx_tests \
        intx_in_plane_test \
        intx_on_sphere_test \
        mgen_test \
        gnom_project_test

if HAVE_NETCDF
   TESTS += intx_rll_cs_sphere_test
endif

if HAVE_HDF5
  TESTS += mbfacet_test \
        gttool_test \
		restore_topo_geom_incl_test \
        cropvol_test \
        mergemesh_test \
        mbground_test \
        lloyd_smoother_test \
        mhdf_public_test \
	test_geom_gqt \
	gqt_simple_test \
	gqt_rayfire_test \
	gqt_pointinvol_test \
	gqt_findvol_test

# imoab tests
if HAVE_32BIT

else
TESTS += imoab_test
if ENABLE_FORTRAN
  TESTS += imoab_testF
endif
endif

endif

check_PROGRAMS = $(TESTS) \
		 kd_tree_tool \
		 kd_tree_time

MOSTLYCLEANFILES = mb_write.g \
                   mb_write2.g \
                   hex_mid_volume_nodes.g \
                   hex_mid_edge_face_vol_nodes.g \
                   hex_mid_edge_nodes.g \
                   tet_mid_edge_nodes.g \
                   tet_mid_face_nodes.g \
                   tet_mid_edge_face_nodes.g \
                   tri_mid_edge_face_nodes.g \
                   mb_stress_out.g \
                   merge_test.g \
                   mm_out.h5m \
                   tree.h5m \
                   TaskMesh_n1.0.h5m \
                   TaskMesh_n1.1.h5m \
                   fIntxOverlap.h5m fIntxTarget.h5m fnew.h5m fnew2.h5m out1.h5m \
                   intx.vtk \
                   polyWithEdges.vtk \
                   intx1.vtk intx.h5m project.vtk

# Tests and such
test_Matrix3_SOURCES = test_Matrix3.cpp
mgen_test_SOURCES = mgen_test.cpp
if ENABLE_FORTRAN
imoab_testF_SOURCES = imoab_testF.F90
imoab_testF_FCLAGS = -DMESHDIR2=\"$(RMESHDIR)\" $(FCLAGS)
#imoab_testF_LDADD = $(LDADD) $(FCLIBS)
if HAVE_TEMPESTREMAP
imoab_testF_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS) $(FCLIBS)
else
imoab_testF_LDADD = $(LDADD) $(FCLIBS)
endif
endif
imoab_test_SOURCES = imoab_test.c
if HAVE_TEMPESTREMAP
imoab_test_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
else
imoab_test_LDADD = $(LDADD)
endif
CCLD = $(CXX)

verdict_test_SOURCES = verdict_test.cpp
range_test_SOURCES = range_test.cpp
moab_test_SOURCES = MBTest.cpp
moab_test_DEPENDENCIES = \
	mb_big_test.g \
	cell1.gen \
	cell2.gen \
	$(LDADD)

homxform_test_SOURCES = homxform_test.cpp
scdseq_test_SOURCES = scdseq_test.cpp
scd_test_partn_SOURCES = scd_test_partn.cpp
#scdseq_timing_SOURCES = scdseq_timing.cpp
#merge_test_SOURCES = merge_test.cpp
#test_exo_SOURCES = test_exo.cpp
#test_tag_server_SOURCES = $(top_srcdir)/src/TagServer.cpp
#test_tag_server_CPPFLAGS = -DTEST $(CPPFLAGS)
test_adj_SOURCES = test_adj.cpp
test_determinism_SOURCES = test_determinism.cpp
seq_man_test_SOURCES = TestUtil.hpp TestTypeSequenceManager.cpp
mesh_set_test_SOURCES = Test_MBMeshSet.cpp TestUtil.hpp
geom_util_test_SOURCES = TestUtil.hpp GeomUtilTests.cpp
xform_test_SOURCES = affinexform_test.cpp
oriented_box_test_SOURCES = oriented_box_test.cpp
adaptive_kd_tree_tests_SOURCES = adaptive_kd_tree_tests.cpp
kd_tree_tool_SOURCES = kd_tree_tool.cpp
kd_tree_time_SOURCES = kd_tree_time.cpp
kd_tree_test_SOURCES = kd_tree_test.cpp
bsp_tree_test_SOURCES = bsp_tree_test.cpp
reorder_test_SOURCES = reorder_test.cpp
elem_eval_test_SOURCES = elem_eval_test.cpp
file_options_test_SOURCES = TestUtil.hpp fileopts_test.cpp
var_len_test_SOURCES = TestUtil.hpp VarLenTagTest.cpp
var_len_test_no_template_SOURCES = $(var_len_test_SOURCES)
var_len_test_no_template_CPPFLAGS = -UTEMPLATE_SPECIALIZATION $(AM_CPPFLAGS) $(CPPFLAGS)
tag_test_SOURCES = TestUtil.hpp TagTest.cpp
spatial_locator_test_SOURCES = TestUtil.hpp spatial_locator_test.cpp
mergemesh_test_SOURCES = mergemesh_test.cpp
# intersection related tests
intx_on_sphere_test_SOURCES = intx_on_sphere_test.cpp
intx_in_plane_test_SOURCES = intx_in_plane_test.cpp
spherical_area_test_SOURCES = spherical_area_test.cpp
arc_intx_tests_SOURCES = arc_intx_tests.cpp
intx_rll_cs_sphere_test_SOURCES = intx_rll_cs_sphere_test.cpp
gnom_project_test_SOURCES = gnom_project_test.cpp

#mbcn_test_SOURCES = $(top_srcdir)/src/moab/CN.hpp \
#                    $(top_srcdir)/src/CN.cpp \
#                    mbcn_test.cpp
mbcn_test_SOURCES = mbcn_test.cpp

bsp_tree_poly_test_SOURCES = bsp_tree_poly_test.cpp
mbfacet_test_SOURCES = mbfacet_test.cpp
mbground_test_SOURCES = mbground_test.cpp
gttool_test_SOURCES = gttool_test.cpp
restore_topo_geom_incl_test_SOURCES = restore_topo_geom_incl_test.cpp
cropvol_test_SOURCES = crop_vol_test.cpp
mhdf_public_test_SOURCES = mhdf_public_test.cpp
mhdf_public_test_LDADD = $(LDADD) $(HDF5_LIBS)

test_prog_opt_SOURCES = test_prog_opt.cpp
test_prog_opt_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS)
moab_test_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS)
kd_tree_test_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS)

coords_connect_iterate_SOURCES = coords_connect_iterate.cpp
coords_connect_iterate_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS)

test_boundbox_SOURCES = test_boundbox.cpp
lloyd_smoother_test_SOURCES = lloyd_smoother_test.cpp

adj_moab_test_SOURCES = adj_moab_test.cpp
uref_mesh_test_SOURCES = urefine_mesh_test.cpp
hireconst_test_SOURCES = hireconst_test.cpp
uref_hirec_test_SOURCES = uref_hirec_test.cpp geomObject.cpp

test_geom_gqt_SOURCES = test_geom_gqt.cpp

gqt_simple_test_SOURCES = gqt_simple_test.cpp
gqt_rayfire_test_SOURCES = gqt_rayfire_test.cpp
gqt_pointinvol_test_SOURCES = gqt_pointinvol_test.cpp
gqt_findvol_test_SOURCES = gqt_findvol_test.cpp

if PARALLEL
moab_test_CPPFLAGS += -I$(top_srcdir)/src/parallel
kd_tree_test_CPPFLAGS += -I$(top_srcdir)/src/parallel
test_prog_opt_CPPFLAGS += -I$(top_srcdir)/src/parallel
endif

# Add thread flags/libraries for the TypeSequenceManagerTest
seq_man_test_CPPFLAGS = $(AM_CPPFLAGS) $(PTHREAD_CFLAGS)
seq_man_test_LDFLAGS = $(AM_LDFLAGS) $(PTHREAD_LIBS)

# Tests for remapping
if ENABLE_mbtempest
if HAVE_TEMPESTREMAP
  test_remapping_SOURCES = test_remapping.cpp
  test_remapping_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/src/Remapping -I$(top_srcdir)/src/IntxMesh $(TEMPESTREMAP_INCLUDE)
  test_remapping_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  imoab_remapping_SOURCES = imoab_remap.cpp
  imoab_remapping_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/src/Remapping -I$(top_srcdir)/src/IntxMesh $(TEMPESTREMAP_INCLUDE)
  imoab_remapping_LDADD = $(LDADD) $(TEMPESTREMAP_LIBS)
  TESTS += test_remapping imoab_remapping
endif
endif

FCDEFS += $(FC_DEFINE)SRCDIR=$(srcdir) \
          $(FC_DEFINE)MESHDIR=$(MESHDIR)

FORTDEFS = $(FCDEFS)
PPFCCOMPILE = $(FC) $(FORTDEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)
LTPPFCCOMPILE = $(LIBTOOL) $(AM_V_lt) --tag=FC $(AM_LIBTOOLFLAGS) \
    $(LIBTOOLFLAGS) --mode=compile $(FC) $(FCDEFS) \
    $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
    $(AM_FCFLAGS) $(FCFLAGS)
PPF77COMPILE = $(F77) $(FORTDEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
   $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FFLAGS) $(FFLAGS)
LTPPF77COMPILE = $(LIBTOOL) $(AM_V_lt) --tag=F77 $(AM_LIBTOOLFLAGS) \
   $(LIBTOOLFLAGS) --mode=compile $(F77) $(FCDEFS) \
   $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
   $(AM_FFLAGS) $(FFLAGS)


CMAKELISTS_FILES = CMakeLists.txt dual/CMakeLists.txt h5file/CMakeLists.txt io/CMakeLists.txt obb/CMakeLists.txt perf/CMakeLists.txt perf/point_location/CMakeLists.txt

EXTRA_DIST = TestUtil.hpp TestRunner.hpp $(CMAKELISTS_FILES)

mb_big_test.g: $(MESHDIR)/mb_big_test.g.gz
	$(AM_V_GEN)$(ZCAT) $< > $@

cell1.gen: $(MESHDIR)/cell1.gen.gz
	$(AM_V_GEN)$(ZCAT) $< > $@

cell2.gen: $(MESHDIR)/cell2.gen.gz
	$(AM_V_GEN)$(ZCAT) $< > $@


CLEANFILES = mb_big_test.g \
             cell1.gen \
             cell2.gen \
             outTempestCS.g \
             outTempestRLL.g \
             outTempestICO.g \
             outTempestICOD.g \
             outTempestOV.g \
             geom_missing_obb.h5m \
             geom_w_obbs.h5m \
             fIntxLndTarget.h5m \
             fIntxAtmTarget.h5m \
             fIntxOcnTarget.h5m \
             atm_ocn_map.nc fOcnOnCpl.h5m fOcnOnCpl2.h5m
DISTCLEANFILES = a.out

# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check
