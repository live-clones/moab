include $(top_srcdir)/src/mesquite/make.inc

MESHDIR = $(abs_top_srcdir)/MeshFiles/mesquite

AM_CPPFLAGS  = -DSRCDIR=$(srcdir) \
               -DMESHDIR=$(MESHDIR) \
               $(MSQ_INCLUDES) \
               -I$(top_srcdir)/src \
               -I$(top_srcdir)/src/parallel \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/test

LDADD = $(top_builddir)/src/mesquite/libmbmesquite.la

check_PROGRAMS = deforming \
								 laplacian \
                 shape_improvement \
                 shape_improver \
                 size_adapt_shape \
                 untangle \
                 viscous_cfd

TESTS = $(check_PROGRAMS)

deforming_SOURCES = deforming.cpp
laplacian_SOURCES = laplacian.cpp
shape_improvement_SOURCES = shape_improvement.cpp
shape_improver_SOURCES = shape_improver.cpp
size_adapt_shape_SOURCES = size_adapt_shape.cpp
untangle_SOURCES = untangle.cpp
viscous_cfd_SOURCES = viscous_cfd.cpp

EXTRA_DIST = CMakeLists.txt \
             sph-10-zsquare.vtk

MOSTLYCLEANFILES = original_mesh.vtk smoothed_mesh.vtk
DISTCLEANFILES = a.out
