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
check_PROGRAMS = par_hex_smooth_laplace
par_hex_smooth_laplace_SOURCES = par_hex_smooth_laplace.cpp
par_hex_smooth_laplace_DEPENDENCIES = $(LDADD)
PARALLEL_TESTS = $(check_PROGRAMS)
EXTRA_DIST = CMakeLists.txt
TESTS = $(PARALLEL_TESTS)

# Run parallel tests with mpiexec
LOG_COMPILER = `if echo "$(PARALLEL_TESTS)" | $(GREP) "$$p" > /dev/null; then echo "${MPIEXEC} ${MPIEXEC_NP} ${NP}"; else echo ""; fi`

MOSTLYCLEANFILES = *.*.vtk
DISTCLEANFILES = a.out
