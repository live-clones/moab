include $(top_srcdir)/src/mesquite/make.inc

LDADD = libmsqutil.la $(top_builddir)/src/mesquite/libmbmesquite.la
bin_PROGRAMS = vtkrandom vtkxform msqshape msqquality
noinst_LTLIBRARIES = libmsqutil.la

AM_CPPFLAGS  = $(MSQ_INCLUDES) \
               -I$(top_srcdir)/itaps \
               -I$(top_srcdir)/src \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/itaps/imesh \
               -I$(top_builddir)/itaps/imesh \
               -I$(top_srcdir)/itaps/fbigeom \
               -I$(top_srcdir)/itaps/irel

libmsqutil_la_SOURCES = CLArgs.cpp CLArgs.hpp \
                CLArgFlag.cpp CLArgFlag.hpp \
                ManPage.cpp ManPage.hpp

vtkrandom_SOURCES = randomize.cpp domain.hpp domain.cpp
msqshape_SOURCES = shape-improve.cpp domain.hpp domain.cpp
vtkxform_SOURCES = transform.cpp
msqquality_SOURCES = quality.cpp domain.hpp domain.cpp

MOSTLYCLEANFILES = $(dist_man1_MANS)
EXTRA_DIST = CMakeLists.txt
DISTCLEANFILES = a.out

