# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = include/Mesquite_all_headers.hpp
lib_LTLIBRARIES = libmbmesquite.la
libmbmesquite_la_DEPENDENCIES=
# Tell libtool that there are no undefined symbols when linking
# the library.  This means that either the symbol is defined in
# the library or we include a -lname argument to specify which
# other library provides the symbol.  Libtool requires this
# guarantee to generate libraries on platforms for which shared
# libs cannot have any unmet depenendencies (windows, and maybe AIX)
libmbmesquite_la_LDFLAGS = $(AM_LDFLAGS) -no-undefined

include make.inc

libmbmesquite_la_includedir = $(includedir)/moab/mesquite
libmbmesquite_la_nodist_includedir = $(includedir)/moab

libmbmesquite_la_SOURCES = $(MSQ_SRCS)

libmbmesquite_la_include_HEADERS = \
																	 include/Mesquite.hpp \
																	 include/MeshInterface.hpp \
																	 include/ParallelMeshInterface.hpp \
																	 include/ParallelHelperInterface.hpp \
																	 include/CurveDomain.hpp \
																	 $(MSQ_HDRS)

# don't put these generated headers into the distribution
nodist_libmbmesquite_la_include_HEADERS = include/Mesquite_all_headers.hpp

libmbmesquite_la_LIBADD = $(MESQUITE_LIBS)
AM_CPPFLAGS  = $(MSQ_INCLUDES) \
               -I$(top_srcdir)/src \
               -I$(top_builddir)/src \
               -I$(top_srcdir)/itaps \
               -I$(top_srcdir)/itaps/imesh \
               -I$(top_builddir)/itaps/imesh \
               -I$(top_srcdir)/itaps/fbigeom \
               -I$(top_srcdir)/itaps/irel

include/Mesquite_all_headers.hpp: $(MSQ_HDRS)
	@echo "Generating $@..."
	@mkdir -p include
	@echo '/* NOTE: This file is generated by make.' > $@
	@echo ' *       Do not manually edit it.'       >> $@
	@echo ' */'                                     >> $@
	@echo '#ifndef MESQUITE_ALL_HEADERS_HPP'        >> $@
	@echo '#define MESQUITE_ALL_HEADERS_HPP'        >> $@
	@echo '#include "Mesquite.hpp"'                 >> $@
	@echo '#include "MeshInterface.hpp"'            >> $@
	@echo '#include "ParallelMeshInterface.hpp"'    >> $@
	@echo '#include "ParallelHelperInterface.hpp"'  >> $@
	@echo '#include "CurveDomain.hpp"'  >> $@
	@for file in $(MSQ_HDRS); do \
	   n=`expr X"$$file" : X".*/\([^/]*\)"`;\
	   test -z "$$n" || echo "#include \"$$n\"" >> $@; \
	done
	@echo '#endif'                                  >> $@

MOSTLYCLEANFILES = $(BUILT_SOURCES)
DISTCLEANFILES = $(OLD_AM_GENERATED_HDRS) a.out

EXTRA_DIST = CMakeLists.txt

