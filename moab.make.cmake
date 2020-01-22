# The values below are for a copy of MOAB used directly
# from its build directory. These values will be overridden below
# for installed copies of MOAB.

# Feature list
MOAB_MPI_ENABLED = @MOAB_HAVE_MPI@
MOAB_FORTRAN_ENABLED = @ENABLE_FORTRAN@
MOAB_HDF5_ENABLED = @MOAB_HAVE_HDF5@
MOAB_NETCDF_ENABLED = @MOAB_HAVE_NETCDF@
MOAB_PNETCDF_ENABLED = @MOAB_HAVE_PNETCDF@
MOAB_IGEOM_ENABLED = @MOAB_HAVE_CGM@
MOAB_IMESH_ENABLED = @ENABLE_IMESH@
MOAB_IREL_ENABLED = @ENABLE_IREL@
MOAB_FBIGEOM_ENABLED = @ENABLE_FBIGEOM@
MOAB_MESQUITE_ENABLED = @ENABLE_MESQUITE@

# Library and Include paths
MOAB_LIBDIR = @abs_builddir@/src/.libs
MOAB_INCLUDES = -I@abs_srcdir@/src \
                -I@abs_builddir@/src \
                -I@abs_srcdir@/src/oldinc \
                -I@abs_srcdir@/src/verdict \
                -I@abs_srcdir@/src/parallel \
                -I@abs_builddir@/src/parallel \
                -I@abs_srcdir@/src/LocalDiscretization \
                -I@abs_srcdir@/src/RefineMesh

ifeq ($(MOAB_IGEOM_ENABLED),yes)
include @CGM_MAKE@
MOAB_INCLUDES += $(IGEOM_INCLUDES)
endif

ifeq ($(MOAB_MESQUITE_ENABLED),yes)
MSQ_LIBS = -L@abs_builddir@/src/mesquite/.libs -lmbmesquite
MSQ_INCLUDES = -I@abs_srcdir@/src/mesquite/include \
	-I@abs_builddir@/src/mesquite/include \
	-I@abs_srcdir@/src/mesquite/Mesh \
	-I@abs_srcdir@/src/mesquite/Control \
	-I@abs_srcdir@/src/mesquite/Wrappers \
	-I@abs_srcdir@/src/mesquite/MappingFunction \
	-I@abs_srcdir@/src/mesquite/MappingFunction/Lagrange \
	-I@abs_srcdir@/src/mesquite/MappingFunction/Linear \
	-I@abs_srcdir@/src/mesquite/Misc \
	-I@abs_srcdir@/src/mesquite/ObjectiveFunction \
	-I@abs_srcdir@/src/mesquite/QualityAssessor \
	-I@abs_srcdir@/src/mesquite/QualityImprover \
	-I@abs_srcdir@/src/mesquite/QualityImprover/OptSolvers \
	-I@abs_srcdir@/src/mesquite/QualityImprover/Relaxation \
	-I@abs_srcdir@/src/mesquite/QualityMetric \
	-I@abs_srcdir@/src/mesquite/QualityMetric/Debug \
	-I@abs_srcdir@/src/mesquite/QualityMetric/Shape \
	-I@abs_srcdir@/src/mesquite/QualityMetric/Smoothness \
	-I@abs_srcdir@/src/mesquite/QualityMetric/TMP \
	-I@abs_srcdir@/src/mesquite/QualityMetric/Untangle \
	-I@abs_srcdir@/src/mesquite/QualityMetric/Volume \
	-I@abs_srcdir@/src/mesquite/TargetCalculator \
	-I@abs_srcdir@/src/mesquite/TargetMetric \
	-I@abs_srcdir@/src/mesquite/TargetMetric/Misc \
	-I@abs_srcdir@/src/mesquite/TargetMetric/Shape \
	-I@abs_srcdir@/src/mesquite/TargetMetric/ShapeOrient \
	-I@abs_srcdir@/src/mesquite/TargetMetric/ShapeSize \
	-I@abs_srcdir@/src/mesquite/TargetMetric/ShapeSizeOrient \
	-I@abs_srcdir@/src/mesquite/TargetMetric/Size \
	-I@abs_srcdir@/src/mesquite/TargetMetric/Untangle

MOAB_INCLUDES += ${MSQ_INCLUDES}
endif

MOAB_INCLUDES += @ZOLTAN_INC_FLAGS@

MOAB_CPPFLAGS = @CPPFLAGS@ 
MOAB_CXXFLAGS = @CXXFLAGS@ 
MOAB_CFLAGS = @CFLAGS@ 
MOAB_FFLAGS = @FFLAGS@
MOAB_FCFLAGS = @FCFLAGS@
MOAB_LDFLAGS = @EXPORT_LDFLAGS@ @CXX_LDFLAGS@ @LDFLAGS@

# missing support for DAMSEL, CCMIO
# PROBABLY FORMATTED INCORRECTLY
MOAB_EXT_LIBS = @NETCDF_LIBRARIES@ @PNETCDF_LIBRARIES@ @CGNS_LIBRARIESS@ @HDF5_LIBRARIES@ \
                @CGM_LIBRARIES@ @TEMPESTREMAP_LIBRARIES@ @ZOLTAN_LIBRARIES@ @PARMETIS_LIBRARIES@ \
                @METIS_LIBRARIES@ @HYPRE_LIBRARIES@ @LAPACK_LIBRARIES@ @BLAS_LIBRARIES@ @LIBS@ 
MOAB_LIBS_LINK = ${MOAB_LDFLAGS} -L${MOAB_LIBDIR} -lMOAB ${MSQ_LIBS} $(MOAB_EXT_LIBS) 
DAGMC_LIBS_LINK = ${MOAB_LDFLAGS} -L${MOAB_LIBDIR} @DAGMC_LIBS@ -lMOAB ${MSQ_LIBS} $(MOAB_EXT_LIBS)

MOAB_CXX = @CXX@
MOAB_CC  = @CC@
MOAB_FC  = @FC@
MOAB_F77  = @F77@

# Override MOAB_LIBDIR and MOAB_INCLUDES from above with the correct
# values for the installed MOAB.

# NEED TO ADD SOMETHING TO MODIFY THIS FILE AT INSTALL TIME (OR FIX ISSUE #30)
