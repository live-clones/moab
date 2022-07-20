dnl Process this file with autoconf 2.60 to produce a configure script
AC_DEFUN([PRINT_CONFIGURATION_SUMMARY],[

  # Some aliases for colors to pretty-print
  NORMAL=$(tput sgr0)
  GREEN=$(tput setaf 2)
  RED=$(tput setaf 1)
  BLUE=$(tput setaf 4)
  BOLD=$(tput bold)
  UNDERLINE_BEG=$(tput smul)
  UNDERLINE_END=$(tput rmul)
  BLINK=$(tput rev)

  MSG_ECHO_SEPARATOR

  printf '%-40s' ""
  MSG_ECHO_CUSTOM([${UNDERLINE_BEG}${BOLD}MOAB Configuration Summary${UNDERLINE_END}${NORMAL} ])
  BLANK_LINE

  SUMMARY_LINE([Installation Prefix], [$prefix])
  SUMMARY_LINE([Debug Mode], [$enable_debug])
  SUMMARY_LINE([Optimized Mode], [$enable_cxx_optimize])
  SUMMARY_LINE([Static Build], [$enable_static])
  SUMMARY_LINE([Shared Build], [$enable_shared])
  SUMMARY_LINE([BLAS/LAPACK support], [$enable_blaslapack])
  SUMMARY_LINE([Python support], [$enable_pymoab])
  SUMMARY_LINE([Eigen3 support], [$enableeigen])
  case "x$CGM_MISSING" in
   "xno") SUMMARY_LINE([CGM support], [yes]) ;;
   *) SUMMARY_LINE([CGM support], [no]) ;;
  esac
  SUMMARY_LINE([MPI parallelism support], [$enablempi])
  
  SUMMARY_LINE([SZip support], [$enableszip])
  SUMMARY_LINE([ZLib support], [$enablezlib])
  case "x$enablehdf5" in
   "xyes") SUMMARY_LINE([HDF5 support], [yes], [(parallel = $enablehdf5parallel)]) ;;
   *) SUMMARY_LINE([HDF5 support], [no], [  ${RED}(No support for native file format)${NORMAL}]) ;;
  esac
  case "x$enablenetcdf" in
   "xyes") SUMMARY_LINE([NetCDF support], [yes]) ;;
   *) SUMMARY_LINE([NetCDF support], [no], [  ${RED}(Support for ExodusII disabled)${NORMAL}]) ;;
  esac
  SUMMARY_LINE([PNetCDF support], [$enablepnetcdf])
  SUMMARY_LINE([Metis support], [$enablemetis])
  SUMMARY_LINE([ParMetis support], [$enableparmetis])
  SUMMARY_LINE([Zoltan support], [$enablezoltan])
  SUMMARY_LINE([TempestRemap support], [$enabletempestremap])

  BLANK_LINE
  MSG_ECHO_SEPARATOR
  BLANK_LINE

  printf '%-42s' ""
  MSG_ECHO_CUSTOM([${UNDERLINE_BEG}${BOLD}Compiler Flag Summary${UNDERLINE_END}${NORMAL} ])
  BLANK_LINE

  SUMMARY_LINE([C compiler], [$CC])
  SUMMARY_LINE([C++ compiler], [$CXX])
  SUMMARY_LINE([Fortran compiler], [$FC])
  BLANK_LINE

  SUMMARY_LINE([C],   [$CC $CFLAGS $CPPFLAGS])
  SUMMARY_LINE([C++], [$CXX $CXXFLAGS $CPPFLAGS])
  if (test "xno" != "x$ENABLE_FORTRAN"); then
    SUMMARY_LINE([Fortran90], [$FC $FCFLAGS $FCPPFLAGS])
  fi
  if (test "xno" != "x$ENABLE_FORTRAN"); then
    SUMMARY_LINE([Fortran77], [$F77 $FFLAGS $FCPPFLAGS])
  fi
  SUMMARY_LINE([Linker], [$LD $EXPORT_LDFLAGS $CXX_LDFLAGS $LDFLAGS -L${PWD}/src/.libs -lMOAB $PNETCDF_LIBS $NETCDF_LIBS $CGNS_LIBS $HDF5_LIBS $CCMIO_LIBS $CGM_LIBS $ZOLTAN_LIBS $PARMETIS_LIBS $METIS_LIBS $HYPRE_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS])

  # Miscellaneous warnings
  if test "x$WARN_PARALLEL_HDF5" = "xyes"; then
    AC_MSG_WARN([
  *************************************************************************
  *        MOAB has been configured with parallel and HDF5 support
  *     but the configured HDF5 library does not support parallel IO.
  *            Some parallel IO capabilities will be disabled.
  *************************************************************************])
  fi

  if test "x$WARN_PARALLEL_HDF5_NO_COMPLEX" = "xyes"; then
    AC_MSG_WARN([
  *************************************************************************
  *     Your parallel HDF5 library is configured without
  *     H5_MPI_COMPLEX_DERIVED_DATATYPE_WORKS .  For the types of IO
  *     patterns MOAB typically does this will result in degrading
  *     collective IO calls to independent IO, which may have a very
  *     significant impact on IO performance.
  *************************************************************************])
  fi

  MSG_ECHO_SEPARATOR
  BLANK_LINE
  # printf '%-80s' "Please email moab-dev@mcs.anl.gov if you have questions/concerns regarding build and usage."
  MSG_ECHO_CUSTOM([${UNDERLINE_BEG}Build:${UNDERLINE_END} ${BOLD}make all${NORMAL}, ${UNDERLINE_BEG}Test:${UNDERLINE_END} ${BOLD}make check${NORMAL}, ${UNDERLINE_BEG}Install:${UNDERLINE_END} ${BOLD}make install${NORMAL}])

])

AC_DEFUN([BLANK_LINE],
[
  MSG_ECHO_CUSTOM([])
])
dnl ------------------------------------------------------------
dnl  Macro to print one summary line at a time 
dnl  Takes 2 arguments: 
dnl   1) Dependency or feature name
dnl   2) Color-coded text message to print as summary item
dnl ------------------------------------------------------------
AC_DEFUN([SUMMARY_LINE],
[
TRIMMEDDESC=`echo "$2" | sed 's/  */ /g'`
case "$2" in
 "yes") SUMSTRING=`printf '%-8s %-35s : %-4s %-35s\n' "" "${BOLD}$1${NORMAL}" "${GREEN}${TRIMMEDDESC}${NORMAL}" "$3"` ;;
 "no") SUMSTRING=`printf '%-8s %-35s :  %-4s %-35s\n' "" "${BOLD}$1${NORMAL}" "${RED}${TRIMMEDDESC}${NORMAL}" "$3"` ;;
 *) SUMSTRING=`printf '%-8s %-35s :  %-4s %-35s\n' "" "${BOLD}$1${NORMAL}" "${BLUE}${TRIMMEDDESC}${NORMAL}" "$3"` ;;
esac
  AS_ECHO([ "$SUMSTRING" ])
  _AS_ECHO_LOG([[[ SUMMARY ]] --   $SUMSTRING ])
])


