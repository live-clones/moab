# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../../makefile.config

default: all

EXAMPLES = spmvApp
ALLEXAMPLES = ${EXAMPLES}

all: $(ALLEXAMPLES)

spmvApp: spmvApp.o netcdf.o ncvalues.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< netcdf.o ncvalues.o ${MOAB_LIBS_LINK}

run-spmvApp: spmvApp
	${RUNSERIAL} ./spmvApp

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

