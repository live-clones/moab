# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../../makefile.config

default: all

EXAMPLES = kdtreeApp
ALLEXAMPLES = ${EXAMPLES}

all: $(ALLEXAMPLES)

kdtreeApp: kdtreeApp.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK} 

run-kdtreeApp: kdtreeApp
	${RUNSERIAL} ./kdtreeApp

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

