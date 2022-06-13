# MOAB_DIR points to top-level install dir, below which MOAB's lib/ and include/ are located
include ../../makefile.config

default: all

ifeq ("$(MOAB_HDF5_ENABLED)","yes")
ALLEXAMPLES += ContinentsOnGlobe BoundaryContinents addPCdata copyPartition AddFieldtoPG2
endif

all: $(ALLEXAMPLES)

ContinentsOnGlobe: ContinentsOnGlobe.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-ContinentsOnGlobe: ContinentsOnGlobe
	${RUNSERIAL} ./ContinentsOnGlobe

BoundaryContinents: BoundaryContinents.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-BoundaryContinents: BoundaryContinents
	${RUNSERIAL} ./BoundaryContinents

addPCdata: addPCdata.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-addPCdata: addPCdata
	${RUNSERIAL} ./addPCdata

copyPartition: copyPartition.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-copyPartition: copyPartition
	${RUNSERIAL} ./copyPartition

AddFieldtoPG2: AddFieldtoPG2.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run-AddFieldtoPG2: AddFieldtoPG2
	${RUNSERIAL} ./AddFieldtoPG2
	
addConnec: addConnec.o ${MOAB_LIBDIR}/libMOAB.la
	@echo "[CXXLD]  $@"
	${VERBOSE}${MOAB_CXX} -o $@ $< ${MOAB_LIBS_LINK}

run: all $(addprefix run-,$(ALLEXAMPLES))

clean: clobber
	rm -rf ${ALLEXAMPLES}

