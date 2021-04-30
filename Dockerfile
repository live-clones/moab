# Builds a Docker image with all the basic dependencies and
# essentials needed to test MOAB library. All TPL have been
# installed using Spack.
#
# Authors:
# Vijay Mahadevan <vijay.m@gmail.com>
FROM sigmaanl/moab:spack
MAINTAINER Vijay Mahadevan <vijay.m@gmail.com>

ENV MOAB_SOURCE_DIR="/home/moabuser/moab" \
    MOAB_BUILD_DIR="/home/moabuser/moab/build" \
    MOAB_INSTALL_DIR="/opt/moab" \
    CC="mpicc" \
    CXX="mpicxx" \
    FC="mpif90" \
    F77="mpif77" \
    F90="mpif90"

ENTRYPOINT ["/bin/bash", "-l"]
