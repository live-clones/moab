# Builds a Docker image with all the necessary compiler and
# build essentials needed to test SIGMA packages and
# computational toolchains.
#
# Authors:
# Vijay Mahadevan <vijay.m@gmail.com>
FROM vijaysm/scientific:latest
MAINTAINER Vijay Mahadevan <vijay.m@gmail.com>

ENV HOME /home/sigma
ENV MOAB_SOURCE_DIR /home/sigma/moab
ENV MOAB_INSTALL_DIR /opt/moab
ENV CC mpicc
ENV CXX=mpicxx
ENV FC mpif90
ENV F77 mpif77
ENV F90 mpif90
ENV OMPI_ALLOW_RUN_AS_ROOT 1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM 1

# Perform a 
# Set up user so that we do not run as root
# See https://github.com/phusion/baseimage-docker/issues/186
RUN git clone https://bitbucket.org/fathomteam/moab.git --depth 1 -b master $MOAB_SOURCE_DIR \
    && cd $MOAB_SOURCE_DIR && autoreconf -fi && mkdir build && cd build \
    && ../configure --with-mpi --enable-shared --enable-static --enable-pymoab PYTHON=python3 \
    --enable-optimize --enable-debug --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/openmpi \
    --with-metis=/usr FLIBS="-lmpi_cxx" --with-eigen3=/usr/include/eigen3 \
    --prefix=$MOAB_INSTALL_DIR \
    && make all install check clean \
    && cp examples/makefile.config $MOAB_SOURCE_DIR/examples/ \
    && cd doc && doxygen user.dox

ENTRYPOINT ["/sbin/setuser","sigmauser","/bin/bash","-c"]
CMD ["/sbin/my_init"]

