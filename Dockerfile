# Builds a Docker image with all the basic dependencies and
# essentials needed to test MOAB library and
#
# Authors:
# Vijay Mahadevan <vijay.m@gmail.com>
FROM vijaysm/scientific:latest
MAINTAINER Vijay Mahadevan <vijay.m@gmail.com>

ENV HOME="/home/sigma" \
    MOAB_SOURCE_DIR="/home/sigma/moab" \
    MOAB_BUILD_DIR="/home/sigma/moab/build" \
    MOAB_INSTALL_DIR="/opt/moab" \
    TERM="xterm" \
    CC="mpicc" \
    CXX="mpicxx" \
    FC="mpif90" \
    F77="mpif77" \
    F90="mpif90"

# Clone the repository, and configure/build the MOAB sources; 
# Then change permissions so that sigmauser can modify/run examples
RUN git clone https://bitbucket.org/fathomteam/moab.git --depth 1 -b master $MOAB_SOURCE_DIR \
    && cd $MOAB_SOURCE_DIR && autoreconf -fi && mkdir $MOAB_BUILD_DIR && cd $MOAB_BUILD_DIR \
    && ../configure --with-mpi --enable-shared --enable-static --enable-pymoab PYTHON=python3 \
    CPPFLAGS="-Wno-cast-function-type -Wno-maybe-uninitialized -Wno-array-bounds -Wno-implicit-fallthrough -Wno-deprecated-copy -Wno-int-in-bool-context -Wno-unused-but-set-variable" \
    --enable-optimize --enable-debug --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/openmpi \
    --with-metis=/usr FLIBS="-lmpi_cxx" --with-eigen3=/usr/include/eigen3 \
    --prefix=$MOAB_INSTALL_DIR \
    && make all install \
    && cp examples/makefile.config $MOAB_SOURCE_DIR/examples/ \
    && chown -R sigmauser:sigmauser $MOAB_SOURCE_DIR \
    && chown -R sigmauser:sigmauser $MOAB_INSTALL_DIR \
    && ln -s /usr/bin/python3 /usr/bin/python \
    && ln -s /usr/bin/pip3 /usr/bin/pip

#RUN make make MPIEXEC="mpiexec --allow-run-as-root" -C $MOAB_BUILD_DIR check \
#    make -C $MOAB_BUILD_DIR clean \
#    && cd $MOAB_BUILD_DIR/doc && doxygen user.dox

ENV LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$MOAB_INSTALL_DIR/lib" \
    PATH="$PATH:$MOAB_INSTALL_DIR/bin" \
    PYTHONPATH="$PYTHONPATH:$MOAB_INSTALL_DIR/lib/python3.8/site-packages"

# Check if our examples build cleanly as well once environment is concretized.
RUN cd $MOAB_SOURCE_DIR/examples/python && python3 laplaciansmoother.py \
    cd $MOAB_SOURCE_DIR/examples/advanced && make GenLargeMesh && ./GenLargeMesh && mpiexec --allow-run-as-root -n 2 ./GenLargeMesh && rm GenLargeMesh

ENTRYPOINT ["/sbin/setuser","sigmauser","/bin/bash","-c"]
CMD ["/sbin/my_init"]

