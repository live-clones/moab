# Builds a Docker image with all the basic dependencies and
# essentials needed to test MOAB library. All TPL have been
# installed using Spack.
#
# Authors:
# Vijay Mahadevan <vijay.m@gmail.com>
# Build stage with Spack pre-installed and ready to be used
FROM spack/ubuntu-bionic:0.16.1 as builder

# What we want to install and how we want to install it
# is specified in a manifest file (spack.yaml)
RUN mkdir /opt/spack-environment \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - netcdf-c+mpi+parallel-netcdf ^mpich" \
&&   echo "  - zoltan+mpi+parmetis ^mpich" \
&&   echo "  - openblas" \
&&   echo "  - eigen" \
&&   echo "  view: /opt/view" \
&&   echo "  concretization: together" \
&&   echo "  config:" \
&&   echo "    install_tree: /opt/software") > /opt/spack-environment/spack.yaml

# Install the software, remove unnecessary deps
RUN cd /opt/spack-environment && spack env activate . && spack install --fail-fast && spack gc -y

# Strip all the binaries
RUN find -L /opt/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-archive\|x-sharedlib' | \
    awk -F: '{print $1}' | xargs strip -s

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . >> /etc/profile.d/z10_spack_environment.sh

# Bare OS image to run the installed executables
FROM ubuntu:18.04
COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software
COPY --from=builder /opt/view /opt/view
COPY --from=builder /etc/profile.d/z10_spack_environment.sh /etc/profile.d/z10_spack_environment.sh

RUN apt -qq update \
    && apt install -y sudo wget pkg-config git \
    python3 python3-pip rsync make cmake \
    autoconf automake libtool \
    clang gcc g++ gfortran vim \
    && apt clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && mkdir -p /home/moabuser/bin /opt/moab \
    && pip3 install cython numpy nose 

ADD bootstrapMOAB /home/moabuser/bin/bootstrapMOAB
ENV LD_LIBRARY_PATH="/opt/view/lib:/opt/view/lib64"
ENV PATH="/opt/view/bin:/opt/spack/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
ENV PKG_CONFIG_PATH="/opt/view/share/pkgconfig:/opt/view/lib/pkgconfig:/opt/view/lib64/pkgconfig"
ENV SPACK_ENV="/opt/spack-environment"

# This makes sure we launch with ENTRYPOINT /bin/bash into the home directory
RUN groupadd -r moabuser -g 433 && useradd -u 431 -r -g moabuser -d /home/moabuser -c "Docker image - MOAB user" moabuser \
    && echo "moabuser:docker" | chpasswd \
    && usermod -a -G sudo moabuser \
    && sh /home/moabuser/bin/bootstrapMOAB \
    && chown -R moabuser:moabuser /opt/moab \
    && chown -R moabuser:moabuser /home/moabuser

USER moabuser
ENV HOME /home/moabuser
ENV BASH_ENV "/home/moabuser/.bash_profile"
WORKDIR /home/moabuser

ENTRYPOINT ["/bin/bash", "-l"]

