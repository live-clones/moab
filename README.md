
# MOAB: Mesh-Oriented datABase

MOAB is a component for representing and evaluating mesh data. MOAB can store structured and unstructured mesh, consisting of elements in the finite element "zoo". The functional interface to MOAB is simple yet powerful, allowing the representation of many types of metadata commonly found on the mesh. MOAB is optimized for efficiency in space and time, based on access to mesh in chunks rather than through individual entities, while also versatile enough to support individual entity access. MOAB can be used in several ways: 

* As the underlying mesh data representation for applications
  - Several computational solvers in various scientific domains (nuclear engineering, climate modeling, nonlinear thermo-mechanics, CFD, etc)
  - Mesh smoothing and optimization
  - Mesh quality computation
* As a mesh input mechanism (using mesh readers included with MOAB),
* As as a translator between mesh formats (using readers and writers included with MOAB).

MOAB was developed originally as part of the CUBIT project at Sandia National Laboratories, and has been partially funded by the DOE SciDAC program (TSTT, ITAPS, FASTMath) and DOE-NE (NEAMS program).
More recently, DOE-BER programs under the E3SM project have provided support for enabling scalable solution transfer techniques for climate applications.

## Dependencies

- **MPI**: MOAB supports usage of MPICH and OpenMPI libraries configured externally in order to enable scalable mesh manipulation algorithms.
- **HDF5**: In order to manage the data dependencies and to natively support parallel I/O, MOAB uses a custom file format that can represent the entire MOAB data model in a native HDF5-based file format. Support for this file format requires version 5 of the HDF library, which can be obtained at [HDF5].
- **NetCDF**: MOAB library optionally depends on the NetCDF libraries (C and C++) to compile the ExodusII reader/writer. To get netcdf, go to [NetCDF].
- **Metis**: MOAB can use the Metis library for partitioning mesh files in serial
- **Zoltan**: Support for online partitioning through Zoltan (and its dependencies on Scotch, Parmetis etc) can be utilized through the partitioner tool
- **TempestRemap**: Provide support for both offline and online remapping of Climate field data on unstructured spherical meshes
- **Eigen3**: A substitute for BLAS/LAPACK interfaces. However if *TempestRemap* tools are to be built, this becomes a required dependency

## Configuration and Build

* Currently, both CMake and Autotools are maintained simultaneously in order to support all platforms (including Windows). Please choose your build system according to your needs and follow instructions below.

* **Autotools based configuration workflow**

  - Please ensure that the autotools toolchain is pre-installed locally. We recommend a minimum autoconf version of *v2.69*.
  - Run `autoreconf -fi` to generate the configure script
  - Run the `configure --help` script in the top source directory to see a list of available configure options.
      - Use `--prefix=INSTALL_DIR` to specify installation directory
      - Override default compilers with environment or user options: `CC, CXX, FC, F77`
      - If you have **MPI** installed, use `--with-mpi=$MPI_DIR`
      - If you have **HDF5** and **NetCDF** installed, use `--with-hdf5=$HDF5_DIR`  and `--with-netcdf=$NETCDF_DIR` to specify external dependencies.
      -  Similarly for **Metis** or **ParMetis** dependencies, use `--with-metis=$METIS_DIR` and `--with-parmetis=$PARMETIS_DIR` respectively
      -  If you have **Zoltan** installed, use `--with-zoltan=$ZOLTAN_DIR`
      - In order to enable Python bindings (**PyMOAB**), you need to use `--enable-shared --enable-pymoab` options
      - **Auto-Download options**: MOAB now supports automatic dependency download and configuration that has been tested on various platforms and architectures.
        + *HDF5*: Use `--download-hdf5` OR `--download-hdf5=TARBALL_PATH`
        + *NetCDF*: Use `--download-netcdf` OR `--download-netcdf=TARBALL_PATH`
        + *Metis*: Use `--download-metis` OR `--download-metis=TARBALL_PATH`
        + *TempestRemap*: Use `--download-tempestremap` OR `--download-tempestremap=master` (to build from Git master)
  - Now run the `configure` script with desired configuration options either in-source or out-of-source (build) directory.

* **CMake based configuration workflow**

  - Please ensure you have CMake (>3.0) available locally.
  - Run `ccmake` visual configuration editor or `cmake` to get a bare configuration of MOAB
  - If you would like to override the compiler used by default, use the following variables to override
    - C: `-DCMAKE_C_COMPILER=mpicc`
    - C++: `-DCMAKE_CXX_COMPILER=mpicxx`
    - Fortran: `-DCMAKE_Fortran_COMPILER=mpif90`
  - If you have **MPI** installed, use `-DENABLE_MPI=ON -DMPI_HOME=$MPI_DIR`
  - If you have **HDF5** installed, use `-DENABLE_HDF5=ON -DHDF5_ROOT=$HDF5_DIR`
  - If you have **NetCDF** installed, use `-DENABLE_NETCDF=ON -DNETCDF_ROOT=$NETCDF_DIR`.
  - Similarly for **Metis** or **ParMetis** dependencies, use `-DENABLE_METIS=ON -DMETIS_DIR=$METIS_DIR` and `-DENABLE_PARMETIS=ON -DPARMETIS_DIR=$PARMETIS_DIR` respectively
  -  In order to enable Python bindings (**PyMOAB**), you need to use `-DBUILD_SHARED_LIBS=ON -DENABLE_PYMOAB=ON` options
  - If you have **Zoltan** installed, use `-DENABLE_ZOLTAN=ON -DZOLTAN_DIR=$ZOLTAN_DIR` options

* Once configuration with autotools or CMake is complete in the build directory, run the following to build the library:
  - Compile MOAB and supported tools: `make -j4`
  - Verify configuration and build setup: `make check`
* Next to install the compiled libraries, headers and tools, run `make install`
* You can now use the `makefile` generated under the `build/examples` folder and modify it to compile downstream code with MOAB dependency

## Continuous Integration

There are several hooks to online continuous integration systems, nightly and commit-based Buildbot/Bitbucket builds that are constantly run during a development day to check the integrity and robustness of the MOAB library.

### Current overall build status

- **Buildbot**: [ ![Buildbot Status](http://gnep.mcs.anl.gov:8010/badges/moab-all.svg)](https://gnep.mcs.anl.gov:8010)
- **CircleCI**: [ ![CircleCI Status](https://circleci.com/bb/fathomteam/moab/tree/master.svg?style=shield)](https://circleci.com/bb/fathomteam/moab)
- **CodeShip**: [ ![Codeship Status](https://codeship.com/projects/286b0e80-5715-0132-1105-0e0cfcc5dfb4/status?branch=master)](https://codeship.com/projects/49743)
- **Code Coverage**: 
  + **Coverity**: [ ![Coverity Scan Build Status](https://scan.coverity.com/projects/6201/badge.svg)](https://scan.coverity.com/projects/moab)
  + **CodeCov**: [![codecov](https://codecov.io/bb/fathomteam/moab/branch/master/graph/badge.svg)](https://codecov.io/bb/fathomteam/moab)

## Bugs, Correspondence, Contributing
MOAB is distributed under LGPL(v3) licensing, and we encourage users to submit bug reports (and, if possible, fixes) directly on the Bitbucket interface. Optionally, users can email and discuss the issue with developers at [moab-dev@mcs.anl.gov](mailto:moab-dev@mcs.anl.gov). Please submit your changes through a pull requests (PR) with a Bitbucket fork (refer to [CONTRIBUTING.md](CONTRIBUTING.md) for further information), or send us patches that you would like merged upstream. Users are also encouraged to check [SIGMA-MOAB] documentation pages for news and updates. 

[NetCDF]: http://www.unidata.ucar.edu/software/netcdf/
[HDF5]: https://www.hdfgroup.org/HDF5/
[SIGMA-MOAB]: http://sigma.mcs.anl.gov/moab-library
