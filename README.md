# MOAB: Mesh-Oriented datABase

MOAB is a software component for representing for creating, storing and accessing mesh data. MOAB can describe structured and unstructured mesh, consisting of elements in the finite element "zoo". The functional interface to MOAB is simple yet powerful, allowing the representation of many types of metadata commonly found on the mesh. MOAB is optimized for efficiency in space and time, based on access to mesh in chunks rather than through individual entities, while also versatile enough to support individual entity access. MOAB can be used in several ways: 

A few highlights of the capabilities in MOAB include:

* Representation of 0-3d elements in the finite element "zoo" (including support for quadratic elements), as well as support for polygon and polyhedron entities
* Highly efficient storage and query of structured and unstructured mesh (e.g. a brick-shaped hex mesh requires approximately 25 and 55 MB per million hex elements in the structured and unstructured representations, respectively)
* Powerful data model allowing representation of various metadata in the form of "sets" (arbitrary groupings of mesh entities and sets) and "tags" (annotation of entities, sets, and entire mesh)
* Open source (LGPL) mesh readers/writers for Sandia ExodusII, CUBIT .cub save/restore, VTK, GMsh, and other mesh formats with capability to translate between them uniformly
* Flexible access to MOAB routines from C and Fortran through *iMesh* and *iMOAB* interfaces are available
* A high level Python interface (*PyMOAB*) based on Cython bindings can also be enabled

Several computational solvers in various scientific domains such as nuclear engineering, climate modeling, nonlinear thermo-mechanics, CFD, etc have been built on top of MOAB. Other common use-cases where MOAB is often applied are:

  - Unstructured mesh generation and manipulation for complex geometry
  - Mesh quality computation along with algorithms for smoothing and optimization
  - Solution field transfers for multi-physics problems
    + High-order interpolation between unstructured grids in two and three dimensions 
    + Conservative remappng between meshes on the sphere for Climate problems

MOAB was developed originally as part of the CUBIT project at Sandia National Laboratories, and has been partially funded by the DOE SciDAC program (TSTT, ITAPS, FASTMath), ASCR (CESAR), and DOE-NE (NEAMS program). More recently, DOE-BER programs under the E3SM project have provided support for enabling scalable solution transfer techniques for climate applications.

MOAB is distributed under an open-source, GNU LGPL licensing agreement. [![LGPL-version3 License](https://img.shields.io/badge/license-LGPLv3-green)](LICENSE)

## Continuous Integration

There are several hooks to online continuous integration systems, nightly and commit-based Buildbot/Bitbucket builds that are constantly run during a development day to check the integrity and robustness of the library.

[![Buildbot Status](http://gnep.mcs.anl.gov:8010/badges/moab-all.svg)](https://gnep.mcs.anl.gov:8010)
[![Bitbucket Pipelines Status](https://img.shields.io/bitbucket/pipelines/fathomteam/moab/master?label=bitbucket&style=plastic)](https://bitbucket.org/fathomteam/moab/addon/pipelines/home#!/results/branch/master/page/1)
[![CircleCI Status](https://circleci.com/bb/fathomteam/moab/tree/master.svg?style=shield)](https://circleci.com/bb/fathomteam/moab) 
[![Codeship Status](https://app.codeship.com/projects/d84b8203-8cf9-4763-ae59-24e623ca5d19/status?branch=master)](https://app.codeship.com/projects/459370)
[![Buddy Pipeline](https://app.buddy.works/vijaym/moab/pipelines/pipeline/398746/badge.svg?token=6b9e040bcd31cd1d8d2f409795f7896548b8604925af3de730aea2499d59520c "buddy pipeline")](https://app.buddy.works/vijaym/moab/pipelines/pipeline/398746)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/6201/badge.svg)](https://scan.coverity.com/projects/moab)
[![codcov](https://codecov.io/bb/fathomteam/moab/branch/master/graph/badge.svg)](https://codecov.io/bb/fathomteam/moab)
[![Codacy grade](https://app.codacy.com/project/badge/Grade/a796f9d9f5f44d628de15ab95717b1d1)](https://www.codacy.com/bb/fathomteam/moab/dashboard?utm_source=vijaysm@bitbucket.org&amp;utm_medium=referral&amp;utm_content=fathomteam/moab&amp;utm_campaign=Badge_Grade)

## Documentation

Detailed API documentation and user/development guides are available for the following repository branches, updated daily.

  - [master](https://ftp.mcs.anl.gov/pub/fathom/moab-docs/index.html)
  - [develop](https://ftp.mcs.anl.gov/pub/fathom/moab-docs-develop/index.html)

## MOAB Pre-installed

  - MOAB pre-installed docker image: [![Docker for MOAB](https://img.shields.io/docker/pulls/vijaysm/moab-root?style=flat-square)](https://hub.docker.com/repository/docker/vijaysm/moab-root)
  - MOAB tools with anaconda: [![conda-forge MOAB](https://img.shields.io/conda/pn/conda-forge/moab?color=blue&label=conda&style=plastic)](https://anaconda.org/conda-forge/moab)

## Optional Dependencies

- **MPI**: MOAB supports usage of MPICH and OpenMPI libraries configured externally in order to enable scalable mesh manipulation algorithms.
- **HDF5**: In order to manage the data dependencies and to natively support parallel I/O, MOAB uses a custom file format that can represent the entire MOAB data model in a native HDF5-based file format. Support for this file format requires version 5 of the HDF library, which can be obtained at [HDF5].
- **NetCDF**: MOAB library optionally depends on the NetCDF libraries (C and C++) to compile the ExodusII reader/writer. To get netcdf, go to [NetCDF].
- **Metis**/**ParMetis**: MOAB can use the Metis or ParMetis library for partitioning mesh files in serial and parallel respectively
- **Zoltan**: Support for online partitioning through Zoltan (and its dependencies on Scotch, ParMetis etc) can be utilized through the partitioner tool
- **TempestRemap**: Provide support for both offline and online remapping of Climate field data on unstructured spherical meshes
- **Eigen3**: A substitute for BLAS/LAPACK interfaces. However if *TempestRemap* tools are to be built, this becomes a required dependency

## Configuration and Build from Source

* Currently, both CMake and Autotools are maintained simultaneously in order to support all platforms (including Windows). Please choose your build system according to your needs and follow instructions below. Both of these workflows follow the same pattern of commands to build and install in your platform.

### **Autotools based configuration workflow**

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

### **CMake based configuration workflow**

  - Please ensure you have CMake (>3.0) available locally.
  - Run `ccmake` visual configuration editor or `cmake` to get a bare configuration of MOAB
  - If you would like to override the compiler used by default, use the following variables to override
    - C: `-DCMAKE_C_COMPILER=mpicc`
    - C++: `-DCMAKE_CXX_COMPILER=mpicxx`
    - Fortran: `-DCMAKE_Fortran_COMPILER=mpif90`
  - Specify your installation directory by using `-DCMAKE_INSTALL_PREFIX=$MOAB_INSTALL_PATH`
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

## Language Bindings

Even though the MOAB library is written in C++ language (conforming to C++11 standard), several partial bindings and interfaces are available for other languages.

  + **C/Fortran**: You can use the iMOAB interface to load, manipulate and query unstructured meshes in memory
    - Supports both serial and parallel invocation under one interface. MOAB needs to be configured using `--with-mpi` option.
    - Supports ability to migrate meshes and tags between processes
    - Supports capability to compute remapping weights for Climate science applications

    Note that the **ITAPS iMesh** interfaces are now deprecated and we encourage users to directly use MOAB or iMOAB interfaces if possible.
    
  + **Python3**: The Python bindings for MOAB can be enabled when configured with `--enable-shared --enable-pymoab` options.
    - Supports access to the structured grid interfaces.
    - Supports queries and access to Entities using a true-Pythonic implementation flavor.
    - Utilizes Cython to provide flexible bindings without sacrificing runtime performance.
    - Only supports serial computations for now. Parallel implementation using MPI4Py is underway.

### Third-party bindings

  + **C#**: This contributed source is developed and maintained by Qingfeng Xia, UKAEA 2021. It is distributed under the same LGPLv3 license as MOAB.
    - The open-source repository containing C# interfaces are available in [MOABSharp Bitbucket repository](https://bitbucket.org/qingfengxia/moabsharp).
    - Supports builds on .NET Core 3.1 under Ubuntu 20.04, and .NET framework 4.x under Windows 10-64bit (target on "netstandard2").
    - Detailed documentation and supported features are listed in the [MOABSharp README file](https://bitbucket.org/qingfengxia/moabsharp/src/dev/README.md).

## Bugs, Correspondence, Contributing

MOAB is distributed under LGPL(v3) licensing, and we encourage users to submit bug reports (and, if possible, fixes) directly on the Bitbucket interface. Optionally, users can email and discuss the issue with developers at [moab-dev@mcs.anl.gov](mailto:moab-dev@mcs.anl.gov). Also refer to [FAQ](FAQ.md) for some commonly asked questions and their resolutions. 

MOAB follows a fully transparent, open-source development workflow (refer to our [Code of Conduct](CODE_OF_CONDUCT.md) for further information), and we welcome all contributions that enhance the feature sets provided by MOAB. If you would like to contribute to MOAB, please submit your changes through pull request (PR) using a Bitbucket fork of MOAB (refer to [CONTRIBUTING.md](CONTRIBUTING.md) for more details), or send us patches that you would like merged upstream. Users are also encouraged to check [SIGMA-MOAB] documentation pages for news and updates. 

## Citing MOAB

If you use MOAB for your research, please use the following bibtex entries for the software and the original report to cite us.

```
@techreport{moab_2004,
  author = {Tautges, T. J. and Meyers, R. and Merkley, K. and Stimpson, C. and Ernst, C.},
  type = {{SAND2004-1592}},
  title = {{MOAB:} A Mesh-Oriented Database},
  institution = {Sandia National Laboratories},
  month = apr,
  year = {2004},
  note = {Report}
}

@software{moab521_2020,
  author       = {Mahadevan, Vijay and
                  Grindeanu, Iulian and
                  Jain, Rajeev and
                  Shriwise, Patrick and
                  Wilson, Paul},
  title        = {MOAB v5.2.1},
  month        = aug,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {5.2.1},
  doi          = {10.5281/zenodo.2584862},
  url          = {https://doi.org/10.5281/zenodo.2584862}
}
```

Additionally, if you would like us to highlight your application, library or tool that uses MOAB, please contact the developers for more information.

[NetCDF]: http://www.unidata.ucar.edu/software/netcdf/
[HDF5]: https://www.hdfgroup.org/HDF5/
[SIGMA-MOAB]: http://sigma.mcs.anl.gov/moab-library

