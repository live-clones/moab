# MOAB Library Release Notes

## Version 5.3.1

### Features

- PR #563: A new tool to visualize maps (mbvisumap) to display DoF coupling between source and target component meshes with a given linear map file. 

### Enhancements

- PR #549: Performance improvement for Eigen3 SparseMatrix insertion
- PR #554: Read map in parallel and apply it on existing field
- PR #561: Fallback to PNetCDF for reading map files in parallel, when NetCDF4/HDF5-parallel interface is unavailable
- PR #566: Consistently decompose polygons in parallel so that the order of operations is preserved. This fix ensures results between parallel mbtempest and serial TempestRemap runs return the same map files
- PR #568: Adding a mirror for TPL tarballs now at ANL FTP site
- PR #569: Configuration fixes for Eigen3 (CMake) and other misc items for autotools

### Fixes

- PR #553: Fix for building on MinGW64
- PR #556: Minor fix for merging by integer tag
- PR #558, #560: iMOAB errors on master due to rebase errors
- PR #559: Fix for iMOAB coupler test and read map test on one task
- PR #565: Fix for correctly computing coordinates of triangle at poles

## Version 5.3.0

### Features

- *PR #502*: Gnomonic projection option for partitioning using an <ins>"inferred"</ins> strategy with Zoltan RCB algorithm for meshes on a sphere
- *PR #513*: New example to produce the dual of a delaunay mesh on a sphere to generate a polygonal grid (MPAS-like)
- *PR #525, PR #526*: CONTRIBUTING.md: Documentation and guidelines for new users contributing to MOAB
- *PR #535*: Docker image of MOAB pre-installed that uses <ins>spack</ins> for all TPL dependencies
    - *PR #551*: Separate docker image for MOAB+TempestRemap for offline climate remappiing workflows
- *PR #544*: MOAB now requires a minimum of **C++11** standard to be supported by the C++ compiler

### Enhancements

- *PR #503*: Download git branches for dependencies directly instead of just tarball. Typical usage: `--download-<package>=git:branchname`
- *PR #523, PR #528*: CI system updates for BitBucket-Pipelines and Shippable
- *PR #531*: Update README
- *PR #530*: Add pkg-config support
- *PR #538*: Print configuration summary after Autotools configure and CMake commands
- *PR #543*: Autodownload a compatible version of Eigen3 header-only library

### Fixes

- *PR #494*: Fully support Python3 and deprecate Python2 usage
- *PR #497*: Factor out repetitive code in many iMOAB tests for Climate remap
- *PR #498*: Improve efficiency of merging for large meshes using MergeMesh
- *PR #504*: Fixing off-by-one error in the check for resizing the root sets data structures in GeomTopoTool
- *PR #505*: Remove -march=native for optimized mode
- *PR #506*: Distclean now removes any orphan a.out files as well
- *PR #507*: Fix sparse tag writing in HDF5 format, when buffer size exceeds first allocated amount (bug in buffered write)
- *PR #508*: Structured mesh output in VTK format - use structured layout rather than serializing explicitly
- *PR #510, PR #517*: Move the handling of the SKIP\_TOPOLOGY option to a more appropriate location in TQDCFR
- *PR #511*: Replace MPI\_INTEGER with MPI\_INT in couple of places that was causing segfaults
- *PR #512, PR #516, PR #518, PR #520, PR #532*: Autotools fixes for general configuration and auto-download of dependencies
- *PR #521, PR #522, PR #524, PR #534, PR #540, PR #542*: Windows build fixes with CMake workflow using MinGW and VS2015 (or higher)
- *PR #529*: Fix to ensure thread-safety in MOAB caused due to a cached variable in TypeSequenceManager::find
- *PR #536*: Several key fixes using clang-tidy for the entire repository
- *PR #545*: Fix MOAB CGNS I/O failure when the underlying CGNS library is configured with `CGNS\_ENABLE\_64BIT` option
- *PR #547*: Some updates for 32-bit builds; This support will get deprecated and dropped in the next major release.
- *PR #550*: Removing deprecated references to std::unary_function and std::binary_function to allow build with C++17 compiler family

### Application Support

- *PR #499*: Capability to read climate SCRIP files in parallel through buffered read mechanism
- *PR #500*: Support instancing LND domain files directly
- *PR #501*: Support both fully parallel and serial builds with TempestRemap
- *PR #509*: Introduce a new tool to compare maps in NetCDF (nc) output files
- *PR #514, PR #519, PR #527*: Capability to read climate SCRIP files in parallel directly (with MPI+NetCDF/HDF5 interfaces)
- *PR #533*: Introduce baseline tests to catch regressions when computing field projections from ATM to OCN component
- *PR #546*: Consolidated set of several fixes for both the offline and online remapping weight generation workflows in E3SM.
    - Updates to ensure that the Kd-tree tolerance computation is more robust, especially for RLL meshes
    - `mbtempest` will now use the same tolerance for intersection as TempestRemap usage for Node tolerances
    - Fix bugs in parallel I/O to get areas, fractions, mesh coordinates and other parameters ordered correctly
    - Add modifications to the mesh reading mechanisms to adapt I/O based on underlying mesh type

## Version 5.2.1

### Fixes

- *PR #492*: Use a tighter boxeps tolerance in Kd-tree searches for performance optimization
- *PR #493*: make check fixes for parallel tests using OpenMPI in Conda
- *PR #489*: New test for 3-component 2-way coupling with the PG2 ATM-OCN-LND cases
- *PR #495*: Minor fix in PyMOAB get\_parent\_meshsets implementation

## Version 5.2.0

### Features

- *PR #424, PR #445, PR #461, PR #462, PR #467, PR #474*: Significant feature additions in `mbtempest` to improve the TempestRemap-based mapping weights generation workflow
    - Expose both Kd-tree and advancing front intersection
    - Verifications of mapping weight generation for very high-res cases at scale
    - Improve `mbconvert` to transform TR supported file formats to MOAB h5m and vice-versa
    - Ability to write out the computed remapping weights in parallel to a h5m file is now supported. A `h5mtoscrip` utility tool to convert the h5m file to SCRIP format is available.
    - Add an extra verification option for maps computed with `mbtempest` to validate projection weights
    - Improve master-slave partitioner tool
    - Performance enhancements: Several optimizations in the iMOAB implementation and OnlineMap computation routines have been implemented
- *PR #437*: Allow ability to support empty parts in parallel mesh handling
- *PR #454, PR #464, PR #475*: Enable coupling and consistent remap of data between discrete meshes and point cloud representations
    - A new ParCommGraph class is now available to store communication patterns between MPI groups
    - Feature additions support rendezvous coupling that has been tested for multi-way transfer between atmosphere, ocean, land components
- *PR #465*: Add a new KDtree-based intersection to replace advancing front intersection when dealing with source grids with holes
- *PR #471*: Introduce a concept of master-slave partitioning method using Zoltan RCB algorithm
    - This feature provided a significant speedup for offline remapping calculations on a sphere
    - *TODO*: Test similar ideas for other coupled computational problems
- *PR #482*: Domain file reader in NetCDF format to support E3SM integration
- *PR #485*: Extensive CMake configuration updates to support Windows platform builds
- *PR #486*: Adding a new option for area computations using Gauss quadrature based integration instead of spherical excess

### Enhancements

- *PR #419*: Cache the latest facet for RayHistory in OBB tree
- *PR #420, PR #435*: Global ID default value (unset) is now -1, while DOFs related to Global ID start from 1
- *PR #429*: Option for consistent (multiple) ghost layers in thin partitions
- *PR #438*: Improve coexistence of LAPACK and Eigen3 in configurations
- *PR #442*: MBTR - Remove assumption that FV DoF numbering is 0-based (now 1-based for both FV/SE)
- *PR #444, PR #449*: CMake build improvements related to MPI, HDF5, OSX
- *PR #446*: Option to reset Global IDs during partitioning
- *PR #455*: Enforce MOAB flags to add c++0x flags if supported
- *PR #463*: Update default download option for NetCDF(4.6.3), HDF5 (1.10.1)
- *PR #468*: PyMOAB bindings to set\_connectivity
- *PR #475*: Improve regridding intersection robustness
- *PR #483*: Add a new download option to sync large files from FTP site to avoid bloating repository with static meshes
- *PR #488*: Enforce uniform formatting of the entire source repository with clang-format tool with updated coding guidelines

### Fixes

- *PR #414, PR #434*: HDF5 writer fixes for polyhedral meshes, to address issues #99 and #110
- *PR #415*: Fixing index overflow in SQIJ partition for structured meshes
- *PR #417, PR #421, PR #469*: Python bindings fixes and compliance to Python 3 standards
- *PR #430*: Fixing HDF5 regression test issues
- *PR #431, PR #453, PR #466, PR #480*: Changes to conform to evolving TempestRemap API co-development
- *PR #450*: Fix seqfault in Range destructor
- *PR #472*: HDF5 file descriptor initialization that resolves some issues in VisIt-MOAB plugin
- *PR #473*: Minor fix for Tqdcfr reader

## Version 5.1.0

### Features

- *PR #399, PR #406, PR #408, PR #410*: Updates for E3SM workflow to support online remapping computation for climate applications through interfaces to TempestRemap library.
    - Major additions to expose online remapping capabilities by computing exact intersection between unstructured meshes on a sphere, and providing projection weight operator computation through interfaces to TempestRemap for both scalar and vector solution fields with conservation and/or monotonicity constraints.
    - The primary application target is E3SM for this particular feature addition. However, the bulk of the interfaces added in this PR can be reused in other coupled physics simulations as well.
- *PR #403*: Updates to Cython file generation process in PyMOAB.
    We've done away with static `.cpp` files in the project and are now relying on users having a Cython installation greater than 0.26 on their machine.  This helps us with python 2/3 compatibility issues as well as superfluous re-compiles of libraries.

### Enhancements

- *PR #376, PR #365*: Updates to PyMOAB to improve the Cython interaction and Python3 compatibility.
- *PR #397*: Consistent CMake configurations for downstream apps
- *PR #404*: Update RTT reader to support v1.0.1 specs
- *PR #391, PR #401*: Delete volume/surface OBB trees and restore OBBs from file
- *PR #370*: Several improvements in the GeomTopoTool implementation

### Fixes

- *PR #363, PR #398*: PyMOAB tag and skinner test fixes
- *PR #364*: CMake configuration updates to not always require Fortran for BLAS/LAPACK linkage
- *PR #362*: Several semantic and memory leak fixes due to Coverity, cppcheck and valgrind
- *PR #364, PR #366, PR #368, PR #373, PR #374, PR #375, PR #379, PR #381, PR #382, PR #384, PR #387, PR #388, PR #389, PR #400, PR #407, PR #409*: Several miscellaneous configuration updates for both autotools and CMake workflows

## Version 5.0.0

### Features

- *PR #248, PR #290, PR #299*: Eigen3, BLAS/LAPACK support; Now BLAS/LAPACK are required dependencies for MOAB in order to support high-performance dense solver support for various algorithms
- *PR #278, PR #286, PR #280, PR #296, PR #300, PR #325, PR #326, PR #332, PR #331, PR #334, PR #335, PR #317, PR #345, PR #346, PR #353, PR #356*: PyMOAB interface now directly exposes Python modules to access common MOAB functionality in order to load, query and manipulate structured/unstructured meshes in memory. It also includes interfaces to some specific classes such as the ScdInterface and its associated classes along with the TopoUtil tool. We have also added support for Tags, Range, Skinner modules that can be used to build more advanced capabilities as needed. Currently, this interface only supports serial configuration and we intend to support parallelism through `mpi4py` dependencies in the next release. Several useful tests and examples have also been added for users to get started.
- *PR #162*: Add capability to load an "*unpartitioned*" mesh file through a trivial partition with the parallel option: `PARALLEL=BCAST_DELETE;PARTITION=TRIVIAL`; This works for any mesh format and so we are no longer restricted to just our parallel h5m native mesh.
- *PR #238*: Discrete Geometry Module -- Higher order curve and surface reconstruction strategies based on local polynomial fittings along with integrated support for recursive uniform refinement.
- *PR #318*: Add Mesquite optimization support to MOAB so that mesh optimization algorithms can be natively run on a MOAB mesh. It can also use the geometry information when MOAB is configured with CGM (through Lasso relational databases in MOAB).
- *PR #330*: Add support for computing high-order, conservative, multi-mesh projection weights with the **TempestRemap** library. Several optimizations have been implemented for the advancing-front algorithm to compute the intersections between two unstructured meshes on a sphere. In addition to this, the integration with Tempest library enables computation of the projection weights between FV/cGLL/dGLL discretizations seamlessly in a scalable fashion, even though Tempest itself is a serial library!
- *PR #288*: HEX8 (linear) to HEX27 (quadratic) mesh conversion utility

### Enhancements

- *PR #263*: Support for building with CGM build directory (instead of installation) to ease SIGMA workflows
- *PR #267*: Modify scale and reflect signature in FBiGeom, now it conforms to iGeom
- *PR #266*: Update Matrix3 constructor
- *PR #143*: The original design of MOAB Core was responsible for initializing MPI (if the user / driver forgot to do it in a parallel setting). However, this is just bad design since after destruction of moab::Core object, and hence call to MPI_Finalize(), no other MPI calls could be made. This design has now been fixed and MPI\_Init and MPI\_Finalize should now be called in user driver code.
- *PR #287*: Code fails gracefully when number of ranks are greater than number of parts in loaded h5m file
- *PR #276*: iMOAB extensions to create and manipulate vertices in addition to querying entities and tags on the mesh
- *PR #291, PR #302*: The autotools parallel build toolchain has been accelerated and simplified
- *PR #298, PR #316*: Read node materials from gmsh file and other I/O issues.
- *PR #285, PR #293, PR #303, PR #306, PR #323, PR #315, PR #327, PR #339, PR #348, PR #355*: Several DagMC related enhancements
    - A major DagMC refactor to move bulk of the capability to OrientedBoxTreeTool, GeomTopoTool, and GeomQueryTool tools
    - various performance and API improvements
    - Completely remove the DagMC tool from MOAB since the bulk of the capabilities are now exposed through other native classes in MOAB directly.
    - Implement context for determining whether ray intersections are counted
    - Add a method to GeomTopoTool to recover the topology of a geometry of nested volumes when the file does not contain this information
- *PR #311*: Update Debian dependencies.
- *PR #308*: `mbconvert` tool can now include the MPAS (graph) partition data during conversion
- *PR #304*: AHF adjacency queries will not be turned on by default until we add support for mixed meshes and poly-meshes
- *PR #201, PR #329*: The examples folder has been fully re-organized to better help users get started with capabilities and understand usage of MOAB for computational use-cases. New examples have also been added.
- *PR #344*: Adding an example to compute exact normals using iRel when a geometry is available.
- *PR #343*: Refactor the GenLargeMesh example to a callable API to generate structured HEX/TET meshes in memory
- *PR #341*: Update the configuration suggestions for LCF machines: ALCF, NERSC
- *PR #359*: Adding an option to explicitly specify NetCDF-C++ interfaces when TempestRemap is enabled

### Fixes

- *PR #261*: Removing several verbose error outputs
- *PR #264*: Fixes for 32-bit architectures
- *PR #250, PR #262, PR #282*: Better support for PGI and Clang compilers
- *PR #268*: Titan cray intel issue for parallel HDF5 I/O
- *PR #250, PR #277, PR #324*: MPAS and NetCDF I/O enhancements and fixes
- *PR #270*: Better tests for ray tracing and OBB
- *PR #273, PR #310, PR #312, PR #322, PR #347*: Several updates and fixes for autotools and CMake build configuration
- *PR #275, PR #297*: Several memory leak fixes for MOAB algorithms
- *PR #283*: Fix for DAGMC looking for OBB\_TAG
- *PR #283*: Several updates to `mbcoupler` tool and adding support for moab::SphericalQuad elements
- *PR #313*: When reading Exodus files, if the entities are higher order, use only the corner nodes for matching
- *PR #320*: Update options to disable fortran (this was complicated due to introduction of BLAS/LAPACK libraries as required dependencies)
- *PR #309, PR #333, PR #354*: Documentation updates
- *PR #349*: Update HDF5 includes - found through failures in Ubuntu installations
- *PR #337*: Fixing several compiler (GNU, Clang, Intel) warnings

## Version 4.9.2

### Features

- *PR #259*: Introduce the new language-agnostic iMOAB interface to MOAB that is oriented towards FEM/FDM/FVM applications codes.
    - This is a ongoing replacement to the ITAPS iMesh interfaces but with much less verbosity and indirection in implementation.
- *PR #196*: Expose "mhdf" publically to be consumed by advanced visualization (VisIt) and parallel I/O use-cases
- *PR #195*: Addition of an optimized parallel mesh resolution algorithm when performing UMR in memory
- *PR #193*: Support for Attila RTT finite element mesh files
- *PR #203*: Support for the general OBJ triangulation files
- Add support for Bitbucket pipelines, Codeship/Drone.io integration for improving CI and testing infrastructure on a per-commit basis

### Enhancements

- *PR #229*: Add support for polygon and polyhedra in VTK and NetCDF formats
- *PR #246*: Improved re-design of the DagMC tool
- *PR #251*: Support 64-bit integers with MetisPartitioner
- *PR #247*: Better support for configuration/build on BG/Q systems (ANL Vesta/Mira with GNU/XLC compilers)
- *PR #258*: Redesign of options to control sequence allocation (previously CoreOptions)

### Fixes

- *PR #234*: Several key fixes found through cppcheck and Coverity (static analysis tools)
- *PR #233*: Parallel h5m file load of recursive sets
- *PR #213, PR #230, PR #255, PR #253*: Autotools/CMake fixes for MPI, HDF5, NetCDF configurations on various architectures
- *PR #231*: Remove explicit dependence on Cubit through MOAB
- *PR #227*: Several fixes for reading NetCDF and MPAS climate files
- *PR #236*: Correctly process mesh files created on Windows architectures
- *PR #242*: Make the iterator over SequenceData deterministic (especially when filling holes in entity sequences)

## Version 4.9.1

### Features

- *PR #200*: MOAB now supports auto-download and configuration/installation of several dependencies (HDF5, NetCDF and Metis/ParMetis) out of the box; Use configuration options:
    - ```--download-hdf5 --download-netcdf --download-metis``` etc
    - We plan to add more dependency support in future releases.
- *PR #187*: Removal `qvdual` tool from MOAB
- *PR #218*: Remove `txt` Reader/Writer support (this was an empty shell class)
- A configuration suggestion script (suggest\_configuration.sh) that is part of MOAB repository.
    This script simplifies options to be used for new users.

### Enhancements

- Better support for OSX, Linux and LCF systems
- *PR #191*: Improve augmenting ghost entities with set information (remove parallel and memory bottlenecks)
- *PR #194*: Generate Fortran prototypes during configuration directly (no auto-generation during build)
- *PR #199*: Enhanced support for dealing with polygons and polyhedra with added functionality for VTK and Exodus formats
- *PR #202*: AHF support for high-order elements
- *PR #208*: The debian support now includes Fortran and Metis dependencies

### Fixes

- Several key fixes found through cppcheck and Coverity (static analysis tools)
- Improvements in both autotools and CMake configuration systems
- *PR #214, PR #217*: Core::connect\_iterate works with non congruently allocated cells of same type.

## Version 4.9.0

### Features

- PR #109: Optimized AHF data structure and scalability of parallel uniform mesh refinement now available (Major)
- PR #140: Merge Lasso sources to MOAB completely (Major)
- PR #137: Debian builds are available for MOAB out-of-the-box now; Check the [link](https://launchpad.net/~nschloe/+archive/ubuntu/moab-nightly/+packages "link") for updates.
- PR #172: Updated support to load geometry representations (CGM 15.0 + Cubit 14.9/15.0)

### Configuration

- CMake
    - PR #166: Missing file in CMakeLists.txt and remove exe flag from sources
    - PR #154: Don't require METIS\_DIR, PARMETIS\_DIR to be set
    - PR #150: Improve CMake HDF5 detection
    - PR #146: Remove usage of MPI\_DIR, remove config/FindMPI.cmake
    - PR #135: CMake cleanup
    - PR #124: CMake support for verdict, Metis/ParMetis

- Autotools
    - PR #168: Remove obsolete code related to netcdf workaround in moab\_mpi.h
    - PR #163: Prepare to use "make run" target in examples folder
    - PR #149: Replace the deprecated finite routine (C99 standard)
    - PR #148: Modify obb tests to not be run with regular suite
    - PR #144: Fortran libs for imeshp test
    - PR #129: Use correct linker flags for Fortran examples in debug mode
    - PR #121: Fixes for examples and upgrade makefile with options to run targets
    - PR #114: Fix compiler warnings related to HAVE\_UNORDERED\_MAP
    - PR #122: Enable MOAB\_HAVE\_VSNPRINTF for ErrorOutput.cpp

### Fixes

- General
    - PR #165: GenLargeExample tetra option fix
    - PR #157: Fix metis partitioner for serial builds
    - PR #147: Rename \_\_SDIR\_\_ macro
    - PR #132: CartVect duplicate removal
    - PR #142: simplify CpuTimer
    - PR #138: Updates to TestUtil and testMatrix3 to make sure eigenvectors (normalized) are checked correctly.
    - PR #136: Add checks for dense double tags
    - PR #130: Fix for the way Eigenvectors are returned in the Eigenvalue decomposition function.
    - PR #133: Face-to-volume Sense Corrections (MOAB)
    - PR #120: Fixed the test\_adj.
    - PR #126: Remove redundant definition of CpuTimer in WriteHDF5Parallel.cpp
    - PR #123: Fix element\_test issue due to missing mesh files
    - PR #125: Bug in Coupler
    - PR #119: Fix dual\_test issue due to missing mesh file
    - PR #169: Fix for quadratic pyramid adjacency creation

- Parallel fixes
    - PR #160: Fix hdf5 parallel reader if number of partitions is too small
    - PR #159: Run more parallel tests when built without HDF5
    - PR #158: Add option to skip augmenting the default sets with ghosts
    - PR #145: Factor in sequence manager for ghosting
    - PR #141: Populate sets with ghost entities
    - PR #131: Merge\_all method added to merge mesh class
    - I/O
        - PR #151: HDF5 parallel reader bug
        - PR #128: HDF5 file with more than 2^32 entities
    - Memory leaks
        - PR #156: Various fixes for memory leaks in serial and parallel examples
        - PR #153: Valgrind fix for parallel\_hdf5\_test
        - PR #127: Gen large mesh memory reduction
        - PR #116: Major changes for large mesh example

## Version 4.8.2

### Fixes

- PR #106: Ghosting polyhedra elements
- PR #107: AHF polyhedra fix (triggered when configured with --enable-ahf)
- PR #110: To fix attribute warning when loading an h5m file
- PR #108: Addresses issues raised by faceting on models that have point curves
- PR #113: Fix a ParallelComm issue in PGI optimization build
- PR #117: Side number should return failure when side is -1

### Enhancement

- PR #89: Writing Dense tags in parallel (minor)
- PR #104: Install configuration header to decipher dependencies (major)
- PR #115: Partitioning with Metis (major)

## Version 4.8.1

### Fixes

- PR #97: Support read and write of point clouds and free nodes in vtk format
- PR #103: Add explicitly created skin entities to the passed non-null meshset. Previously, this was only available from root set
- GMSH and VTK reader fixes
- Several minor compiler warnings removed

### Enhancement

- PR #101: Support computation of surface normals for different entity types

## Version 4.8.0

### Features

- UMR: Implemented the uniform mesh refinement capabilities for different entity types for degrees (2,3,5). Tools to load a MOAB mesh and uniformly refine it based on user inputs is now available
- Coupler: Improvements in scalability of the spatial coupler and improved capabilities to perform global and subset based conservation/normalization
- Verdict: A new set of API functions to measure quality metrics for MOAB entities (based on Verdict library)
- Enhanced error handling support with clean stack traces in both serial and parallel are now provided natively (examples included)

### Build

- Considerably improved CMake-based build support for linux/OSX/Windows
- Updated autoconf based build system to address all warnings in different architectures and more robust configuration
    - Improved support for OSX Mavericks
    - Improved 32-bit handle support
    - Support for configuration on ALCF machines out of the box
- Moved tools/mbzoltan/MBZoltan *to src/ZoltanPartitioner*
- Several bug fixes and warning removals based on GNU, Clang, PGI and Intel compilers
    - PR#59: ParallelComm - update the correct partition number when creating a part
    - PR#54: WriteNCDF - minor bug when writing out side sets
    - PR#94: Exodus Writer - fixes for variable length tags
- Additional unit tests for testing several finer grained APIs
- Several updates to the User and Developer guide in documentation to detail aspects of the new features (UMR, Verdict, Error handling)
- Overall enhanced Windows support (VS2008) (contributions by Kitware)

### Lasso

- Several code restructuring to conform library better to SIGMA stack of tools
- Updates to the autoconf build system to make it more robust
- Warning fixes for Clang/PGI/GNU/Intel
- Build in optimized mode by default

## Version 4.7.0

- Implemented the Array-based Half Facet data structures in MOAB to improve adjacency querying and efficient traversals between related entities (use `--enable-ahf` during configure)
- Co-developed (MOAB and PETSc) to introduce a new native MOAB-mesh based Discretization Manager (DM) implementation in PETSc that supports parallel solution to PDE problems on unstructured grids.
- Remove memory leak for point-to-point communication (issue #7) caused due to reuse of MPIRequest structures
- Reading sets in parallel (ticket 273)
- Parallel mbconvert bug (ticket 274, material sets compromised before writing)
- Fix memory leak for tracked sets (ticket 285)
- Fix issue #8 for merging higher dimensional entities
- Fix edge adjacency and skinner to operate on Polygon and Polyhedra
- Eliminated miscellaneous compiler warnings (tested on GCC-4.6.3, clang-3.0, Intel-13.1, PGI-13.6)
- New reader/writer formats (CGNS, MPAS, GCRM); Improve support for several existing formats (Climate-NC, HOMME, CAM-Euler, CAM-FV)
- Introduce padded polygons option for lower memory fragmentation in MPAS and GCRM climate readers
- Add climate data NC writer, in serial and parallel, for all data formats supported by the NC reader with options to append, read time-series with concatenation
- Add Zoltan repartitioning for MPAS reader
- Improvements in CUBIT reader for boundary conditions data, support CUBIT 14.0 format changes, endianness for BG/Q, variable length tags for saving block attributes and headers
- Add TET4 as supported element in Exodus reader
- Enhanced Doxygen-based documentation support for API, user guide, developer guide and examples.
- Added several new examples to demonstrate the usage of the API for different applications; More to be added here.
- Added robust CMake build support for MOAB; Preliminary fixes for Windows systems.
- Enhancements and refactoring to support reading the CGM geometry models for DagMC
- Parallelization of the MOAB based conservative intersection (CSLAM) algorithm for advection/transport applications
- Refactoring of MBCoupler and refactoring to support dynamic usage of new search tree algorithms (Kd, BVH)
- Update VTKReader Paraview plugin for MOAB; Removed the old vtkMOABReader and replaced it with vtkMOABReaderNew.
- Added meshset argument for methods in ParallelComm and Skinner
- New options for NC reader
    - rename PARTITION\_METHOD=TRIVIAL\_PARTITION to PARTITION\_METHOD=TRIVIAL
    - add NO\_MIXED\_ELEMENTS option fo reading polygon data (it means padded)
    - add NO\_EDGES option for MPAS reader
    - RENAME option for renaming a variable at writing an nc file
    - TIMESTEP option for reading nc files
- Changes in Core::list\_entities to print dense tag information
- Changes in memory evaluators, move from type long to type long long
- Changes to ScdInterface, mostly to support new partitioning method (sqijk) but also to allow for periodic meshes in all 3 directions
- Install TupleList for serial builds also
- New method in MergeMesh that allows stitching a mesh based on an integer tag defined on vertices (GLOBAL\_ID)
- Introduced interface to delete entities in ParallelComm

### Lasso

- Miscellaneous warning fixes exposed through Buildbot
- Fix some stack smashing bugs
- Added [Doxygen support][http://ftp.mcs.anl.gov/pub/fathom/lasso-docs] with nightly API built from repository
- Better flag determination during configuration

## Version 4.6

- Removed deprecated functions from the Interface: (some overloaded variants of) query\_interface, release\_interface,
    tag\_get\_data, tag\_set\_data, tag\_get\_size.
- Added Interface::adjacencies\_iterate, whose behavior is much like tag\_iterate.
- Extended Interface::get\_coords, to return entity centroids (avg vertex position) for non-vertex entities instead of failing
- Added new iMeshP extension function iMeshP\_getCommunicator, converting fortran integer to MPI C communicator
- Various bug fixes and improvements of mbcoupler (parallel mesh to mesh transfer): handling of spectral elements, works in
    serial now, better tolerancing for point searches
- New conventions for handling spectral element types; see doc/metadata\_info.doc for details
- New options for mbsize, similar to those for mbconvert (allows user to specify parallel reading options, so mbsize can be
    run as a parallel reading test)
- Initial implementation of Damsel reader/writer
- Major enhancements, efficiency improvements for NC data reader, including going mostly to a 2.5D representation, addition
    of a HOMME data reader, bug fixes in representing periodic meshes in parallel
- Initial implementation of a GCRM data reader, and a better ReadTemplate class (skeleton for writing new readers)
- Add new metadata and properties handling to DAGMC
- More extensive warning reporting when using GCC compilers, and many warning eliminations
- Support for 1D in structured mesh interface
- Improved doxygen developer documentation
- Alternative vtkMOABReaderNew & paraview plugin (in tools directory)
- Additions/improvements to examples
- New performance tests that include direct tag access tests (in test/perf/perf.cpp)
- Undeprecated one of the "tag\_get\_handle" functions for returning tag handle given tag name
- Several conventional tags (GLOBAL\_ID, MATERIAL\_SET, others) now have conventional default
    values; see src/MBTagConventions.hpp for a list of default values and other conventional tag
    characteristics.

## Version 4.5

- ITAPS: added iMesh\_stepEntIter, iMesh\_stepEntArrIter, iMesh\_tagIterate, iMesh\_createTagWithOptions (see iMesh\_extensions.h)
- More partitioning options for structured mesh (see ScdInterface::compute\_partition functions, and doc/metadata\_info.doc)
- Error class exposed in API, and query\_interface now supports getting access to the moab::Core member of that class.
- Added Interface::coords\_iterate and Interface::connect\_iterate, analogous to tag\_iterate (allows direct access to
    coordinate and connectivity memory for blocks of entity handles)
- Added new iMeshP extension tag\_reduce

## Version 4.1

- Structured mesh API (see src/moab/ScdInterface.hpp)
- Parallel read of netcdf-based .nc files using pnetcdf (see doc/metadata\_info.pdf)
- Updated ParaView plugin (see tools/vtkMOABReader/README for details)
- Direct access to dense tag storage (see tag\_iterate function in src/moab/Interface.hpp)
- Add set iterators (see src/moab/SetIterator.hpp and usage in test/perf/perf.cpp)
- Fix zoltan build on case-insensitive file systems (e.g. MacOS)
- Fix netcdf build w/ required HDF5 in non-system path

## Version 4.0.1

- Compatible with iMesh 1.2 (see README.IMESH for details on compliance testing)

## Version 4.0

- Many improvements to parallel correctness and efficiency
- Use of MPIO for parallel read of HDF5-based file format
- Parallel import of file formats utilizing internal communication and/or partial read of files.
- Partial read of HDF5-based files
- Import files from: ABAQUS, IDEAS, MCNP5, NASTRAN, Sms, TetGen, Star-CCM+
- String-based file options for controlling format-specific file options (see README.IO for a
    list of options.)
- Mesh refinement tool
- Compact storage of structured mesh data
- Variable-length tag data
- Alternate, cmake-based build system
- Support for most recent ITAPS APIs
- New data coupling tool
- Python API based on ITAPS APIs
- Many performance improvements (both runtime and memory), including entity sets, dense tag data,
    bit tag data, skinning, and entity deletion.
- MOAB namespace
- Fixed bug in get\_entities\_by\_type\_and\_tag for cases with non-zero
    input set which has or doesn't have entities

## Version 3.0.0 (SVN tag 3.0.0)

- Updated QVDual to work with new versions of VTK and removed dependence on graphviz
- Move ITAPS/TSTT interface implementation into tools/iMesh and make it work with configure system
- Implement new version number system
- Incorporate DagMC library (does fast facet-based ray tracing) into tools/dagmc
- Prefer g++ to pgCC in configure scripts
- Many improvements in kd tree functionality, performance
- Move entity sets to be stored in sequences, better performance & functionality
- Improved various file writers/readers, including:
- Better performance of HDF5 file reader
- Configuration with/without HDF5, netcdf
- Vtk writer improvements
- Added functions to MBInterface to get memory usage of MOAB
- Fixes to various MBCN functions for 64-bit builds
- Small changes to #defines of some tag names in MBTagConventions.hpp

## Version 2.00 (CVS tag version\_200)

- New MBInterface method to get blocked coordinate data for vertices.
- Speed up reading of entity sets in .h5m files.
- Store sets in entity sequences
- Remove use of virtual functions from MBMeshSet
- Add API for quering total and break-down of memory use by MOAB.
- Add initial Adaptive kD-tree implementation.
- Add upper\_bound method for MBRange.
- Make parallel configuration (MPI rank and size) per-MOAB-instance values,
    and add utility methods for querying/manipulating processor ID portion
    of entity handles.
- Fix allocation of handles such that they are allocated with the
    correct processor ID for parallel
- Remove MPI calls from MOAB constructor.  Make paralle config (MPI
    rank and size) arguments to the MOAB constuctor.
- Separate type definitions from interface definitions.  MBEntityType.h
    contains the definition of MBEntityType and MBTypes.h contains the
    definitions of all other types and includes MBEntityType.h.  MBInterface
    now includes MBTypes.h rather than MBCN.hpp so some applications using
    MBCN.hpp may have to add an explicit include.
- Add methods to MBRange to check if all contained entities are of a given
    type or dimension
- Change internal storage of entity set parent/child lists so that we have
    better behavior parent/child links become stale (don't try to delete/deref
    stale pointers.)
- Add lower\_bound, upper\_bound, equal\_range methods that accept an MBEntityType
    as input to MBRange.
- Add front, back, pop\_front and pop\_back methods to MBRange
- Change internal MBRange::PairNode definition so that the
    MBRange::const\_pair\_iterator::operator->() works correctly.
- Added 'bool' value to tag\_create, defaulting to false.  When true, tag\_create
    will return MB\_SUCCESS if the tag already exists and matches the tag
    creation parameters.
- Fixed bugs saving/restoring of mesh and default values for tags containing
    MBEntityHandles to/from HDF5 files.
- Allow special case null (zero) value for MBEntityHandle tags in HDF5 files
- Added processor rank to entity handle, right below entity type and above id
  fields; width of this field is computed at initialization, and depends on the
  number of processors being used.  On serial versions, zero bits are used so
  handles are as before.
- Added option to specify requested start id and processor id when creating
  a mesh set.
- Added functionality (in MBParallelComm class) for passing mesh between processors.
- Corrected set-related functions when inputting '0' (which is taken to mean
  the interface set, i.e. the whole mesh); in this case, one can't add parent/child
  sets, but one can request them (return no sets in that case)
- Added functions to support tag semantics
- Added num\_hops argument for num\_child\_meshsets and num\_parent\_meshsets.
- Removed default value for default value in tag\_create function (this
  argument was making the choice between overloaded versions of this
  function ambiguous)
- Fixed bug in MBCN::NumSubEntities when input d=0 (i.e. vertices)
- Changed arguments to get\_connectivity to take const MBEntityHandle\* and size
  instead of std::vector, so single MBEntityHandle can be used as input
- Added version of get\_connectivity which returns results in an MBRange,
  for convenience of calling code needing to do range-based logic
- Added std::string MBInterface::get\_error\_string(const MBErrorCode) const, which
  returns a string for the error code passed in (usually just a string representation
  of the error code enum)
- Added MBRange variants of get\_parent\_meshsets, get\_child\_meshsets
- Added list\_entity function to MBInterface
- Fix bug writing global/default values for tags containing entity handles
    to HDF5 files when using 64-bit handles.
- Fix bugs in VTK I/O for structured mesh, polygons, quadratic elements, and
    bit tags.

## Version 1.01 (CVS tag version\_101)

### New Capabilities

- Added support for polygons/polyhedra; polyhedra represented by
  storing polygon handles as connectivity array, otherwise poly elements
  similar to other entities in MOAB
- Added DualTool, to compute mesh dual, and for hex meshes, dual
  surfaces/curves
- Added support for new HDF5-based native MOAB reader/writer; this is
  the only format capable of storing any data represented in MOAB
- Added writers for GMV, SLAC, Vtk (limited), and also a template for
  constructing new mesh writers WriteTemplate
- Added tools/converter tool for mesh format conversion to/from any of
  the formats supported by MOAB
- Added support for dynamically adding readers/writers and dynamically
  testing whether any in the list can read/write a given file; required
  substantial additions to MBWriteUtil
- Added MBInterface::tag\_get\_default\_value
- Added MBRange functions subtract, lowerBound, operator+=, operator-=
- Added rudimentary mesh joining capability, and a test for that to
  MBTest
- Added "categories" tag, which represent broad category types on
  entity sets; used e.g. to indicate a set represents geometric
  vertex/edge/face/region, dual surface/curve, etc.; currently only
  supported by .cub file reader

### Bug Fixes/Rearrangements

- Fixed bug getting up-adjacencies from entities having equivalent entities;
  some adjacencies were being missed.
- Fixed a bug in set\_connectivity, where old vertices were put on the
  end of a scratch array instead of the beginning; symptom showed up as
  old vertices still being adjacent to the element whose connectivity
  was being set.
- Changed error returned when tag\_delete\_data called for entity which can't be found,
  from MB\_TAG\_NOT\_FOUND to MB\_ENTITY\_NOT\_FOUND
- Fixed bug in tag\_get\_tags, where tag handles weren't passed back properly
- Improved efficiency of MOAB's TSTT mesh interface implementation in
  various ways
- Extensive changes to performance tests (in test/perf) to test MOAB
  performance compared to that of cubit and of MOAB TSTT mesh interface
- When requesting entities with a tag value equal to the (defined)
  default value for that tag, entities not having a tag are returned in
  the list
- Move conventional tag names from MBInterface.hpp into
  MBTagConventions.hpp
- renamed MBCN::SubEntityConn to MBCN::SubEntityVertexIndices, and added
  new function MBCN::SubEntityConn, which returns the actual
  connectivity of a subentity given the parent connectivity, and the
  subentity dimension and index

## Version 1.00

Initial release (woo-hoo!)
