#ifndef IMOAB_H
#define IMOAB_H
/**
 * \file iMOAB.h
 * iMOAB: a language-agnostic, lightweight interface to MOAB.
 *
 * Supports usage from C/C++, Fortran (77/90/2003), Python.
 *
 * \remark 1) All data in the interface are exposed via POD-types.
 * \remark 2) Pass everything by reference, so we do not have to use %VAL()
 * \remark 3) Arrays are allocated by the user code. No concerns about
 *            de-allocation of the data will be taken up by the interface.
 * \remark 4) Always pass the pointer to the start of array along with the
 *            total allocated size for the array.
 * \remark 5) Return the filled array requested by user along with
 *            optionally the actual length of the array that was filled.
 *            (for typical cases, should be the allocated length)
 */

/**
 * \Notes
 * 1) Fortran MPI_Comm won't work. Take an integer argument and use MPI_F2C calls to get the C-Comm object
 * 2) ReadHeaderInfo - Does it need the pid ?
 * 3) Reuse the comm object from the registration for both load and write operations.
 *    Do not take comm objects again to avoid confusion and retain consistency.
 * 4) Decipher the global comm object and the subset partioning for each application based on the comm object
 * 5) GetMeshInfo - return separately the owned and ghosted vertices/elements -- not together in visible_** but rather
 *    owned_** and ghosted_**. Make these arrays of size 2.
 * 6) Should we sort the vertices after ghosting, such that the locally owned is first, and the ghosted appended next.
 * 7) RCM only for the owned part of the mesh -- do not screw with the ghosted layers
 * 8) GetBlockID - identical to GetVertexID -- return the global numbering for block
 * 9) GetVertexID -- remember that the order of vertices returned have an implicit numbering embedded in it.
 *    DO NOT CHANGE THIS ORDERING...
 * 10) GetBlockInfo takes global Block ID; Remove blockname unless there is a separate use case for it..
 * 11) GetElementConnectivity - clarify whether we return global or local vertex numbering.
 *     Preferably local numbering else lot of deciphering for global.
 */
#include "moab/MOABConfig.h"
#include "moab/Types.hpp"

#ifdef __cplusplus
#include <cstdio>
#include <cassert>
#include <cstdlib>
#else
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#endif

#define iMOAB_AppID    int*
#define iMOAB_String   char*
#define iMOAB_GlobalID int
#define iMOAB_LocalID  int
#ifdef __cplusplus
#define ErrCode moab::ErrorCode
#else
#define ErrCode int
#endif

#ifdef _MSC_VER
#define __PRETTY_FUNCTION__ __FUNCSIG__
#else
#if !defined( __cplusplus ) || !defined( __PRETTY_FUNCTION__ )
#define __PRETTY_FUNCTION__ __func__
#endif
#endif

/**
 * @brief MOAB tag types can be: dense/sparse, and of int/double/EntityHandle types.
 * They can also be defined on both elements and vertices of the mesh.
 */
enum MOAB_TAG_TYPE
{
    DENSE_INTEGER       = 0,
    DENSE_DOUBLE        = 1,
    DENSE_ENTITYHANDLE  = 2,
    SPARSE_INTEGER      = 3,
    SPARSE_DOUBLE       = 4,
    SPARSE_ENTITYHANDLE = 5
};

/**
 * @brief Define MOAB tag ownership information i.e., the entity where the tag is primarily defined
 * This is typically used to query the tag data and to set it back.
 */
enum MOAB_TAG_OWNER_TYPE
{
    TAG_VERTEX  = 0,
    TAG_EDGE    = 1,
    TAG_FACE    = 2,
    TAG_ELEMENT = 3
};

#define CHK_MPI_ERR( ierr )                          \
    {                                                \
        if( 0 != ( ierr ) ) return moab::MB_FAILURE; \
    }

#define IMOAB_CHECKPOINTER( prmObj, position )                                                 \
    {                                                                                          \
        if( prmObj == nullptr )                                                                \
        {                                                                                      \
            printf( "InputParamError at %d: '%s' is invalid and null.\n", position, #prmObj ); \
            return moab::MB_UNHANDLED_OPTION;                                                  \
        }                                                                                      \
    }

#define IMOAB_THROW_ERROR( message, retval )                                                                       \
    {                                                                                                              \
        printf( "iMOAB Error: %s.\n Location: %s in %s:%d.\n", message, __PRETTY_FUNCTION__, __FILE__, __LINE__ ); \
        return retval;                                                                                             \
    }

#ifndef NDEBUG
#define IMOAB_ASSERT_RET( condition, message, retval )                                                         \
    {                                                                                                          \
        if( !( condition ) )                                                                                   \
        {                                                                                                      \
            printf( "iMOAB Error: %s.\n Failed condition: %s\n Location: %s in %s:%d.\n", message, #condition, \
                    __PRETTY_FUNCTION__, __FILE__, __LINE__ );                                                 \
            return retval;                                                                                     \
        }                                                                                                      \
    }

#define IMOAB_ASSERT( condition, message ) IMOAB_ASSERT_RET( condition, message, moab::MB_UNHANDLED_OPTION )

#else
#define IMOAB_ASSERT( condition, message )
#define IMOAB_ASSERT_RET( condition, message, retval )
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Initialize the iMOAB interface implementation.
 *
 * \note Will create the MOAB instance, if not created already (reference counted).
 *
 * <B>Operations:</B> Collective
 *
 * \param[in] argc (int)           Number of command line arguments
 * \param[in] argv (iMOAB_String*) Command line arguments
 * @return ErrCode
 */
ErrCode iMOAB_Initialize( int argc, iMOAB_String* argv );

/**
 * \brief Initialize the iMOAB interface implementation from Fortran driver.
 *
 * It will create the MOAB instance, if not created already (reference counted).
 *
 * <B>Operations:</B> Collective
 *
 * \return ErrCode    The error code indicating success or failure.
 */
ErrCode iMOAB_InitializeFortran();

/**
 * \brief Finalize the iMOAB interface implementation.
 *
 * It will delete the internally reference counted MOAB instance, if the reference count reaches 0.
 *
 * <B>Operations:</B> Collective
 *
 * \return ErrCode    The error code indicating success or failure.
 */
ErrCode iMOAB_Finalize();

/**
 * \brief Register application - Create a unique application ID and bootstrap interfaces for further queries.
 *
 * \note
 * Internally, a mesh set will be associated with the application ID and all subsequent queries on the MOAB
 * instance will be directed to this mesh/file set.
 *
 * <B>Operations:</B> Collective
 *
 * \param[in]  app_name (iMOAB_String) Application name (PROTEUS, NEK5000, etc)
 * \param[in]  comm (MPI_Comm*)        MPI communicator to be used for all mesh-related queries originating
 *                                     from this application
 * \param[in]  compid (int*)           The unique external component identifier
 * \param[out] pid (iMOAB_AppID)       The unique pointer to the application ID
 * \return ErrCode                     The error code indicating success or failure.
 */
ErrCode iMOAB_RegisterApplication( const iMOAB_String app_name,
#ifdef MOAB_HAVE_MPI
                                   MPI_Comm* comm,
#endif
                                   int* compid,
                                   iMOAB_AppID pid );

/**
 * \brief De-Register application: delete mesh (set) associated with the application ID.
 *
 * \note The associated communicator will be released, and all associated mesh entities and sets will be deleted from the
 * mesh data structure. Associated tag storage data will be freed too.
 *
 * <B>Operations:</B> Collective
 *
 * \param[in] pid (iMOAB_AppID) The unique pointer to the application ID
 * \return ErrCode                     The error code indicating success or failure.
 */
ErrCode iMOAB_DeregisterApplication( iMOAB_AppID pid );

/**
 * \brief Register a Fortran-based application - Create a unique application ID and bootstrap interfaces
 * for further queries.
 *
 * \note
 * Internally, the comm object, usually a 32 bit integer in Fortran, will be converted using MPI_Comm_f2c and stored as
 * MPI_Comm. A mesh set will be associated with the application ID and all subsequent queries on the MOAB instance will
 * be directed to this mesh/file set.
 *
 * <B>Operations:</B> Collective
 *
 * \param[in]  app_name (iMOAB_String) Application name ("MySolver", "E3SM", "DMMoab", "PROTEUS", "NEK5000", etc)
 * \param[in]  communicator (int*)     Fortran MPI communicator to be used for all mesh-related queries originating from this application.
 * \param[in]  compid (int*)           External component id, (A *const* unique identifier).
 * \param[out] pid (iMOAB_AppID)       The unique pointer to the application ID.
 * \return ErrCode                     The error code indicating success or failure.
 */
ErrCode iMOAB_RegisterApplicationFortran( const iMOAB_String app_name,
#ifdef MOAB_HAVE_MPI
                                          int* communicator,
#endif
                                          int* compid,
                                          iMOAB_AppID pid );

/**
 * \brief De-Register the Fortran based application: delete mesh (set) associated with the application ID.
 *
 * \note The associated communicator will be released, and all associated mesh entities and sets will be
 * deleted from the mesh data structure. Associated tag storage data will be freed too.
 *
 * <B>Operations:</B> Collective
 *
 * \param[in] pid (iMOAB_AppID)   The unique pointer to the application ID
 * \return ErrCode                The error code indicating success or failure.
 */
ErrCode iMOAB_DeregisterApplicationFortran( iMOAB_AppID pid );

/**
 * \brief It should be called on master task only, and information obtained could be broadcasted by the user.
 * It is a fast lookup in the header of the file.
 *
 * \note
 * <B>Operations:</B> Not collective
 *
 * \param[in]  filename (iMOAB_String)    The MOAB mesh file (H5M) to probe for header information.
 * \param[out] num_global_vertices (int*) The total number of vertices in the mesh file.
 * \param[out] num_global_elements (int*) The total number of elements (of highest dimension only).
 * \param[out] num_dimension (int*)       The highest dimension of elements in the mesh
 *                                        (Edge=1, Tri/Quad=2, Tet/Hex/Prism/Pyramid=3).
 * \param[out] num_parts (int*)           The total number of partitions available in the mesh file, typically
 *                                        partitioned with mbpart during pre-processing.
 * \return ErrCode                        The error code indicating success or failure.
 */
ErrCode iMOAB_ReadHeaderInfo( const iMOAB_String filename,
                              int* num_global_vertices,
                              int* num_global_elements,
                              int* num_dimension,
                              int* num_parts );

/**
 * \brief Load a MOAB mesh file in parallel and exchange ghost layers as requested.
 *
 * \note All communication is MPI-based, and read options include parallel loading information, resolving
 * shared entities. Local MOAB instance is populated with mesh cells and vertices in the corresponding
 * local partitions.
 *
 * This will also exchange ghost cells and vertices, as requested. The default bridge dimension is 0 (vertices),
 * and all additional lower dimensional sub-entities are exchanged (mesh edges and faces). The tags in the file are not
 * exchanged by default. Default tag information for GLOBAL_ID, MATERIAL_SET, NEUMANN_SET and DIRICHLET_SET is
 * exchanged. Global ID tag is exchanged for all cells and vertices. Material sets, Neumann sets and Dirichlet sets
 * are all augmented with the ghost entities.
 *
 * <B>Operations:</B> Collective
 *
 * \param[in] pid (iMOAB_AppID)            The unique pointer to the application ID.
 * \param[in] filename (iMOAB_String)      The MOAB mesh file (H5M) to load onto the internal application mesh set.
 * \param[in] read_options (iMOAB_String)  Additional options for reading the MOAB mesh file in parallel.
 * \param[in] num_ghost_layers (int*)      The total number of ghost layers to exchange during mesh loading.
 * \return ErrCode                         The error code indicating success or failure.
 */
ErrCode iMOAB_LoadMesh( iMOAB_AppID pid,
                        const iMOAB_String filename,
                        const iMOAB_String read_options,
                        int* num_ghost_layers );

/**
 * \brief Create vertices for an app; it assumes no other vertices
 *
 * \note <B>Operations:</B> Not collective
 *
 * \param[in]  pid (iMOAB_AppID)           The unique pointer to the application ID.
 * \param[in]  coords_len (int*)           Size of the \p coordinates array (nverts * dim).
 * \param[in]  dim (int*)                  Dimension of the vertex coordinates (usually 3).
 * \param[in]  coordinates (double*)       Coordinates of all vertices, interleaved.
 * \return ErrCode                         The error code indicating success or failure.
 */
ErrCode iMOAB_CreateVertices( iMOAB_AppID pid, int* coords_len, int* dim, double* coordinates );

/**
 * \brief Create elements for an app; it assumes connectivity from local vertices, in order
 *
 * \note <B>Operations:</B> Not collective
 *
 * \param[in]  pid (iMOAB_AppID)               The unique pointer to the application ID.
 * \param[in]  num_elem (int*)                 Number of elements.
 * \param[in]  type (int*)                     Type of element (moab type).
 * \param[in]  num_nodes_per_element (int*)    Number of nodes per element.
 * \param[in]  connectivity (int *)            Connectivity array, with respect to vertices (assumed contiguous).
 * \param[in]  block_ID (int *)                The local block identifier, which will now contain the elements.
 * \return ErrCode                             The error code indicating success or failure.
 */
ErrCode iMOAB_CreateElements( iMOAB_AppID pid,
                              int* num_elem,
                              int* type,
                              int* num_nodes_per_element,
                              int* connectivity,
                              int* block_ID );

/**
 * \brief Resolve shared entities using global markers on shared vertices.
 *
 * \note  Global markers can be a global node id, for example, or a global DoF number (as for HOMME)
 *
 * <B>Operations:</B> Collective .
 *
 * \param[in] pid (iMOAB_AppID)            The unique pointer to the application ID.
 * \param[in] num_verts (int*)             Number of vertices.
 * \param[in] marker (int*)                Resolving marker (global id marker).
 * \return ErrCode                         The error code indicating success or failure.
 */
ErrCode iMOAB_ResolveSharedEntities( iMOAB_AppID pid, int* num_verts, int* marker );

/**
 * \brief Create the requested number of ghost layers for the parallel mesh.
 *
 * \note The routine assumes that the shared entities were resolved successfully, and that
 * the mesh is properly distributed on separate tasks.
 *
 * <B>Operations:</B> Collective.
 *
 * \param[in] pid (iMOAB_AppID)            The unique pointer to the application ID.
 * \param[in] ghost_dim (int*)             Desired ghost dimension (2 or 3, most of the time).
 * \param[in] num_ghost_layers (int*)      Number of ghost layers requested.
 * \param[in] bridge_dim (int*)            Bridge dimension (0 for vertices, 1 for edges, 2 for faces).
 * \return ErrCode                         The error code indicating success or failure.
 */
ErrCode iMOAB_DetermineGhostEntities( iMOAB_AppID pid, int* ghost_dim, int* num_ghost_layers, int* bridge_dim );

/**
 * \brief Write a MOAB mesh along with the solution tags to a file.
 *
 * \note The interface will write one single file (H5M) and for serial files (VTK/Exodus), it will write one
 * file per task. Write options include parallel write options, if needed. Only the mesh set and solution data
 * associated to the application will be written to the file.
 *
 * <B>Operations:</B> Collective for parallel write, non collective for serial write.
 *
 * \param[in] pid (iMOAB_AppID)            The unique pointer to the application ID.
 * \param[in] filename (iMOAB_String)      The MOAB mesh file (H5M) to write all the entities contained in the
 *                                         internal application mesh set.
 * \param[in] write_options (iMOAB_String) Additional options for writing the MOAB mesh in parallel.
 * \return ErrCode                         The error code indicating success or failure.
 */
ErrCode iMOAB_WriteMesh( iMOAB_AppID pid, const iMOAB_String filename, const iMOAB_String write_options );

/**
 * \brief Write a local MOAB mesh copy.
 *
 * \note The interface will write one file per task. Only the local mesh set and solution data associated
 * to the application will be written to the file. Very convenient method used for debugging.
 *
 * <B>Operations:</B> Not Collective (independent operation).
 *
 * \param[in] pid (iMOAB_AppID)       The unique pointer to the application ID.
 * \param[in] prefix (iMOAB_String)   The MOAB file prefix that will be used with task ID in parallel.
 * \return ErrCode                    The error code indicating success or failure.
 */
ErrCode iMOAB_WriteLocalMesh( iMOAB_AppID pid, iMOAB_String prefix );

/**
 * \brief Update local mesh data structure, from file information.
 *
 * \note The method should be called after mesh modifications, for example reading a file or creating mesh in memory
 *
 * <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)       The unique pointer to the application ID.
 * \return ErrCode                     The error code indicating success or failure.
 */
ErrCode iMOAB_UpdateMeshInfo( iMOAB_AppID pid );

/**
 * \brief Retrieve all the important mesh information related to vertices, elements, vertex and surface boundary conditions.
 *
 * \note Number of visible vertices and cells include ghost entities. All arrays returned have size 3.
 * Local entities are first, then ghost entities are next. Shared vertices can be owned in MOAB sense by
 * different tasks. Ghost vertices and cells are always owned by other tasks.
 *
 * <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)            The unique pointer to the application ID.
 * \param[out] num_visible_vertices (int*)  The number of vertices in the current partition/process arranged.
 *                                          as: owned/shared only, ghosted, total_visible (array allocated by
 *                                          user, <TT>size := 3</TT>).
 * \param[out] num_visible_elements (int*)  The number of elements in current partition/process arranged as: owned only,
 *                                          ghosted/shared, total_visible (array allocated by user, <TT>size := 3</TT>).
 * \param[out] num_visible_blocks (int*)    The number of material sets in local mesh in current partition/process
 *                                          arranged as: owned only, ghosted/shared, total_visible (array allocated by
 *                                          user, <TT>size := 3</TT>).
 * \param[out] num_visible_surfaceBC (int*) The number of mesh surfaces that have a NEUMANN_SET BC defined in local mesh
 *                                          in current partition/process arranged as: owned only, ghosted/shared,
 *                                          total_visible (array allocated by user, <TT>size := 3</TT>).
 * \param[out] num_visible_vertexBC (int*)  The number of vertices that have a DIRICHLET_SET BC defined in local mesh in
 *                                          current partition/process arranged as: owned only, ghosted/shared,
 *                                          total_visible (array allocated by user, <TT>size := 3</TT>).
 * \return ErrCode                          The error code indicating success or failure.
 */
ErrCode iMOAB_GetMeshInfo( iMOAB_AppID pid,
                           int* num_visible_vertices,
                           int* num_visible_elements,
                           int* num_visible_blocks,
                           int* num_visible_surfaceBC,
                           int* num_visible_vertexBC );

/**
 * \brief Get the global vertex IDs for all locally visible (owned and shared/ghosted) vertices.
 *
 * \note The array should be allocated by the user, sized with the total number of visible vertices
 * from iMOAB_GetMeshInfo method.
 *
 * <B>Operations:</B> Not collective
 *
 * \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID.
 * \param[in]  vertices_length (int*)              The allocated size of array
 *                                                 (typical <TT>size := num_visible_vertices</TT>).
 * \param[out] global_vertex_ID (iMOAB_GlobalID*)  The global IDs for all locally visible
 *                                                 vertices (array allocated by user).
 * \return ErrCode                                 The error code indicating success or failure.
 */
ErrCode iMOAB_GetVertexID( iMOAB_AppID pid, int* vertices_length, iMOAB_GlobalID* global_vertex_ID );

/**
 * \brief Get vertex ownership information.
 *
 * For each vertex based on the local ID, return the process that owns the vertex (local, shared or ghost)
 *
 * \note Shared vertices could be owned by different tasks. Local and shared vertices are first, ghost
 * vertices are next in the array. Ghost vertices are always owned by a different process ID. Array allocated
 * by the user with total size of visible vertices.
 *
 * <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)             The unique pointer to the application ID.
 * \param[in]  vertices_length (int*)        The allocated size of array
 *                                           (typically <TT>size := num_visible_vertices</TT>).
 * \param[out] visible_global_rank_ID (int*) The processor rank owning each of the local vertices.
 * \return ErrCode                           The error code indicating success or failure.
 */
ErrCode iMOAB_GetVertexOwnership( iMOAB_AppID pid, int* vertices_length, int* visible_global_rank_ID );

/**
 * \brief Get vertex coordinates for all local (owned and ghosted) vertices.
 *
 * \note coordinates are returned in an array allocated by user, interleaved. (do need an option for blocked
 * coordinates ?) size of the array is dimension times number of visible vertices. The local ordering is implicit,
 * owned/shared vertices are first, then ghosts.
 *
 * <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)     The unique pointer to the application ID.
 * \param[in]  coords_length (int*)  The size of the allocated coordinate array (array allocated by
 *                                   user, <TT>size := 3*num_visible_vertices</TT>).
 * \param[out] coordinates (double*) The pointer to user allocated memory that will be filled with
 *                                   interleaved coordinates.
 * \return ErrCode                   The error code indicating success or failure.
 */
ErrCode iMOAB_GetVisibleVerticesCoordinates( iMOAB_AppID pid, int* coords_length, double* coordinates );

/**
 * \brief Get the global block IDs for all locally visible (owned and shared/ghosted) blocks.
 *
 * \note Block IDs are corresponding to MATERIAL_SET tags for material sets. Usually the block ID is exported
 * from Cubit as a unique integer value. First blocks are local, and next blocks are fully ghosted. First blocks
 * have at least one owned cell/element, ghost blocks have only ghost cells. Internally, a block corresponds to a
 * mesh set with a MATERIAL_SET tag value equal to the block ID.
 *
 * <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)                  The unique pointer to the application ID.
 * \param[in]  block_length (int*)                The allocated size of array
 *                                                (typical <TT>size := num_visible_blocks</TT>).
 * \param[out] global_block_IDs (iMOAB_GlobalID*) The global IDs for all locally visible blocks
 *                                                (array allocated by user).
 * \return ErrCode                                The error code indicating success or failure.
 */
ErrCode iMOAB_GetBlockID( iMOAB_AppID pid, int* block_length, iMOAB_GlobalID* global_block_IDs );

/**
 * \brief Get the global block information and number of visible elements of belonging to a
 * block (MATERIAL SET).
 *
 * \note A block has to be homogeneous, it can contain elements of a single type.
 *
 * <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID.
 * \param[in]  global_block_ID (iMOAB_GlobalID)  The global block ID of the set to be queried.
 * \param[out] vertices_per_element (int*)       The number of vertices per element.
 * \param[out] num_elements_in_block (int*)      The number of elements in block.
 * \return ErrCode                               The error code indicating success or failure.
 */
ErrCode iMOAB_GetBlockInfo( iMOAB_AppID pid,
                            iMOAB_GlobalID* global_block_ID,
                            int* vertices_per_element,
                            int* num_elements_in_block );

/**
 * \brief Get the visible elements information.
 *
 * Return for all visible elements the global IDs, ranks they belong to, block ids they belong to.
 *
 * \todo : Can we also return the index for each block ID?
 *
 * \note <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)                     The unique pointer to the application ID.
 * \param[in]  num_visible_elements (int*)           The number of visible elements (returned by \c GetMeshInfo).
 * \param[out] element_global_IDs (iMOAB_GlobalID*)  Element global ids to be added to the block.
 * \param[out] ranks (int*)                          The owning ranks of elements.
 * \param[out] block_IDs (iMOAB_GlobalID*)           The block ids containing the elements.
 * \return ErrCode                                   The error code indicating success or failure.
 */
ErrCode iMOAB_GetVisibleElementsInfo( iMOAB_AppID pid,
                                      int* num_visible_elements,
                                      iMOAB_GlobalID* element_global_IDs,
                                      int* ranks,
                                      iMOAB_GlobalID* block_IDs );

/**
 * \brief Get the connectivities for elements contained within a certain block.
 *
 * \note input is the block ID. Should we change to visible block local index?
 *
 * <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID.
 * \param[in]  global_block_ID (iMOAB_GlobalID*) The global block ID of the set being queried.
 * \param[in]  connectivity_length (int*)        The allocated size of array.
 *                                               (typical <TT>size := vertices_per_element*num_elements_in_block</TT>)
 * \param[out] element_connectivity (int*)       The connectivity array to store element ordering in MOAB canonical
 *                                               numbering scheme (array allocated by user); array contains vertex
 *                                               indices in the local numbering order for vertices elements are in
 *                                               the same order as provided by GetElementOwnership and GetElementID.
 * \return ErrCode                               The error code indicating success or failure.
 */
ErrCode iMOAB_GetBlockElementConnectivities( iMOAB_AppID pid,
                                             iMOAB_GlobalID* global_block_ID,
                                             int* connectivity_length,
                                             int* element_connectivity );

/**
 * \brief Get the connectivity for one element only.
 *
 * \note This is a convenience method and in general should not be used unless necessary.
 * You should use colaesced calls for multiple elements for efficiency.
 *
 * <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID.
 * \param[in]  elem_index (iMOAB_LocalID *)      Local element index.
 * \param[in,out] connectivity_length (int *)    On input, maximum length of connectivity. On output, actual length.
 * \param[out] element_connectivity (int*)       The connectivity array to store connectivity in MOAB canonical numbering
 *                                               scheme. Array contains vertex indices in the local numbering order for
 *                                               vertices.
 * \return ErrCode                               The error code indicating success or failure.
 */
ErrCode iMOAB_GetElementConnectivity( iMOAB_AppID pid,
                                      iMOAB_LocalID* elem_index,
                                      int* connectivity_length,
                                      int* element_connectivity );

/**
 * \brief Get the element ownership within a certain block i.e., processor ID of the element owner.
 *
 * \note : Should we use local block index for input, instead of the global block ID ?
 *
 * <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)                 The unique pointer to the application ID.
 * \param[in]  global_block_ID (iMOAB_GlobalID)  The global block ID of the set being queried.
 * \param[in]  num_elements_in_block (int*)      The allocated size of ownership array, same as
 *                                               <TT>num_elements_in_block</TT> returned from GetBlockInfo().
 * \param[out] element_ownership (int*)          The ownership array to store processor ID for all elements
 *                                               (array allocated by user).
 * \return ErrCode                               The error code indicating success or failure.
 */
ErrCode iMOAB_GetElementOwnership( iMOAB_AppID pid,
                                   iMOAB_GlobalID* global_block_ID,
                                   int* num_elements_in_block,
                                   int* element_ownership );

/**
 * \brief Get the global IDs for all locally visible elements belonging to a particular block.
 *
 * \note The method will return also the local index of each element, in the local range that contains
 * all visible elements.
 *
 * <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID.
 * \param[in]  global_block_ID (iMOAB_GlobalID*)   The global block ID of the set being queried.
 * \param[in]  num_elements_in_block (int*)        The allocated size of global element ID array, same as
 *                                                 <TT>num_elements_in_block</TT> returned from GetBlockInfo().
 * \param[out] global_element_ID (iMOAB_GlobalID*) The global IDs for all locally visible elements
 *                                                 (array allocated by user).
 * \param[out] local_element_ID (iMOAB_LocalID*)   (<I><TT>Optional</TT></I>) The local IDs for all locally visible
 *                                                 elements (index in the range of all primary elements in the rank)
 * \return ErrCode                                 The error code indicating success or failure.
 */
ErrCode iMOAB_GetElementID( iMOAB_AppID pid,
                            iMOAB_GlobalID* global_block_ID,
                            int* num_elements_in_block,
                            iMOAB_GlobalID* global_element_ID,
                            iMOAB_LocalID* local_element_ID );

/**
 * \brief Get the surface boundary condition information.
 *
 * \note <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)                  The unique pointer to the application ID.
 * \param[in]  surface_BC_length (int*)           The allocated size of surface boundary condition array, same
 *                                                as <TT>num_visible_surfaceBC</TT> returned by GetMeshInfo().
 * \param[out] local_element_ID (iMOAB_LocalID*)  The local element IDs that contains the side with the surface BC.
 * \param[out] reference_surface_ID (int*)        The surface number with the BC in the canonical reference
 *                                                element (e.g., 1 to 6 for HEX, 1-4 for TET).
 * \param[out] boundary_condition_value (int*)    The boundary condition type as obtained from the mesh description
 *                                                (value of the NeumannSet defined on the element).
 * \return ErrCode                                The error code indicating success or failure.
 */
ErrCode iMOAB_GetPointerToSurfaceBC( iMOAB_AppID pid,
                                     int* surface_BC_length,
                                     iMOAB_LocalID* local_element_ID,
                                     int* reference_surface_ID,
                                     int* boundary_condition_value );

/**
 * \brief Get the vertex boundary condition information.
 *
 * \note <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)                   The unique pointer to the application ID.
 * \param[in]  vertex_BC_length (int)              The allocated size of vertex boundary condition array, same as
 *                                                 <TT>num_visible_vertexBC</TT> returned by GetMeshInfo().
 * \param[out] local_vertex_ID (iMOAB_LocalID*)    The local vertex ID that has Dirichlet BC defined.
 * \param[out] boundary_condition_value (int*)     The boundary condition type as obtained from the mesh description
 *                                                 (value of the Dirichlet_Set tag defined on the vertex).
 * \return ErrCode                                 The error code indicating success or failure.
 */
ErrCode iMOAB_GetPointerToVertexBC( iMOAB_AppID pid,
                                    int* vertex_BC_length,
                                    iMOAB_LocalID* local_vertex_ID,
                                    int* boundary_condition_value );

/**
 * \brief Define a MOAB Tag corresponding to the application depending on requested types.
 *
 * \note In MOAB, for most solution vectors, we only need to create a "Dense", "Double" Tag. A sparse tag can
 * be created too. If the tag is already existing in the file, it will not be created. If it is a new tag,
 * memory will be allocated when setting the values. Default values are 0 for for integer tags, 0.0 for double tags,
 * 0 for entity handle tags.
 *
 * <B>Operations:</B> Collective
 *
 * \param[in] pid (iMOAB_AppID)               The unique pointer to the application ID.
 * \param[in] tag_storage_name (iMOAB_String) The tag name to store/retrieve the data in MOAB.
 * \param[in] tag_type (int*)                 The type of MOAB tag (Dense/Sparse, Double/Int/EntityHandle), enum MOAB_TAG_TYPE.
 * \param[in] components_per_entity (int*)    The total size of vector dimension per entity for the tag (e.g., number of doubles per entity).
 * \param[out] tag_index (int*)               The tag index which can be used as identifier in synchronize methods.
 * \return ErrCode                            The error code indicating success or failure.
 */
ErrCode iMOAB_DefineTagStorage( iMOAB_AppID pid,
                                const iMOAB_String tag_storage_name,
                                int* tag_type,
                                int* components_per_entity,
                                int* tag_index );

/**
 * \brief Store the specified values in a MOAB integer Tag.
 *
 * Values are set on vertices or elements, depending on entity_type
 *
 * \note we could change the api to accept as input the tag index, as returned by iMOAB_DefineTagStorage.
 *
 * <B>Operations:</B> Collective
 *
 * \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID.
 * \param[in]  tag_storage_name (iMOAB_String)         The tag name to store/retreive the data in MOAB.
 * \param[in]  num_tag_storage_length (int*)           The size of tag storage data (e.g., num_visible_vertices*components_per_entity or
 *                                                     num_visible_elements*components_per_entity).
 * \param[in]  entity_type (int*)                      Type=0 for vertices, and Type=1 for primary elements.
 * \param[out] tag_storage_data (int*)                 The array data of type <I>int</I> to replace the internal tag memory;
 *                                                     The data is assumed to be contiguous over the local set of visible entities
 *                                                     (either vertices or elements).
 * \return ErrCode                                     The error code indicating success or failure.
 */
ErrCode iMOAB_SetIntTagStorage( iMOAB_AppID pid,
                                const iMOAB_String tag_storage_name,
                                int* num_tag_storage_length,
                                int* entity_type,
                                int* tag_storage_data );

/**
 * \brief Get the specified values in a MOAB integer Tag.
 *
 * \note <B>Operations:</B> Collective
 *
 * \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID.
 * \param[in]  tag_storage_name (iMOAB_String)         The tag name to store/retreive the data in MOAB.
 * \param[in]  num_tag_storage_length (int*)           The size of tag storage data (e.g., num_visible_vertices*components_per_entity or
 *                                                     num_visible_elements*components_per_entity).
 * \param[in]  entity_type (int*)                      Type=0 for vertices, and Type=1 for primary elements.
 * \param[out] tag_storage_data (int*)                 The array data of type <I>int</I> to be copied from the internal tag memory;
 *                                                     The data is assumed to be contiguous over the local set of visible
 *                                                     entities (either vertices or elements).
 * \return ErrCode                                     The error code indicating success or failure.
 */
ErrCode iMOAB_GetIntTagStorage( iMOAB_AppID pid,
                                const iMOAB_String tag_storage_name,
                                int* num_tag_storage_length,
                                int* entity_type,
                                int* tag_storage_data );

/**
 * \brief Store the specified values in a MOAB double Tag.
 *
 * \note <B>Operations:</B> Collective
 *
 * \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID.
 * \param[in]  tag_storage_name (iMOAB_String)         The tag names, separated, to store the data.
 * \param[in]  num_tag_storage_length (int*)           The size of total tag storage data (e.g., num_visible_vertices*components_per_entity
 *                                                     or num_visible_elements*components_per_entity*num_tags).
 * \param[in]  entity_type (int*)                      Type=0 for vertices, and Type=1 for primary elements.
 * \param[out] tag_storage_data (double*)              The array data of type <I>double</I> to replace the internal tag memory;
 *                                                     The data is assumed to be contiguous over the local set of visible
 *                                                     entities (either vertices or elements). unrolled by tags
 * \return ErrCode                                     The error code indicating success or failure.
 */
ErrCode iMOAB_SetDoubleTagStorage( iMOAB_AppID pid,
                                   const iMOAB_String tag_storage_name,
                                   int* num_tag_storage_length,
                                   int* entity_type,
                                   double* tag_storage_data );

/**
 * \brief Store the specified values in a MOAB double Tag, for given ids.
 *
 * \note <B>Operations:</B> Collective
 *
 * \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID.
 * \param[in]  tag_storage_name (iMOAB_String)         The tag names, separated, to store the data.
 * \param[in]  num_tag_storage_length (int*)           The size of total tag storage data (e.g., num_visible_vertices*components_per_entity
 *                                                     or num_visible_elements*components_per_entity*num_tags).
 * \param[in]  entity_type (int*)                      Type=0 for vertices, and Type=1 for primary elements.
 * \param[in]  tag_storage_data (double*)              The array data of type <I>double</I> to replace the internal tag memory;
 *                                                     The data is assumed to be permuted over the local set of visible
 *                                                     entities (either vertices or elements). unrolled by tags
 *                                                     in parallel, might be different order, or different entities on the task
 * \param[in]  globalIds                               global ids of the cells to be set;
 * \return ErrCode                                     The error code indicating success or failure.
 */

ErrCode iMOAB_SetDoubleTagStorageWithGid( iMOAB_AppID pid,
                                          const iMOAB_String tag_storage_names,
                                          int* num_tag_storage_length,
                                          int* ent_type,
                                          double* tag_storage_data,
                                          int* globalIds );
/**
 * \brief Retrieve the specified values in a MOAB double Tag.
 *
 * \note <B>Operations:</B> Collective
 *
 * \param[in]  pid (iMOAB_AppID)                       The unique pointer to the application ID.
 * \param[in]  tag_storage_name (iMOAB_String)         The tag names to retrieve the data in MOAB.
 * \param[in]  num_tag_storage_length (int)            The size of total tag storage data (e.g., num_visible_vertices*components_per_entity*num_tags
 *                                                     or num_visible_elements*components_per_entity*num_tags)
 * \param[in]  entity_type (int*)                      Type=0 for vertices, and Type=1 for primary elements.
 * \param[out] tag_storage_data (double*)              The array data of type <I>double</I> to replace the internal tag memory;
 *                                                     The data is assumed to be contiguous over the local set of visible
 *                                                     entities (either vertices or elements). unrolled by tags
 * \return ErrCode                                     The error code indicating success or failure.
 */
ErrCode iMOAB_GetDoubleTagStorage( iMOAB_AppID pid,
                                   const iMOAB_String tag_storage_name,
                                   int* num_tag_storage_length,
                                   int* entity_type,
                                   double* tag_storage_data );

/**
 * \brief Exchange tag values for the given tags across process boundaries.
 *
 * \note <B>Operations:</B> Collective
 *
 * \param[in]  pid (iMOAB_AppID)                The unique pointer to the application ID.
 * \param[in]  num_tags (int*)                  Number of tags to exchange.
 * \param[in]  tag_indices (int*)               Array with tag indices of interest (size  = *num_tags).
 * \param[in]  ent_type (int*)                  The type of entity for tag exchange.
 * \return ErrCode                              The error code indicating success or failure.
 */
ErrCode iMOAB_SynchronizeTags( iMOAB_AppID pid, int* num_tag, int* tag_indices, int* ent_type );

/**
 * \brief Perform global reductions with the processes in the current applications communicator.
 * Specifically this routine performs reductions on the maximum value of the tag indicated by \ref tag_index.
 *
 * \note <B>Operations:</B> Collective
 *
 * \param[in]  pid (iMOAB_AppID)                The unique pointer to the application ID.
 * \param[in]  tag_index   (int*)               The tag index of interest.
 * \param[in]  ent_type (int*)                  The type of entity for tag reduction operation
 *                                              (vertices = 0, elements = 1)
 * \return ErrCode                              The error code indicating success or failure.
 */
ErrCode iMOAB_ReduceTagsMax( iMOAB_AppID pid, int* tag_index, int* ent_type );

/**
 * \brief Retrieve the adjacencies for the element entities.
 *
 * \note <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)                      The unique pointer to the application ID.
 * \param[in]  local_index (iMOAB_LocalID*)           The local element ID for which adjacency information is needed.
 * \param[out] num_adjacent_elements (int*)           The total number of adjacent elements.
 * \param[out] adjacent_element_IDs (iMOAB_LocalID*)  The local element IDs of all adjacent elements to the current one
 *                                                    (typically, num_total_sides for internal elements or
 *                                                    num_total_sides-num_sides_on_boundary for boundary elements).
 * \return ErrCode                              The error code indicating success or failure.
 */
ErrCode iMOAB_GetNeighborElements( iMOAB_AppID pid,
                                   iMOAB_LocalID* local_index,
                                   int* num_adjacent_elements,
                                   iMOAB_LocalID* adjacent_element_IDs );

/**
 * \brief Get the adjacencies for the vertex entities.
 *
 * \todo The implementation may not be fully complete
 *
 * \note <B>Operations:</B> Not Collective
 *
 * \param[in]  pid (iMOAB_AppID)                      The unique pointer to the application ID.
 * \param[in]  local_vertex_ID (iMOAB_LocalID*)       The local vertex ID for which adjacency information is needed.
 * \param[out] num_adjacent_vertices (int*)           The total number of adjacent vertices.
 * \param[out] adjacent_vertex_IDs (iMOAB_LocalID*)   The local element IDs of all adjacent vertices to the current one
 *                                                    (typically, num_total_sides for internal elements or
 *                                                    num_total_sides-num_sides_on_boundary for boundary elements).
 * \return ErrCode                                    The error code indicating success or failure.
 */
ErrCode iMOAB_GetNeighborVertices( iMOAB_AppID pid,
                                   iMOAB_LocalID* local_vertex_ID,
                                   int* num_adjacent_vertices,
                                   iMOAB_LocalID* adjacent_vertex_IDs );

/**
 * \brief Set global information for number of vertices and number of elements;
 * It is usually available from the h5m file, or it can be computed with a MPI_Reduce.
 *
 * \param[in]  pid (iMOAB_AppID)               The unique pointer to the application ID.
 * \param[in]  num_global_verts (int*)         The number of total vertices in the mesh.
 * \param[in]  num_global_elems (MPI_Comm)     The number of total elements in the mesh.
 * \return ErrCode                             The error code indicating success or failure.
 */
ErrCode iMOAB_SetGlobalInfo( iMOAB_AppID pid, int* num_global_verts, int* num_global_elems );

/**
 * \brief Get global information about number of vertices and number of elements.
 *
 * \param[in]  pid (iMOAB_AppID)             The unique pointer to the application ID.
 * \param[in]  num_global_verts (int*)       The number of total vertices in the mesh.
 * \param[in]  num_global_elems (MPI_Comm)   The number of total elements in the mesh.
 * \return ErrCode                           The error code indicating success or failure.
 */
ErrCode iMOAB_GetGlobalInfo( iMOAB_AppID pid, int* num_global_verts, int* num_global_elems );

#ifdef MOAB_HAVE_MPI

/**
 * \brief migrate (send) a set of elements from a group of tasks (senders) to another group of tasks (receivers)
 *
 * \note <B>Operations:</B>  Collective on sender group
 *
 * \param[in]  pid (iMOAB_AppID)                  The unique pointer to the application ID of the sender mesh.
 * \param[in]  joint_communicator (MPI_Comm)      The joint communicator that overlaps both the sender and receiver groups.
 * \param[in]  receivingGroup (MPI_Group *)       Receiving group information.
 * \param[in]  rcompid  (int*)                    External id of application that receives the mesh
 * \param[in]  method (int*)                      Partitioning ethod to use when sending the mesh;
 *                                                0=trivial, 1=Zoltan Graph (PHG), 2=Zoltan Geometric (RCB).
 * \return ErrCode                                The error code indicating success or failure.
 */
ErrCode iMOAB_SendMesh( iMOAB_AppID pid,
                        MPI_Comm* joint_communicator,
                        MPI_Group* receivingGroup,
                        int* rcompid,
                        int* method );

/**
 * \brief Free up all of the buffers that were allocated when data transfer between
 * components are performed. The nonblocking send and receive call buffers will be de-allocated
 * once all the data is received by all involved tasks.
 *
 * \param[in]  pid (iMOAB_AppID)                      The unique pointer to the application ID sender mesh.
 * \param[in]  context_id  (int*)                     The context used for sending, to identify the communication graph.
 * \return ErrCode                                    The error code indicating success of failure.
 */
ErrCode iMOAB_FreeSenderBuffers( iMOAB_AppID pid, int* context_d );

/**
 * \brief Migrate (receive) a set of elements from a sender group of tasks
 *
 * \note <B>Operations:</B>  Collective on receiver group
 *
 * \param[in]  pid (iMOAB_AppID)                      The unique pointer to the application ID  mesh (receiver).
 * \param[in]  joint_communicator (MPI_Comm*)         The joint communicator that overlaps both groups.
 * \param[in]  sending_group (MPI_Group *)            The group of processes that are sending the mesh.
 * \param[in]  sender_comp_id ( int *)                The unique application identifier that is sending the mesh
 * \return ErrCode                                    The error code indicating success or failure.
 */
ErrCode iMOAB_ReceiveMesh( iMOAB_AppID pid,
                           MPI_Comm* joint_communicator,
                           MPI_Group* sending_group,
                           int* sender_comp_id );

/**
 * \brief migrate (send) a list of tags, from a sender group of tasks to a receiver group of tasks
 *
 * \note <B>Operations:</B> Collective over the sender group, nonblocking sends
 *
 * \param[in]  pid (iMOAB_AppID)                      The unique pointer to the application ID of sender mesh.
 * \param[in]  tag_storage_name (const iMOAB_String)  The name of the tags to be sent between the processes.
 *                                                    Generally multiple tag names are concatenated, separated by ";".
 * \param[in]  joint_communicator (MPI_Comm)          The joint communicator that overlaps both groups.
 * \param[in]  context_id (int*)                      The identifier for the other participating component in the
 *                                                    intersection context (typically target); -1 if this refers to
 *                                                    the original migration of meshes (internal detail).
 * \return ErrCode                                    The error code indicating success or failure.
 */
ErrCode iMOAB_SendElementTag( iMOAB_AppID pid,
                              const iMOAB_String tag_storage_name,
                              MPI_Comm* joint_communicator,
                              int* context_id );

/**
 * \brief migrate (receive) a list of tags, from a sender group of tasks to a receiver group of tasks
 *
 * \note <B>Operations:</B> Collective over the receiver group, blocking receives
 *
 * \param[in]  pid (iMOAB_AppID)                      The unique pointer to the application ID of receiver mesh.
 * \param[in]  tag_storage_name (const iMOAB_String)  The name of the tags to be sent between the processes.
 *                                                    Generally multiple tag names are concatenated, separated by ";".
 * \param[in]  joint_communicator (MPI_Comm)          The joint communicator that overlaps both groups.
 * \param[in]  context_id (int*)                      The identifier for the other participating component in the
 *                                                    intersection context (typically target); -1 if this refers to
 *                                                    the original migration of meshes (internal detail).
 * \return ErrCode                                    The error code indicating success or failure.
 */
ErrCode iMOAB_ReceiveElementTag( iMOAB_AppID pid,
                                 const iMOAB_String tag_storage_name,
                                 MPI_Comm* joint_communicator,
                                 int* context_id );

/**
 * \brief Compute a communication graph between 2 iMOAB applications, based on ID matching of key data.
 *
 * \note <B>Operations:</B> Collective
 *
 * \param[in]  pid1 (iMOAB_AppID)                     The unique pointer to the first application ID.
 * \param[in]  pid2 (iMOAB_AppID)                     The unique pointer to the second application ID.
 * \param[in]  joint_communicator (MPI_Comm*)         The MPI communicator that overlaps both groups.
 * \param[in]  group1 (MPI_Group *)                   MPI group for first component.
 * \param[in]  group2 (MPI_Group *)                   MPI group for second component.
 * \param[in]  type1 (int *)                          The type of mesh used by pid1 (spectral with GLOBAL_DOFS, etc).
 * \param[in]  type2 (int *)                          The type of mesh used by pid2 (point cloud with GLOBAL_ID, etc).
 * \param[in]  comp1 (int*)                           The unique identifier of the first component.
 * \param[in]  comp2 (int*)                           The unique identifier of the second component.
 * \return ErrCode                                    The error code indicating success or failure.
 */
ErrCode iMOAB_ComputeCommGraph( iMOAB_AppID pid1,
                                iMOAB_AppID pid2,
                                MPI_Comm* joint_communicator,
                                MPI_Group* group1,
                                MPI_Group* group2,
                                int* type1,
                                int* type2,
                                int* comp1,
                                int* comp2 );

/**
 * \brief Recompute the communication graph between component and coupler, considering intersection coverage.
 *
 * \note Original communication graph for source used an initial partition, while during intersection some of the source
 * elements were sent to multiple tasks; send back the intersection coverage information for a direct communication
 * between source cx mesh on coupler tasks and source cc mesh on interested tasks on the component.
 * The intersection tasks will send to the original source component tasks, in a nonblocking way, the ids of all the
 * cells involved in intersection with the target cells. The new ParCommGraph between cc source mesh and cx source mesh
 * will be used just for tag migration, later on; The original ParCommGraph will stay unchanged, because this source mesh
 * could be used for other intersection (atm with lnd) ? on component source tasks, we will wait for information; from each
 * intersection task, will receive cells ids involved in intersection.
 *
 * \param[in]  joint_communicator (MPI_Comm *)     The joint communicator that overlaps component PEs and coupler PEs.
 * \param[in]  pid_src (iMOAB_AppID)               The unique application identifier for the component mesh on component PEs.
 * \param[in]  pid_migr (iMOAB_AppID)              The unique application identifier for the coupler mesh on coupler PEs.
 * \param[in]  pid_intx (iMOAB_AppID)              The unique application identifier representing the intersection context on coupler PEs.
 * \param[in]  src_id (int*)                       The external id for the component mesh on component PE.
 * \param[in]  migr_id (int*)                      The external id for the migrated mesh on coupler PEs.
 * \param[in]  context_id (int*)                   The unique identifier of the other participating component in intersection (target).
 * \return ErrCode                                 The error code indicating success or failure.
 */
ErrCode iMOAB_CoverageGraph( MPI_Comm* joint_communicator,
                             iMOAB_AppID pid_src,
                             iMOAB_AppID pid_migr,
                             iMOAB_AppID pid_intx,
                             int* src_id,
                             int* migr_id,
                             int* context_id );

/**
 * \brief Dump info about communication graph.
 *
 * \note <B>Operations:</B> Collective per sender or receiver group
 *
 * \param[in] pid  (iMOAB_AppID)          The unique pointer to the application ID.
 * \param[in] context_id  (int*)          The context_id names are separated by a semi-colon (";");
 *                                        This is similar to the specifications in tag migration.
 * \param[in] is_sender (int*)            Specify whether it is called from sender or receiver side.
 * \param[in] prefix  (iMOAB_String)      The prefix for file names in order to differentiate different stages.
 * \return ErrCode                        The error code indicating success or failure.
 */
ErrCode iMOAB_DumpCommGraph( iMOAB_AppID pid, int* context_id, int* is_sender, const iMOAB_String prefix );

/**
 * \brief merge vertices in an explicit, parallel mesh; it will also reassign global IDs on vertices,
 * and resolve parallel sharing of vertices.
 *
 * \note <B>Operations:</B> Collective
 *
 * \param[in]  pid (iMOAB_AppID)    The unique pointer to the application ID.
 * \return ErrCode                        The error code indicating success or failure.
 */
ErrCode iMOAB_MergeVertices( iMOAB_AppID pid );

#endif /* #ifdef MOAB_HAVE_MPI */

#ifdef MOAB_HAVE_TEMPESTREMAP

/**
 * @brief Compute intersection of the surface meshes defined on a sphere. The resulting intersected mesh consists
 * of (convex) polygons with 1-1 associativity with both the source and destination meshes provided.
 *
 * \note  The mesh intersection between two heterogeneous resolutions will be computed and stored in the fileset
 * corresponding to the \p pid_intersection application. This intersection data can be used to compute solution
 * projection weights between these meshes.
 *
 * <B>Operations:</B> Collective on coupler tasks
 *
 * \param[in]  pid_source (iMOAB_AppID)            The unique pointer to the source application ID.
 * \param[in]  pid_target (iMOAB_AppID)            The unique pointer to the destination application ID.
 * \param[out] pid_intersection (iMOAB_AppID)      The unique pointer to the intersection application ID.
 * \return ErrCode                                    The error code indicating success or failure.
 */
ErrCode iMOAB_ComputeMeshIntersectionOnSphere( iMOAB_AppID pid_source,
                                               iMOAB_AppID pid_target,
                                               iMOAB_AppID pid_intersection );

/**
 * \brief Compute the intersection of DoFs corresponding to surface meshes defined on a sphere. The resulting intersected
 * mesh essentially contains a communication pattern or a permutation matrix that couples both the source and destination
 * DoFs.
 *
 * \note <B>Operations:</B> Collective on coupler tasks
 *
 * \param[in]  pid_source (iMOAB_AppID)            The unique pointer to the source application ID.
 * \param[in]  pid_target (iMOAB_AppID)            The unique pointer to the destination application ID.
 * \param[out] pid_intersection (iMOAB_AppID)      The unique pointer to the intersection application ID.
 * \return ErrCode                                 The error code indicating success or failure.
 */
ErrCode iMOAB_ComputePointDoFIntersection( iMOAB_AppID pid_source,
                                           iMOAB_AppID pid_target,
                                           iMOAB_AppID pid_intersection );

#ifdef MOAB_HAVE_NETCDF

#ifdef MOAB_HAVE_MPI

/**
 * \brief Compute a communication graph between 2 moab apps, based on ID matching, between a component and map that
 * was read on coupler.
 *
 * \note Component can be target or source; also migrate meshes to coupler, from the components;
 * the mesh on coupler will be either source-like == coverage mesh or target-like the map is usually read
 * with rows fully distributed on tasks, so the target mesh will be a proper partition, while source mesh
 * coverage is dictated by the active columns on respective tasks.
 *
 * <B>Operations:</B> Collective
 *
 * \param[in]  pid1 (iMOAB_AppID)         The unique pointer to the first component.
 * \param[in]  pid2 (iMOAB_AppID)         The unique pointer to the second component.
 * \param[in]  pid3 (iMOAB_AppID)         The unique pointer to the coupler instance of the mesh (component 2).
 * \param[in]  join (MPI_Comm)            The joint communicator that overlaps both groups.
 * \param[in]  group1 (MPI_Group *)       The MPI group for the first component.
 * \param[in]  group2 (MPI_Group *)       The MPI group for the second component.
 * \param[in]  type1 (int *)              The type of mesh being migrated;
 *                                        (1) spectral with GLOBAL_DOFS, (2) Point Cloud (3) FV cell.
 * \param[in]  comp1 (int*)               The universally unique identifier of first component.
 * \param[in]  comp2 (int*)               The universally unique identifier of second component.
 * \param[in]  direction (int*)           A parameter indicating direction of the mesh migration.
 *                                        i.e., whether it is from source to coupler (1), or from coupler to target (2).
 * \return ErrCode                        The error code indicating success or failure.
*/
ErrCode iMOAB_MigrateMapMesh( iMOAB_AppID pid1,
                              iMOAB_AppID pid2,
                              iMOAB_AppID pid3,
                              MPI_Comm* join,
                              MPI_Group* group1,
                              MPI_Group* group2,
                              int* type,
                              int* comp1,
                              int* comp2,
                              int* direction );

#endif /* #ifdef MOAB_HAVE_MPI */

/**
 * \brief Load the projection weights from disk to transfer a solution from a source surface mesh to a destination mesh defined on a sphere.
 * The intersection of the mesh should be computed a-priori.
 *
 * <B>Operations:</B> Collective
 *
 * \param[in] pid_intersection (iMOAB_AppID)               The unique pointer to the application ID to store the map
 * \param[in] pid_cpl (iMOAB_AppID)                        The unique pointer to coupler instance of component; (-1) for old load
 * \param[in] col_or_row (int *)                           The flag to indicate whether distribution is according to source (0) or target grid (1)
 * \param[in] type (int *)                                 type of mesh (1) spectral with GLOBAL_DOFS, (2) Point Cloud (3) FV cell
 * \param[in] solution_weights_identifier  (iMOAB_String)  The unique identifier used to store the computed projection weights locally. Typically,
 *                                                         values could be identifiers such as "scalar", "flux" or "custom".
 * \param[in] remap_weights_filename  (iMOAB_String)       The filename path to the mapping file to load in memory.
*/
ErrCode iMOAB_LoadMappingWeightsFromFile(
    iMOAB_AppID pid_intersection,
    iMOAB_AppID pid_cpl,
    int* col_or_row,
    int* type,
    const iMOAB_String solution_weights_identifier, /* "scalar", "flux", "custom" */
    const iMOAB_String remap_weights_filename );

/**
 * \brief Write the projection weights to disk in order to transfer a solution from a source surface mesh to a destination
 * mesh defined on a sphere.
 *
 * \note  <B>Operations:</B> Collective
 *
 * \param[in/out] pid_intersection (iMOAB_AppID)           The unique pointer to the application ID to store the map.
 * \param[in] solution_weights_identifier  (iMOAB_String)  The unique identifier used to store the computed projection weights locally. Typically,
 *                                                         values could be identifiers such as "scalar", "flux" or "custom".
 * \param[in] remap_weights_filename  (iMOAB_String)       The filename path to the mapping file to load in memory.
 * \return ErrCode                                         The error code indicating success or failure.
*/
ErrCode iMOAB_WriteMappingWeightsToFile(
    iMOAB_AppID pid_intersection,
    const iMOAB_String solution_weights_identifier, /* "scalar", "flux", "custom" */
    const iMOAB_String remap_weights_filename );

#endif /* #ifdef MOAB_HAVE_NETCDF */

/**
 * \brief Compute the projection weights to transfer a solution from a source surface mesh to a destination mesh defined
 * on a sphere. The intersection of the mesh should be computed a-priori.
 *
 * \note
 * The mesh intersection data-structure is used along with conservative remapping techniques in TempestRemap to compute
 * the projection weights for transferring the tag (\p soln_tag_name) from \p pid_src to \p pid_target applications.
 *
 * <B>Operations:</B> Collective
 *
 * \param[in] pid_intersection (iMOAB_AppID)                The unique pointer to the intersection application ID.
 * \param[in] solution_weights_identifier  (iMOAB_String)   The unique identifier used to store the computed projection weights locally. Typically,
 *                                                          values could be identifiers such as "scalar", "flux" or "custom".
 * \param[in] disc_method_source  (iMOAB_String)            The discretization type ("fv", "cgll", "dgll") for the solution field on the source grid.
 * \param[in] disc_order_source   (int *)                   The discretization order for the solution field on the source grid.
 * \param[in] disc_method_target  (iMOAB_String)            The discretization type ("fv", "cgll", "dgll") for the solution field on the source grid.
 * \param[in] disc_order_target   (int *)                   The discretization order for the solution field on the source grid.
 * \param[in] fNoBubble           (int *)                   The flag to indicate whether to use bubbles in the interior of SE nodes.
 *                                                          Default: true i.e., no bubbles used.
 * \param[in] fMonotoneTypeID     (int *)                   The flag to indicate whether solution monotonicity is to be preserved. 0: none, 1: mono.
 * \param[in] fVolumetric         (int *)                   The flag to indicate whether we need to compute volumetric projection weights.
 * \param[in] fNoConservation     (int *)                   The flag to indicate whether to ignore conservation of solution field during projection.
 * \param[in] fValidate           (int *)                   The flag to indicate whether to validate the consistency and conservation of solution field during projection;
 *                                                          Production runs should not have this flag enabled to minimize collective operations.
 * \param[in] source_solution_tag_dof_name   (iMOAB_String) The global DoF IDs corresponding to participating degrees-of-freedom for the source discretization.
 * \param[in] target_solution_tag_dof_name   (iMOAB_String) The global DoF IDs corresponding to participating degrees-of-freedom for the target discretization.
 * \return ErrCode                                          The error code indicating success or failure.
*/
ErrCode iMOAB_ComputeScalarProjectionWeights(
    iMOAB_AppID pid_intersection,
    const iMOAB_String solution_weights_identifier, /* "scalar", "flux", "custom" */
    const iMOAB_String disc_method_source,
    int* disc_order_source,
    const iMOAB_String disc_method_target,
    int* disc_order_target,
    int* fNoBubble,
    int* fMonotoneTypeID,
    int* fVolumetric,
    int* fInverseDistanceMap,
    int* fNoConservation,
    int* fValidate,
    const iMOAB_String source_solution_tag_dof_name,
    const iMOAB_String target_solution_tag_dof_name );

/**
 * \brief Apply the projection weights matrix operator onto the source tag in order to compute the solution (tag)
 * represented on the target grid. This operation can be understood as the application of a matrix vector product
 * (Y=P*X).
 *
 * \note <B>Operations:</B> Collective
 *
 * \param[in] pid_intersection (iMOAB_AppID)                The unique pointer to the intersection application ID.
 * \param[in] solution_weights_identifier  (iMOAB_String)   The unique identifier used to store the computed projection weights locally. Typically,
 *                                                          values could be identifiers such as "scalar", "flux" or "custom".
 * \param[in] source_solution_tag_name   (iMOAB_String)     list of tag names corresponding to participating degrees-of-freedom for the source discretization;
 *                                                          names are separated by ";", the same way as for tag migration.
 * \param[in] target_solution_tag_name   (iMOAB_String)     list of tag names corresponding to participating degrees-of-freedom for the target discretization;
 *                                                          names are separated by ";", the same way as for tag migration.
 * \return ErrCode                                          The error code indicating success or failure.
*/
ErrCode iMOAB_ApplyScalarProjectionWeights(
    iMOAB_AppID pid_intersection,
    const iMOAB_String solution_weights_identifier, /* "scalar", "flux", "custom" */
    const iMOAB_String source_solution_tag_name,
    const iMOAB_String target_solution_tag_name );

#endif /* #ifdef MOAB_HAVE_TEMPESTREMAP */

#ifdef MOAB_HAVE_MPI

/**
 * Helper functions for MPI type conversions between Fortran and C++
 */

/**
 * \brief MOAB MPI helper function to convert from a C-based MPI_Comm object
 * to a Fortran compatible MPI_Fint (integer) value.
 *
 * \param comm       Input MPI_Comm C-object pointer.
 * \return MPI_Fint  Output Fortran integer representing the comm object.
 */
inline MPI_Fint MOAB_MPI_Comm_c2f( MPI_Comm* comm )
{
    return MPI_Comm_c2f( *comm );
}

/**
 * \brief MOAB MPI helper function to convert from a Fortran-based MPI_Fint (integer)
 * value to a C/C++ compatible MPI_Comm object.
 *
 * \param fgroup (MPI_Fint)  Input Fortran integer representing the MPI_Comm object.
 * \return MPI_Comm*         Output C/C++-compatible MPI_Comm object pointer.
 */
inline MPI_Comm* MOAB_MPI_Comm_f2c( MPI_Fint fcomm )
{
    MPI_Comm* ccomm;
    ccomm  = (MPI_Comm*)( malloc( sizeof( MPI_Comm ) ) );  // memory leak ?!
    *ccomm = MPI_Comm_f2c( fcomm );
    if( *ccomm == MPI_COMM_NULL )
    {
        free( ccomm );
        IMOAB_THROW_ERROR( "The MPI_Comm conversion from Fortran to C failed.", NULL );
    }
    else
        return ccomm;
}

/**
 * \brief MOAB MPI helper functions to convert from a C-based MPI_Group object
 * to a Fortran compatible MPI_Fint (integer) value.
 *
 * \param group      Input MPI_Comm C-object pointer.
 * \return MPI_Fint  Output Fortran integer representing the comm object.
 */
inline MPI_Fint MOAB_MPI_Group_c2f( MPI_Group* group )
{
    return MPI_Group_c2f( *group );
}

/**
 * \brief MOAB MPI helper function to convert from a Fortran-based MPI_Fint (integer)
 * value to a C/C++ compatible MPI_Group object.
 *
 * \param fgroup (MPI_Fint)   Input Fortran integer representing the MPIGroup object.
 * \return MPI_Group*         Output C/C++-compatible MPI_Group object pointer.
 */
inline MPI_Group* MOAB_MPI_Group_f2c( MPI_Fint fgroup )
{
    MPI_Group* cgroup;
    cgroup  = (MPI_Group*)( malloc( sizeof( MPI_Group ) ) );  // memory leak ?!
    *cgroup = MPI_Group_f2c( fgroup );
    if( *cgroup == MPI_GROUP_NULL )
    {
        free( cgroup );
        IMOAB_THROW_ERROR( "The MPI_Group conversion from Fortran to C failed.", NULL );
    }
    else
        return cgroup;
}

#endif /* #ifdef MOAB_HAVE_MPI */

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef IMOAB_H */
