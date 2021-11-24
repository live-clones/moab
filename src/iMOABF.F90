! iMOAB Interface for Fortran 90/2003/2008 compatible compilers

#include "moab/MOABConfig.h"

module iMOAB

use iso_c_binding
implicit none
  
! Interface to all the API routines
interface

  integer(c_int) function iMOAB_Initialize() bind(C, name='iMOAB_InitializeFortran')
    ! Directly forward the call
    use, intrinsic :: iso_c_binding, only : c_int
  end function iMOAB_Initialize

  integer(c_int) function iMOAB_Finalize() bind(C, name='iMOAB_Finalize')
    ! Directly forward the call
    use, intrinsic :: iso_c_binding, only : c_int
  end function iMOAB_Finalize

#ifdef MOAB_HAVE_MPI
  integer(c_int) function iMOAB_RegisterApplication(app_name, fcomm, compid, pid) bind(C, name='iMOAB_RegisterApplicationFortran')
#else
  integer(c_int) function iMOAB_RegisterApplication(app_name, compid, pid) bind(C, name='iMOAB_RegisterApplicationFortran')
#endif
    ! Interface blocks don't know about their context,
    ! so we need to use iso_c_binding inside all interface functions again
    use, intrinsic :: iso_c_binding, only : c_int, c_char, c_ptr
    character(kind=c_char), intent(in) :: app_name(*)
#ifdef MOAB_HAVE_MPI
    integer, intent(in) :: fcomm
    ! type(c_ptr), value, intent(in) :: fcomm
#endif
    integer(c_int), intent(in) :: compid
    integer(c_int), intent(out) :: pid
  end function iMOAB_RegisterApplication

  integer(c_int) function iMOAB_DeregisterApplication(pid) bind(C, name='iMOAB_DeregisterApplicationFortran')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
  end function iMOAB_DeregisterApplication

  integer(c_int) function iMOAB_ReadHeaderInfo(filename, num_global_vertices, num_global_elements, num_dimension, num_parts) bind(C, name='iMOAB_ReadHeaderInfo')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    character(kind=c_char), intent(in) :: filename(*)
    integer(c_int), intent(out) :: num_global_vertices
    integer(c_int), intent(out) :: num_global_elements
    integer(c_int), intent(out) :: num_dimension
    integer(c_int), intent(out) :: num_parts
  end function iMOAB_ReadHeaderInfo


  integer(c_int) function iMOAB_LoadMesh(pid, filename, read_options, num_ghost_layers) bind(C, name='iMOAB_LoadMesh')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    integer(c_int), intent(in) :: pid
    character(kind=c_char), intent(in) :: filename(*)
    character(kind=c_char), intent(in) :: read_options(*)
    integer(c_int), intent(in) :: num_ghost_layers
  end function iMOAB_LoadMesh


  integer(c_int) function iMOAB_CreateVertices(pid, coords_len, dim, coordinates) bind(C, name='iMOAB_CreateVertices')
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: coords_len
    integer(c_int), intent(in) :: dim
    real(c_double), intent(in) :: coordinates(*)
  end function iMOAB_CreateVertices

  integer(c_int) function iMOAB_CreateElements(pid, num_elem, type, num_nodes_per_element, connectivity, block_ID) bind(C, name='iMOAB_CreateElements')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: num_elem
    integer(c_int), intent(in) :: type
    integer(c_int), intent(in) :: num_nodes_per_element
    integer(c_int), intent(in) :: connectivity(*)
    integer(c_int), intent(in) :: block_ID
  end function iMOAB_CreateElements

  integer(c_int) function iMOAB_ResolveSharedEntities(pid, num_verts, marker) bind(C, name='iMOAB_ResolveSharedEntities')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: num_verts
    integer(c_int), intent(in) :: marker(*)
  end function iMOAB_ResolveSharedEntities

  integer(c_int) function iMOAB_DetermineGhostEntities(pid, ghost_dim, num_ghost_layers, bridge_dim) bind(C, name='iMOAB_DetermineGhostEntities')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: ghost_dim
    integer(c_int), intent(in) :: num_ghost_layers
    integer(c_int), intent(in) :: bridge_dim
  end function iMOAB_DetermineGhostEntities

  integer(c_int) function iMOAB_WriteMesh(pid, filename, write_options) bind(C, name='iMOAB_WriteMesh')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    integer(c_int), intent(in) :: pid
    character(kind=c_char), intent(in) :: filename(*)
    character(kind=c_char), intent(in) :: write_options(*)
  end function iMOAB_WriteMesh

  integer(c_int) function iMOAB_UpdateMeshInfo(pid) bind(C, name='iMOAB_UpdateMeshInfo')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
  end function iMOAB_UpdateMeshInfo

  integer(c_int) function iMOAB_GetMeshInfo(pid, num_visible_vertices, num_visible_elements, &
                                              num_visible_blocks, num_visible_surfaceBC, num_visible_vertexBC) bind(C, name='iMOAB_GetMeshInfo')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(out) :: num_visible_vertices(3)
    integer(c_int), intent(out) :: num_visible_elements(3)
    integer(c_int), intent(out) :: num_visible_blocks(3)
    integer(c_int), intent(out) :: num_visible_surfaceBC(3)
    integer(c_int), intent(out) :: num_visible_vertexBC(3)
  end function iMOAB_GetMeshInfo

  integer(c_int) function iMOAB_GetVertexID(pid, vertices_length, global_vertex_ID) bind(C, name='iMOAB_GetVertexID')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: vertices_length
    integer(c_int), intent(out) :: global_vertex_ID(*)
  end function iMOAB_GetVertexID

  integer(c_int) function iMOAB_GetVertexOwnership(pid, vertices_length, visible_global_rank_ID) bind(C, name='iMOAB_GetVertexOwnership')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: vertices_length
    integer(c_int), intent(out) :: visible_global_rank_ID(*)
  end function iMOAB_GetVertexOwnership


  integer(c_int) function iMOAB_GetVisibleVerticesCoordinates(pid, coords_length, coordinates) bind(C, name='iMOAB_GetVisibleVerticesCoordinates')
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: coords_length
    real(c_double), intent(out) :: coordinates(*)
  end function iMOAB_GetVisibleVerticesCoordinates


  integer(c_int) function iMOAB_GetBlockID(pid, block_length, global_block_IDs) bind(C, name='iMOAB_GetBlockID')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: block_length
    integer(c_int), intent(out) :: global_block_IDs(*)
  end function iMOAB_GetBlockID


  integer(c_int) function iMOAB_GetBlockInfo(pid, global_block_ID, vertices_per_element, &
                                              num_elements_in_block) bind(C, name='iMOAB_GetBlockInfo')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: global_block_ID
    integer(c_int), intent(out) :: vertices_per_element
    integer(c_int), intent(out) :: num_elements_in_block
  end function iMOAB_GetBlockInfo

  integer(c_int) function iMOAB_GetVisibleElementsInfo(pid, num_visible_elements, element_global_IDs, &
                                                        ranks, block_IDs) bind(C, name='iMOAB_GetVisibleElementsInfo')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(out) :: num_visible_elements
    integer(c_int), intent(out) :: element_global_IDs(*)
    integer(c_int), intent(out) :: ranks(*)
    integer(c_int), intent(out) :: block_IDs(*)
  end function iMOAB_GetVisibleElementsInfo


  integer(c_int) function iMOAB_GetBlockElementConnectivities(pid, global_block_ID, connectivity_length, &
                                                              element_connectivity) bind(C, name='iMOAB_GetBlockElementConnectivities')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: global_block_ID
    integer(c_int), intent(in) :: connectivity_length
    integer(c_int), intent(out) :: element_connectivity(*)
  end function iMOAB_GetBlockElementConnectivities

  integer(c_int) function iMOAB_GetElementConnectivity(pid, elem_index, connectivity_length, &
                                                        element_connectivity) bind(C, name='iMOAB_GetElementConnectivity')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: elem_index
    integer(c_int), intent(in) :: connectivity_length
    integer(c_int), intent(out) :: element_connectivity(*)
  end function iMOAB_GetElementConnectivity

  integer(c_int) function iMOAB_GetElementOwnership(pid, global_block_ID, num_elements_in_block, &
                                                      element_ownership) bind(C, name='iMOAB_GetElementOwnership')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: global_block_ID
    integer(c_int), intent(out) :: num_elements_in_block
    integer(c_int), intent(out) :: element_ownership(*)
  end function iMOAB_GetElementOwnership

  integer(c_int) function iMOAB_GetElementID(pid, global_block_ID, num_elements_in_block, &
                                              global_element_ID, local_element_ID) bind(C, name='iMOAB_GetElementID')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: global_block_ID
    integer(c_int), intent(in) :: num_elements_in_block
    integer(c_int), intent(in) :: global_element_ID
    integer(c_int), intent(out) :: local_element_ID
  end function iMOAB_GetElementID

  integer(c_int) function iMOAB_GetPointerToSurfaceBC(pid, surface_BC_length, local_element_ID, &
                                                      reference_surface_ID, boundary_condition_value) bind(C, name='iMOAB_GetPointerToSurfaceBC')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: surface_BC_length
    integer(c_int), intent(in) :: local_element_ID
    integer(c_int), intent(in) :: reference_surface_ID
    integer(c_int), intent(out) :: boundary_condition_value
  end function iMOAB_GetPointerToSurfaceBC

  integer(c_int) function iMOAB_GetPointerToVertexBC(pid, vertex_BC_length, local_vertex_ID, &
                                                      boundary_condition_value) bind(C, name='iMOAB_GetPointerToVertexBC')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: vertex_BC_length
    integer(c_int), intent(in) :: local_vertex_ID
    integer(c_int), intent(out) :: boundary_condition_value
  end function iMOAB_GetPointerToVertexBC

  integer(c_int) function iMOAB_DefineTagStorage(pid, tag_storage_name, tag_type, components_per_entity, &
                                                  tag_index) bind(C, name='iMOAB_DefineTagStorage')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    integer(c_int), intent(in) :: pid
    character(kind=c_char), intent(in) :: tag_storage_name(*)
    integer(c_int), intent(in) :: tag_type
    integer(c_int), intent(in) :: components_per_entity
    integer(c_int), intent(in) :: tag_index
  end function iMOAB_DefineTagStorage

  integer(c_int) function iMOAB_SetIntTagStorage(pid, tag_storage_name, num_tag_storage_length, entity_type, &
                                                  tag_storage_data) bind(C, name='iMOAB_SetIntTagStorage')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    integer(c_int), intent(in) :: pid
    character(kind=c_char), intent(in) :: tag_storage_name(*)
    integer(c_int), intent(in) :: num_tag_storage_length
    integer(c_int), intent(in) :: entity_type
    integer(c_int), intent(in) :: tag_storage_data(*)
  end function iMOAB_SetIntTagStorage

  integer(c_int) function iMOAB_GetIntTagStorage(pid, tag_storage_name, num_tag_storage_length, entity_type, tag_storage_data) bind(C, name='iMOAB_GetIntTagStorage')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    integer(c_int), intent(in) :: pid
    character(kind=c_char), intent(in) :: tag_storage_name(*)
    integer(c_int), intent(in) :: num_tag_storage_length
    integer(c_int), intent(in) :: entity_type
    integer(c_int), intent(out) :: tag_storage_data(*)
  end function iMOAB_GetIntTagStorage

  integer(c_int) function iMOAB_SetDoubleTagStorage(pid, tag_storage_name, num_tag_storage_length, entity_type, tag_storage_data) bind(C, name='iMOAB_SetDoubleTagStorage')
    use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double
    integer(c_int), intent(in) :: pid
    character(kind=c_char), intent(in) :: tag_storage_name(*)
    integer(c_int), intent(in) :: num_tag_storage_length
    integer(c_int), intent(in) :: entity_type
    real(c_double), intent(in) :: tag_storage_data(*)
  end function iMOAB_SetDoubleTagStorage

  integer(c_int) function iMOAB_GetDoubleTagStorage(pid, tag_storage_name, num_tag_storage_length, entity_type, &
                                                      tag_storage_data) bind(C, name='iMOAB_GetDoubleTagStorage')
    use, intrinsic :: iso_c_binding, only : c_int, c_char, c_double
    integer(c_int), intent(in) :: pid
    character(kind=c_char), intent(in) :: tag_storage_name(*)
    integer(c_int), intent(in) :: num_tag_storage_length
    integer(c_int), intent(in) :: entity_type
    real(c_double), intent(out) :: tag_storage_data(*)
  end function iMOAB_GetDoubleTagStorage

  integer(c_int) function iMOAB_SynchronizeTags(pid, num_tag, tag_indices, entity_type) bind(C, name='iMOAB_SynchronizeTags')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: num_tag
    integer(c_int), intent(in) :: tag_indices(*)
    integer(c_int), intent(in) :: entity_type(*)
  end function iMOAB_SynchronizeTags

  integer(c_int) function iMOAB_ReduceTagsMax(pid, tag_index, entity_type) bind(C, name='iMOAB_ReduceTagsMax')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: tag_index
    integer(c_int), intent(in) :: entity_type
  end function iMOAB_ReduceTagsMax

  integer(c_int) function iMOAB_GetNeighborElements(pid, local_index, num_adjacent_elements, adjacent_element_IDs) bind(C, name='iMOAB_GetNeighborElements')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: local_index
    integer(c_int), intent(out) :: num_adjacent_elements
    integer(c_int), intent(out) :: adjacent_element_IDs
  end function iMOAB_GetNeighborElements

  integer(c_int) function iMOAB_GetNeighborVertices(pid, local_index, num_adjacent_vertices, adjacent_vertex_IDs) bind(C, name='iMOAB_GetNeighborVertices')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: local_index
    integer(c_int), intent(out) :: num_adjacent_vertices
    integer(c_int), intent(out) :: adjacent_vertex_IDs
  end function iMOAB_GetNeighborVertices

  integer(c_int) function iMOAB_SetGlobalInfo(pid, num_global_verts, num_global_elems) bind(C, name='iMOAB_SetGlobalInfo')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: num_global_verts
    integer(c_int), intent(in) :: num_global_elems
  end function iMOAB_SetGlobalInfo

  integer(c_int) function iMOAB_GetGlobalInfo(pid, num_global_verts, num_global_elems) bind(C, name='iMOAB_GetGlobalInfo')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(out) :: num_global_verts
    integer(c_int), intent(out) :: num_global_elems
  end function iMOAB_GetGlobalInfo

  integer(c_int) function iMOAB_WriteLocalMesh(pid, prefix) bind(C, name='iMOAB_WriteLocalMesh')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    integer(c_int), intent(in) :: pid
    character(kind=c_char), intent(in) :: prefix(*)
  end function iMOAB_WriteLocalMesh

#ifdef MOAB_HAVE_MPI

  integer(c_int) function iMOAB_SendMesh(pid, joint_comm, receivingGroup, rcompid, method) bind(C, name='iMOAB_SendMesh')
    use, intrinsic :: iso_c_binding, only : c_int, c_ptr
    integer(c_int), intent(in) :: pid
    integer, intent(in) :: joint_comm ! MPI_Comm
    integer, intent(in) :: receivingGroup ! MPI_Group
    integer(c_int), intent(in) :: rcompid
    integer(c_int), intent(in) :: method
  end function iMOAB_SendMesh

  integer(c_int) function iMOAB_ReceiveMesh(pid, joint_comm, sendingGroup, scompid) bind(C, name='iMOAB_ReceiveMesh')
    use, intrinsic :: iso_c_binding, only : c_int, c_ptr
    integer(c_int), intent(in) :: pid
    integer, intent(in) :: joint_comm ! MPI_Comm
    integer, intent(in) :: sendingGroup ! MPI_Group
    integer(c_int), intent(in) :: scompid
  end function iMOAB_ReceiveMesh

  integer(c_int) function iMOAB_FreeSenderBuffers(pid, context_id) bind(C, name='iMOAB_FreeSenderBuffers')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: context_id
  end function iMOAB_FreeSenderBuffers

  integer(c_int) function iMOAB_SendElementTag(pid, tag_storage_name, joint_comm, context_id) bind(C, name='iMOAB_SendElementTag')
    use, intrinsic :: iso_c_binding, only : c_int, c_char, c_ptr
    integer(c_int), intent(in) :: pid
    character(kind=c_char), intent(in) :: tag_storage_name(*)
    integer, intent(in) :: joint_comm ! MPI_Comm
    integer(c_int), intent(in) :: context_id
  end function iMOAB_SendElementTag

  integer(c_int) function iMOAB_ReceiveElementTag(pid, tag_storage_name, joint_comm, context_id) bind(C, name='iMOAB_ReceiveElementTag')
    use, intrinsic :: iso_c_binding, only : c_int, c_char, c_ptr
    integer(c_int), intent(in) :: pid
    character(kind=c_char), intent(in) :: tag_storage_name(*)
    integer, intent(in) :: joint_comm ! MPI_Comm
    integer(c_int), intent(in) :: context_id
  end function iMOAB_ReceiveElementTag

  integer(c_int) function iMOAB_ComputeCommGraph(pid1, pid2, joint_comm, group1, group2, type1, type2, comp1, comp2) bind(C, name='iMOAB_ComputeCommGraph')
    use, intrinsic :: iso_c_binding, only : c_int, c_ptr
    integer(c_int), intent(in) :: pid1
    integer(c_int), intent(in) :: pid2
    integer, intent(in) :: joint_comm ! MPI_Comm
    integer, intent(in) :: group1     ! MPI_Group
    integer, intent(in) :: group2     ! MPI_Group
    integer(c_int), intent(in) :: type1
    integer(c_int), intent(in) :: type2
    integer(c_int), intent(in) :: comp1
    integer(c_int), intent(in) :: comp2
  end function iMOAB_ComputeCommGraph

  integer(c_int) function iMOAB_CoverageGraph(joint_comm, pid_source, pid_migration, pid_intx, source_id, migration_id, context_id) bind(C, name='iMOAB_CoverageGraph')
    use, intrinsic :: iso_c_binding, only : c_int, c_ptr
    integer, intent(in) :: joint_comm ! MPI_Comm
    integer(c_int), intent(in) :: pid_source
    integer(c_int), intent(in) :: pid_migration
    integer(c_int), intent(in) :: pid_intx
    integer(c_int), intent(in) :: source_id
    integer(c_int), intent(in) :: migration_id
    integer(c_int), intent(in) :: context_id
  end function iMOAB_CoverageGraph

  integer(c_int) function iMOAB_DumpCommGraph(pid, context_id, is_sender, prefix) bind(C, name='iMOAB_DumpCommGraph')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    integer(c_int), intent(in) :: pid
    integer(c_int), intent(in) :: context_id
    integer(c_int), intent(in) :: is_sender
    character(kind=c_char), intent(in) :: prefix(*)
  end function iMOAB_DumpCommGraph

  integer(c_int) function iMOAB_MergeVertices(pid) bind(C, name='iMOAB_MergeVertices')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid
  end function iMOAB_MergeVertices

  integer(c_int) function iMOAB_MigrateMapMesh( pid1, pid2, pid3, jointcomm, groupA, groupB, type, comp1, comp2, direction) bind(C, name='iMOAB_MigrateMapMesh')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid1
    integer(c_int), intent(in) :: pid2
    integer(c_int), intent(in) :: pid3
    integer, intent(in) :: jointcomm  ! MPI_Comm
    integer, intent(in) :: groupA     ! MPI_Group
    integer, intent(in) :: groupB     ! MPI_Group
    integer(c_int), intent(in) :: type
    integer(c_int), intent(in) :: comp1
    integer(c_int), intent(in) :: comp2
    integer(c_int), intent(in) :: direction
  end function iMOAB_MigrateMapMesh
! closing endif: MOAB_HAVE_MPI
#endif

#ifdef MOAB_HAVE_TEMPESTREMAP

  integer(c_int) function iMOAB_ComputeMeshIntersectionOnSphere(pid_source, pid_target, pid_intersection) bind(C, name='iMOAB_ComputeMeshIntersectionOnSphere')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid_source
    integer(c_int), intent(in) :: pid_target
    integer(c_int), intent(in) :: pid_intersection
  end function iMOAB_ComputeMeshIntersectionOnSphere

  integer(c_int) function iMOAB_ComputePointDoFIntersection(pid_source, pid_target, pid_intersection) bind(C, name='iMOAB_ComputePointDoFIntersection')
    use, intrinsic :: iso_c_binding, only : c_int
    integer(c_int), intent(in) :: pid_source
    integer(c_int), intent(in) :: pid_target
    integer(c_int), intent(in) :: pid_intersection
  end function iMOAB_ComputePointDoFIntersection

#ifdef MOAB_HAVE_NETCDF

  integer(c_int) function iMOAB_LoadMappingWeightsFromFile(pid_intersection, solution_weights_identifier, remap_weights_filename) &
                                                            bind(C, name='iMOAB_LoadMappingWeightsFromFile')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    integer(c_int), intent(in) :: pid_intersection
    character(kind=c_char), intent(in) :: solution_weights_identifier(*)
    character(kind=c_char), intent(in) :: remap_weights_filename(*)
  end function iMOAB_LoadMappingWeightsFromFile

  integer(c_int) function iMOAB_LoadMappingWeights(pid_intersection, pid_cpl, col_or_row, mtype, solution_weights_identifier, remap_weights_filename) &
                                                            bind(C, name='iMOAB_LoadMappingWeights')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    integer(c_int), intent(in) :: pid_intersection
    integer(c_int), intent(in) :: pid_cpl
    integer(c_int), intent(in) :: col_or_row
    integer(c_int), intent(in) :: mtype
    character(kind=c_char), intent(in) :: solution_weights_identifier(*)
    character(kind=c_char), intent(in) :: remap_weights_filename(*)
  end function iMOAB_LoadMappingWeights

  integer(c_int) function iMOAB_WriteMappingWeightsToFile(pid_intersection, solution_weights_identifier, remap_weights_filename) bind(C, name='iMOAB_WriteMappingWeightsToFile')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    integer(c_int), intent(in) :: pid_intersection
    character(kind=c_char), intent(in) :: solution_weights_identifier(*)
    character(kind=c_char), intent(in) :: remap_weights_filename(*)
  end function iMOAB_WriteMappingWeightsToFile

! closing endif: MOAB_HAVE_NETCDF
#endif

  integer(c_int) function iMOAB_ComputeScalarProjectionWeights(pid_intersection, solution_weights_identifier, disc_method_source, disc_order_source, disc_method_target, disc_order_target, &
                                                                fNoBubble, fMonotoneTypeID, fVolumetric, fNoConservation, fValidate, &
                                                                source_solution_tag_dof_name, target_solution_tag_dof_name) bind(C, name='iMOAB_ComputeScalarProjectionWeights')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    integer(c_int), intent(in) :: pid_intersection
    character(kind=c_char), intent(in) :: solution_weights_identifier(*)
    character(kind=c_char), intent(in) :: disc_method_source(*)
    integer(c_int), intent(in) :: disc_order_source
    character(kind=c_char), intent(in) :: disc_method_target(*)
    integer(c_int), intent(in) :: disc_order_target
    integer(c_int), intent(in) :: fNoBubble
    integer(c_int), intent(in) :: fMonotoneTypeID
    integer(c_int), intent(in) :: fVolumetric
    integer(c_int), intent(in) :: fNoConservation
    integer(c_int), intent(in) :: fValidate
    character(kind=c_char), intent(in) :: source_solution_tag_dof_name(*)
    character(kind=c_char), intent(in) :: target_solution_tag_dof_name(*)
  end function iMOAB_ComputeScalarProjectionWeights

  integer(c_int) function iMOAB_ApplyScalarProjectionWeights(pid_intersection, solution_weights_identifier, source_solution_tag_name, &
                                                              target_solution_tag_name ) bind(C, name='iMOAB_ApplyScalarProjectionWeights')
    use, intrinsic :: iso_c_binding, only : c_int, c_char
    integer(c_int), intent(in) :: pid_intersection
    character(kind=c_char), intent(in) :: solution_weights_identifier(*)
    character(kind=c_char), intent(in) :: source_solution_tag_name(*)
    character(kind=c_char), intent(in) :: target_solution_tag_name(*)
  end function iMOAB_ApplyScalarProjectionWeights

! closing endif: MOAB_HAVE_TEMPESTREMAP
#endif

#ifdef MOAB_HAVE_MPI

! Add some helper MPI converter functions

  ! Convert a C-based MPI-Comm handle to fortran
  integer function iMOAB_MPI_Comm_c2f(c_handle) bind(C, name="MOAB_MPI_Comm_c2f")
    use, intrinsic :: iso_c_binding, only : c_ptr
    type(c_ptr), value :: c_handle
  end function

  ! Convert a Fortran-based MPI-integer handle to MPI-Comm pointer
  type(c_ptr) function iMOAB_MPI_Comm_f2c(f_handle) bind(C, name="MOAB_MPI_Comm_f2c")
    use, intrinsic :: iso_c_binding, only : c_ptr
    integer, intent(in) :: f_handle
  end function

  ! Convert a C-based MPI-Group handle to fortran
  integer function iMOAB_MPI_Group_c2f(c_handle) bind(C, name="MOAB_MPI_Group_c2f")
    use, intrinsic :: iso_c_binding, only : c_ptr
    type(c_ptr), value :: c_handle
  end function

  ! Convert a Fortran-based MPI-integer handle to MPI-Group pointer
  type(c_ptr) function iMOAB_MPI_Group_f2c(f_handle) bind(C, name="MOAB_MPI_Group_f2c")
    use, intrinsic :: iso_c_binding, only : c_ptr
    integer, intent(in) :: f_handle
  end function

#endif

end interface

end module iMOAB
