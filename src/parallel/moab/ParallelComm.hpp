/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

/**
 * \class ParallelComm
 * \brief Parallel communications in MOAB
 * \author Tim Tautges
 *
 *  This class implements methods to communicate mesh between processors
 *
 */

#ifndef MOAB_PARALLEL_COMM_HPP
#define MOAB_PARALLEL_COMM_HPP

#include "moab/Forward.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/ProcConfig.hpp"
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <stdlib.h>
#include "math.h"
#include "moab_mpi.h"

extern "C" {
  struct tuple_list;
}

static inline 
void ASSERT_ALIGNED(unsigned char *buff) 
{
#ifdef MEM_ALIGN
  assert(!((size_t)buff%MEM_ALIGN));
#endif
}

static inline 
void ALIGN_BUFFER(unsigned char *&buff, size_t extra) 
{
#ifdef MEM_ALIGN
  extra = (size_t)buff%MEM_ALIGN;
  if (extra) extra = MEM_ALIGN - extra;
  buff += extra;
  ASSERT_ALIGNED(buff);
#endif
}

template<typename T> static inline
void ALIGN_BY_SIZE(unsigned char *&buff, size_t &extra) 
{
#ifdef MEM_ALIGN
  extra = (size_t)buff%sizeof(T);
  if (extra) extra = sizeof(T) - extra;
  buff += extra;
  assert(!((size_t)buff%sizeof(T)));
#endif
}

namespace moab {

class TagServer;
class SequenceManager;
template <typename KeyType, typename ValType, ValType NullVal> class RangeMap;
typedef RangeMap<EntityHandle, EntityHandle, 0> HandleMap;

#define MAX_SHARING_PROCS 64

class ParallelComm 
{
public:

    // ==================================
    // \section CONSTRUCTORS/DESTRUCTORS/PCOMM MANAGEMENT
    // ==================================

    //! constructor
  ParallelComm(Interface *impl,
                 MPI_Comm comm = MPI_COMM_WORLD,
                 int* pcomm_id_out = 0);

    //! constructor taking packed buffer, for testing
  ParallelComm(Interface *impl,
                 std::vector<unsigned char> &tmp_buff,
                 MPI_Comm comm = MPI_COMM_WORLD,
                 int* pcomm_id_out = 0);

    //! Get ID used to reference this PCOMM instance
  int get_id() const { return pcommID; }

    //! get the indexed pcomm object from the interface
  static ParallelComm *get_pcomm(Interface *impl, const int index);
  
    //! Get ParallelComm instance associated with partition handle
    //! Will create ParallelComm instance if a) one does not already
    //! exist and b) a valid value for MPI_Comm is passed.
  static ParallelComm *get_pcomm( Interface* impl, 
                                    EntityHandle partitioning,
                                    const MPI_Comm* comm = 0 );

  static ErrorCode get_all_pcomm( Interface* impl,
                                    std::vector<ParallelComm*>& list );

    //! destructor
  ~ParallelComm();
  
  static unsigned char PROC_SHARED, PROC_OWNER;
  
    // ==================================
    // \section GLOBAL IDS
    // ==================================

    //! assign a global id space, for largest-dimension or all entities (and
    //! in either case for vertices too)
  ErrorCode assign_global_ids(EntityHandle this_set,
                                const int dimension,
                                const int start_id = 1,
                                const bool largest_dim_only = true,
                                const bool parallel = true);

    //! check for global ids; based only on tag handle being there or not;
    //! if it's not there, create them for the specified dimensions
  ErrorCode check_global_ids(EntityHandle this_set,
                               const int dimension, 
                               const int start_id = 1,
                               const bool largest_dim_only = true,
                               const bool parallel = true);
  
    // ==================================
    // \section HIGH-LEVEL COMMUNICATION (send/recv/bcast ents, exchange tags)
    // ==================================

    /** \brief send entities to another processor, optionally waiting until it's done
     *
     * Send entities to another processor, with adjs, sets, and tags.  
     * If store_remote_handles is true, this call receives back handles assigned to
     * entities sent to destination processor and stores them in sharedh_tag or 
     * sharedhs_tag.
     * \param to_proc Destination processor
     * \param orig_ents Entities requested to send
     * \param adjs If true, send adjacencies for equiv entities (currently unsupported)
     * \param tags If true, send tag values for all tags assigned to entities
     * \param store_remote_handles If true, also recv message with handles on destination processor
     * \param final_ents Range containing all entities sent
     * \param wait_all If true, wait until all messages received/sent complete
     */
  ErrorCode send_entities(const int to_proc,
                            Range &orig_ents,
                            const bool adjs,
                            const bool tags,
                            const bool store_remote_handles,
                            Range &final_ents,
                            bool wait_all = true);
  
    /** \brief Receive entities from another processor, optionally waiting until it's done
     *
     * Receive entities from another processor, with adjs, sets, and tags.  
     * If store_remote_handles is true, this call sends back handles assigned to
     * the entities received.
     * \param from_proc Source processor
     * \param store_remote_handles If true, send message with new entity handles to source processor
     * \param final_ents Range containing all entities received
     * \param wait_all If true, wait until all messages received/sent complete
     */
  ErrorCode recv_entities(const int from_proc,
                            const bool store_remote_handles,
                            Range &final_ents,
                            bool wait_all = true);
  
    /** \brief Exchange ghost cells with neighboring procs
     * Neighboring processors are those sharing an interface 
     * with this processor.  All entities of dimension ghost_dim
     * within num_layers of interface, measured going through bridge_dim,
     * are exchanged.  See MeshTopoUtil::get_bridge_adjacencies for description
     * of bridge adjacencies.  If wait_all is false and store_remote_handles
     * is true, MPI_Request objects are available in the sendReqs[2*MAX_SHARING_PROCS] 
     * member array, with inactive requests marked as MPI_REQUEST_NULL.  If
     * store_remote_handles or wait_all is false, this function returns after 
     * all entities have been received and processed.
     * \param ghost_dim Dimension of ghost entities to be exchanged
     * \param bridge_dim Dimension of entities used to measure layers from interface
     * \param num_layers Number of layers of ghosts requested
     * \param store_remote_handles If true, send message with new entity handles to source processor
     * \param wait_all If true, function does not return until all send buffers
     *       are cleared.
     */
  ErrorCode exchange_ghost_cells(int ghost_dim, int bridge_dim, 
                                   int num_layers,
                                   bool store_remote_handles,
                                   bool wait_all = true);

    /** \brief Static version of exchange_ghost_cells, exchanging info through
     * buffers rather than messages
     */
  static ErrorCode exchange_ghost_cells(ParallelComm **pc,
                                          unsigned int num_procs,
                                          int ghost_dim, int bridge_dim,
                                          int num_layers,
                                          bool store_remote_handles);
  
    /** \brief Exchange tags for all shared and ghosted entities
     * This function should be called collectively over the communicator for this ParallelComm.
     * If this version is called, all ghosted/shared entities should have a value for this
     * tag (or the tag should have a default value).
     * \param tags Vector of tag handles to be exchanged
     */
  ErrorCode exchange_tags( const std::vector<Tag> &src_tags,
                           const  std::vector<Tag> &dst_tags,
                           const Range &entities);
  
    /** \brief Exchange tags for all shared and ghosted entities
     * This function should be called collectively over the communicator for this ParallelComm
     * \param tag_name Name of tag to be exchanged
     */
  ErrorCode exchange_tags( const char *tag_name,
                           const Range &entities);
  
    /** \brief Exchange tags for all shared and ghosted entities
     * This function should be called collectively over the communicator for this ParallelComm
     * \param tagh Handle of tag to be exchanged
     */
  ErrorCode exchange_tags( Tag tagh,
                           const Range &entities);
  
    /** \brief Broadcast all entities resident on from_proc to other processors
     * This function assumes remote handles are *not* being stored, since (usually)
     * every processor will know about the whole mesh.
     * \param from_proc Processor having the mesh to be broadcast
     * \param entities On return, the entities sent or received in this call
     * \param adjacencies If true, adjacencies are sent for equiv entities (currently unsupported)
     * \param tags If true, all non-default-valued tags are sent for sent entities
     */
  ErrorCode broadcast_entities(const int from_proc,
                                 Range& entities,
                                 const bool adjacencies = false,
                                 const bool tags = true );

    // ==================================
    // \section INITIALIZATION OF PARALLEL DATA (resolve_shared_ents, etc.)
    // ==================================

    /** \brief Resolve shared entities between processors
     *
     * Resolve shared entities between processors for entities in proc_ents,
     * by comparing global id tag values on vertices on skin of elements in
     * proc_ents.  Shared entities are assigned a tag that's either
     * PARALLEL_SHARED_PROC_TAG_NAME, which is 1 integer in length, or 
     * PARALLEL_SHARED_PROCS_TAG_NAME, whose length depends on the maximum
     * number of sharing processors.  Values in these tags denote the ranks
     * of sharing processors, and the list ends with the value -1.
     *
     * If shared_dim is input as -1 or not input, a value one less than the
     * maximum dimension of entities in proc_ents is used.
     *
     * \param proc_ents Entities for which to resolve shared entities
     * \param shared_dim Maximum dimension of shared entities to look for
     */
  ErrorCode resolve_shared_ents(EntityHandle this_set,
                                  Range &proc_ents, 
                                  int resolve_dim = -1,
                                  int shared_dim = -1,
                                  const Tag* id_tag = 0);
  
    /** \brief Resolve shared entities between processors
     *
     * Same as resolve_shared_ents(Range&), except works for
     * all entities in instance of dimension dim.  
     *
     * If shared_dim is input as -1 or not input, a value one less than the
     * maximum dimension of entities is used.

     * \param dim Dimension of entities in the partition
     * \param shared_dim Maximum dimension of shared entities to look for
     */
  ErrorCode resolve_shared_ents(EntityHandle this_set,
                                  int resolve_dim = 3, 
                                  int shared_dim = -1,
                                  const Tag* id_tag = 0);
    
  static ErrorCode resolve_shared_ents(ParallelComm **pc, 
                                         const unsigned int np, 
                                         const int to_dim);
  
    // ==================================
    // \section GET PARALLEL DATA (shared/owned/iface entities, etc.)
    // ==================================

    /** \brief Get entities with the given pstatus bit(s) set
     * Returns any entities whose pstatus tag value v satisfies (v & pstatus_val)
     *
     * \param dim Dimension of entities to be returned, or -1 if any
     * \param pstatus_val pstatus value of desired entities
     * \param pstatus_ents Entities returned from function
     */
  ErrorCode get_pstatus_entities(int dim,
                                   unsigned char pstatus_val,
                                   Range &pstatus_ents);
  
    /** \brief Return the rank of the entity owner
     */
  ErrorCode get_owner(EntityHandle entity, int &owner);
  
    /** \brief Return the owner processor and handle of a given entity
     */
  ErrorCode get_owner_handle(EntityHandle entity,
                               int &owner,
                               EntityHandle &handle);

    /** \brief Get the shared processors/handles for an entity
     * Get the shared processors/handles for an entity.  Arrays must
     * be large enough to receive data for all sharing procs.
     * \param entity Entity being queried
     * \param ps Pointer to sharing proc data
     * \param hs Pointer to shared proc handle data
     * \param pstat Reference to pstatus data returned from this function
     */
  ErrorCode get_sharing_data(const EntityHandle entity,
                               int *ps, 
                               EntityHandle *hs,
                               unsigned char &pstat,
                               unsigned int &num_ps);

    /** \brief Get the shared processors/handles for an entity
     * Same as other version but with int num_ps
     * \param entity Entity being queried
     * \param ps Pointer to sharing proc data
     * \param hs Pointer to shared proc handle data
     * \param pstat Reference to pstatus data returned from this function
     */
  ErrorCode get_sharing_data(const EntityHandle entity,
                               int *ps, 
                               EntityHandle *hs,
                               unsigned char &pstat,
                               int &num_ps);

    /** \brief Get the intersection or union of all sharing processors
     * Get the intersection or union of all sharing processors.  Processor set
     * is cleared as part of this function.
     * \param entities Entity list ptr
     * \param num_entities Number of entities
     * \param procs Processors returned
     * \param op Either Interface::UNION or Interface::INTERSECT
     */
  ErrorCode get_sharing_data(const EntityHandle *entities,
                               int num_entities,
                               std::set<int> &procs,
                               int op = Interface::INTERSECT);
  
    /** \brief Get the intersection or union of all sharing processors
     * Same as previous variant but with range as input
     */
  ErrorCode get_sharing_data(const Range &entities,
                               std::set<int> &procs,
                               int op = Interface::INTERSECT);
  
    /** \brief Get entities on an inter-processor interface and of specified dimension
     * If other_proc is -1, any interface entities are returned.  If dim is -1,
     * entities of all dimensions on interface are returned.
     * \param other_proc Rank of processor for which interface entities are requested
     * \param shared_ents Entities returned from function
     * \param dim Dimension of interface entities requested
     * \param iface If true, return only entities on the interface
     */
  ErrorCode get_shared_entities(int other_proc,
                                  Range &shared_ents,
                                  int dim = -1,
                                  const bool iface = false,
                                  const bool owned_filter = false);
/*  
    //! return partition sets; if tag_name is input, gets sets with
    //! that tag name, otherwise uses PARALLEL_PARTITION tag
  ErrorCode get_partition_sets(EntityHandle this_set,
                                 Range &part_sets,
                                 const char *tag_name = NULL);
*/
    //! get processors with which this processor shares an interface
  ErrorCode get_interface_procs(std::set<unsigned int> &iface_procs,
                                  const bool get_buffs = false);

    //! get processors with which this processor communicates
  ErrorCode get_comm_procs(std::set<unsigned int> &procs);
  
    // ==================================
    // \section LOW-LEVEL DATA (tags, sets on interface/partition, etc.)
    // ==================================

    //! Get proc config for this communication object
  const ProcConfig &proc_config() const {return procConfig;}
  
    //! Get proc config for this communication object
  ProcConfig &proc_config() {return procConfig;}
  
  unsigned rank() const { return proc_config().proc_rank(); }
  unsigned size() const { return proc_config().proc_size(); }
  MPI_Comm comm() const { return proc_config().proc_comm(); }
  
    //! return the tags used to indicate shared procs and handles
  ErrorCode get_shared_proc_tags(Tag &sharedp_tag,
                                   Tag &sharedps_tag,
                                   Tag &sharedh_tag,
                                   Tag &sharedhs_tag,
                                   Tag &pstatus_tag);

    //! return partition, interface set ranges
  Range &partition_sets() {return partitionSets;}
  const Range &partition_sets() const {return partitionSets;}
  Range &interface_sets() {return interfaceSets;}
  const Range &interface_sets() const {return interfaceSets;}
      
    //! return sharedp tag
  Tag sharedp_tag();
  
    //! return sharedps tag
  Tag sharedps_tag();
  
    //! return sharedh tag
  Tag sharedh_tag();
  
    //! return sharedhs tag
  Tag sharedhs_tag();
  
    //! return pstatus tag
  Tag pstatus_tag();

    //! return pcomm tag; static because might not have a pcomm before going
    //! to look for one on the interface
  static Tag pcomm_tag(Interface *impl,
                         bool create_if_missing = true);
  
    //! return partitions set tag
  Tag partition_tag();
  Tag part_tag() { return partition_tag(); }

    // ==================================
    // \section IMESHP-RELATED FUNCTIONS
    // ==================================

    //! return all the entities in parts owned locally
  ErrorCode get_part_entities(Range &ents, int dim = -1);
  
  EntityHandle get_partitioning() const { return partitioningSet; }
  ErrorCode set_partitioning( EntityHandle h );
  ErrorCode get_global_part_count( int& count_out ) const;
  ErrorCode get_part_owner( int part_id, int& owner_out ) const;
  ErrorCode get_part_id( EntityHandle part, int& id_out ) const;
  ErrorCode get_part_handle( int id, EntityHandle& handle_out ) const;
  ErrorCode create_part( EntityHandle& part_out );
  ErrorCode destroy_part( EntityHandle part ) ;
  ErrorCode collective_sync_partition();
  ErrorCode get_part_neighbor_ids( EntityHandle part, 
                                     int neighbors_out[MAX_SHARING_PROCS],
                                     int& num_neighbors_out );
  ErrorCode get_interface_sets( EntityHandle part, 
                                  Range& iface_sets_out,
                                  int* adj_part_id = 0 );
  ErrorCode get_owning_part( EntityHandle entity, 
                               int& owning_part_id_out,
                               EntityHandle* owning_handle = 0 );
  ErrorCode get_sharing_parts( EntityHandle entity,
                                 int part_ids_out[MAX_SHARING_PROCS],
                                 int& num_part_ids_out,
                                 EntityHandle remote_handles[MAX_SHARING_PROCS] = 0);
  
    // Propogate mesh modification amongst shared entities
    // from the onwing processor to any procs with copies.
  ErrorCode update_shared_mesh();

    /** Filter the entities by pstatus tag.  
     * op is one of PSTATUS_ AND, OR, NOT; an entity is output if:
     * AND: all bits set in pstatus_val are also set on entity
     * OR: any bits set in pstatus_val also set on entity
     * NOT: any bits set in pstatus_val are not set on entity
     *
     * Results returned in input list, unless result_ents is passed in non-null,
     * in which case results are returned in result_ents.
     *
     * If ents is passed in empty, filter is done on shared entities in this
     * pcomm instance, i.e. contents of sharedEnts.
     *
     *\param ents       Input entities to filter
     *\param pstatus_val pstatus value to which entities are compared
     *\param op Bitwise operation performed between pstatus values
     *\param to_proc If non-negative and PSTATUS_SHARED is set on pstatus_val,
     *               only entities shared with to_proc are returned
     *\param result_ents If non-null, results of filter are put in the 
     *       pointed-to range
     */
  ErrorCode filter_pstatus( Range &ents,
                              const unsigned char pstatus_val,
                              const unsigned char op,
                              int to_proc = -1,
                              Range *returned_ents = NULL);

    /** \brief Get entities on interfaces shared with another proc
     *
     * \param other_proc Other proc sharing the interface
     * \param dim Dimension of entities to return, -1 if all dims
     * \param iface_ents Returned entities
     */
  ErrorCode get_iface_entities(int other_proc,
                                 int dim,
                                 Range &iface_ents);
  
  Interface* get_moab() const { return mbImpl; }

  class Buffer {
  public:
    unsigned char *mem_ptr;
    unsigned char *buff_ptr;
    unsigned int alloc_size;
    
    Buffer(unsigned int sz = 0);
    Buffer(const Buffer &);
    ~Buffer();
    void reset_buffer(size_t buff_pos = 0) {reset_ptr(buff_pos); reserve(INITIAL_BUFF_SIZE);}
    void reset_ptr(size_t buff_pos = 0) {assert((!mem_ptr && !buff_pos)|| (alloc_size >= buff_pos)); buff_ptr = mem_ptr + buff_pos;}
    void reserve(unsigned int new_size);
    void set_stored_size() {*((int*)mem_ptr) = (int)(buff_ptr - mem_ptr);}
    int get_stored_size() {return *((int*)mem_ptr);}
    unsigned long get_size() 
        {return (unsigned long)buff_ptr - (unsigned long)mem_ptr;}
          
    void check_space(unsigned int addl_space, unsigned int num_aligns = 1);
  };

    //! public 'cuz we want to unit test these externally
  ErrorCode pack_buffer(Range &orig_ents, 
                          const bool adjacencies,
                          const bool tags,
                          const bool store_remote_handles,
                          const int to_proc,
                          Buffer *buff);
  
  ErrorCode unpack_buffer(Buffer *buff,
                            const bool store_remote_handles,
                            const int from_proc,
                            const int ind,
                            std::vector<std::vector<EntityHandle> > &L1hloc,
                            std::vector<std::vector<EntityHandle> > &L1hrem,
                            std::vector<std::vector<int> > &L1p,
                            std::vector<EntityHandle> &L2hloc, 
                            std::vector<EntityHandle> &L2hrem,
                            std::vector<unsigned int> &L2p,
                            Range &new_ents);
  
  ErrorCode pack_entities(Range &entities,
                          Buffer *buff,
                            const bool store_remote_handles,
                            const int to_proc,
                            const bool is_iface,
                            std::vector<std::set<unsigned int> > *entprocs = NULL,
                            Range *allsent = NULL);

    //! unpack entities in buff_ptr
  ErrorCode unpack_entities(Buffer *buff,
                              const bool store_remote_handles,
                              const int from_ind,
                              const bool is_iface,
                              std::vector<std::vector<EntityHandle> > &L1hloc,
                              std::vector<std::vector<EntityHandle> > &L1hrem,
                              std::vector<std::vector<int> > &L1p,
                              std::vector<EntityHandle> &L2hloc, 
                              std::vector<EntityHandle> &L2hrem,
                              std::vector<unsigned int> &L2p,
                              Range &new_ents);
  
    //! Call exchange_all_shared_handles, then compare the results with tag data
    //! on local shared entities.
  ErrorCode check_all_shared_handles(bool print_em = false);

  static ErrorCode check_all_shared_handles(ParallelComm **pcs,
                                              int num_pcs);
  
  struct SharedEntityData {
    EntityHandle local;
    EntityHandle remote;
    EntityID owner;
  };

  ErrorCode pack_shared_handles(
      std::vector<std::vector<SharedEntityData> > &send_data);

    // check consistency of sharedEnts against their tags and their
    // vertices' tags
  ErrorCode check_local_shared();
  
    // check contents of communicated shared entity data against tags
  ErrorCode check_my_shared_handles(
      std::vector<std::vector<SharedEntityData> > &shents,
                                      const char *prefix = NULL);
  
    //! set rank for this pcomm; USED FOR TESTING ONLY!
  void set_rank(unsigned int r);
  
    //! set rank for this pcomm; USED FOR TESTING ONLY!
  void set_size(unsigned int r);
  
    //! get (and possibly allocate) buffers for messages to/from to_proc; returns
    //! index of to_proc in buffProcs vector; if is_new is non-NULL, sets to
    //! whether new buffer was allocated
    //! PUBLIC ONLY FOR TESTING!
  int get_buffers(int to_proc, bool *is_new = NULL);

    /* \brief Unpack message with remote handles
     * PUBLIC ONLY FOR TESTING!
     */
  ErrorCode unpack_remote_handles(unsigned int from_proc,
                                  Buffer *buff,
                                    std::vector<EntityHandle> &L2hloc,
                                    std::vector<EntityHandle> &L2hrem,
                                    std::vector<unsigned int> &L2p);
  
    /* \brief Pack message with remote handles
     * PUBLIC ONLY FOR TESTING!
     */
  ErrorCode pack_remote_handles(std::vector<EntityHandle> &L1hloc,
                                  std::vector<EntityHandle> &L1hrem,
                                  std::vector<int> &procs,
                                  unsigned int to_proc,
                                  Buffer *buff);
  
  ErrorCode list_entities(const EntityHandle *ents, int num_ents);
  
  ErrorCode list_entities(const Range &ents);
  
  static const unsigned int INITIAL_BUFF_SIZE;

private:

    // common initialization code, called from various constructors
  void initialize();
  
  ErrorCode set_sharing_data(EntityHandle ent, unsigned char pstatus,
                               int old_nump, int new_nump,
                               int *ps, EntityHandle *hs);
  
  ErrorCode check_clean_iface(Range &allsent);
  
  void define_mpe();

  ErrorCode get_sent_ents(const bool is_iface,
                            const int bridge_dim, const int ghost_dim,
                            const int num_layers,
                            Range *sent_ents, Range &allsent,
                            std::vector<std::set<unsigned int> > &entprocs);
  
    /** \brief Set pstatus values on entities
     *
     * \param pstatus_ents Entities to be set
     * \param pstatus_val Pstatus value to be set
     * \param lower_dim_ents If true, lower-dimensional ents (incl. vertices) set too
     *        (and created if they don't exist)
     * \param verts_too If true, vertices also set
     * \param operation If UNION, pstatus_val is OR-d with existing value, otherwise
     *        existing value is over-written
     */
  ErrorCode set_pstatus_entities(Range &pstatus_ents,
                                   unsigned char pstatus_val,
                                   bool lower_dim_ents = false,
                                   bool verts_too = true,
                                   int operation = Interface::UNION);

    /** \brief Set pstatus values on entities (vector-based function)
     *
     * \param pstatus_ents Entities to be set
     * \param pstatus_val Pstatus value to be set
     * \param lower_dim_ents If true, lower-dimensional ents (incl. vertices) set too
     *        (and created if they don't exist)
     * \param verts_too If true, vertices also set
     * \param operation If UNION, pstatus_val is OR-d with existing value, otherwise
     *        existing value is over-written
     */
  ErrorCode set_pstatus_entities(EntityHandle *pstatus_ents,
                                   int num_ents,
                                   unsigned char pstatus_val,
                                   bool lower_dim_ents = false,
                                   bool verts_too = true,
                                   int operation = Interface::UNION);

  int num_subranges(const Range &this_range);

    //! estimate size required to pack entities
  int estimate_ents_buffer_size(Range &entities,
                                const bool store_remote_handles,
                                unsigned int &num_aligns);
  
    //! estimate size required to pack sets
  int estimate_sets_buffer_size(Range &entities,
                                const bool store_remote_handles,
                                unsigned int &num_aligns);
  
    //! send the indicated buffer, possibly sending size first
  ErrorCode send_buffer(const unsigned int to_proc,
                          Buffer *send_buff,
                          const int msg_tag,
                          MPI_Request &send_req,
                          MPI_Request &ack_recv_req,
                          int *ack_buff,
                          int &this_incoming,
                          int next_mesg_tag = -1,
                          Buffer *next_recv_buff = NULL,
                          MPI_Request *next_recv_req = NULL,
                          int *next_incoming = NULL);
  
    //! process incoming message; if longer than the initial size, post
    //! recv for next part then send ack; if ack, send second part; else
    //! indicate that we're done and buffer is ready for processing
  ErrorCode recv_buffer(int mesg_tag_expected,
                          const MPI_Status &mpi_status,
                          Buffer *recv_buff,
                          MPI_Request &recv_2nd_req,
                          MPI_Request &ack_req,
                          int &this_incoming,
                          Buffer *send_buff,
                          MPI_Request &send_req,
                          MPI_Request &sent_ack_req,
                          bool &done,
                          Buffer *next_buff = NULL,
                          int next_tag = -1,
                          MPI_Request *next_req = NULL,
                          int *next_incoming = NULL);
  
    //! pack a range of entities with equal # verts per entity, along with
    //! the range on the sending proc
  ErrorCode pack_entity_seq(const int nodes_per_entity,
                              const bool store_remote_handles,
                              const int to_proc,
                              Range &these_ents,
                              Range &entities,
                              Buffer *buff);
  
  ErrorCode print_buffer(unsigned char *buff_ptr, int mesg_type, int from_proc,
                           bool sent);
  
    //! for all the entities in the received buffer; for each, save
    //! entities in this instance which match connectivity, or zero if none found
  ErrorCode unpack_iface_entities(Buffer *buff,
                                    const int from_proc,
                                    const int ind,
                                    std::vector<EntityHandle> &recd_ents);
  
  ErrorCode pack_sets(Range &entities,
                        Buffer *buff,
                        const bool store_handles,
                        const int to_proc);
  
  ErrorCode unpack_sets(Buffer *buff,
                          Range &entities,
                          const bool store_handles,
                          const int to_proc);
  
  ErrorCode pack_adjacencies(Range &entities,
                               Range::const_iterator &start_rit,
                               Range &whole_range,
                               Buffer *buff,
                               int &count,
                               const bool just_count,
                               const bool store_handles,
                               const int to_proc);

  ErrorCode unpack_adjacencies(Buffer *buff,
                                 Range &entities,
                                 const bool store_handles,
                                 const int from_proc);
  

    /* \brief Unpack message with remote handles (const pointer to buffer)
     */
  ErrorCode unpack_remote_handles(unsigned int from_proc,
                                  const Buffer *buff,
                                  std::vector<EntityHandle> &L2hloc,
                                  std::vector<EntityHandle> &L2hrem,
                                  std::vector<unsigned int> &L2p);
  
    //! given connectivity and type, find an existing entity, if there is one
  ErrorCode find_existing_entity(const bool is_iface,
                                   const int owner_p,
                                   const EntityHandle owner_h,
                                   const int num_ents,
                                   const EntityHandle *connect,
                                   const int num_connect,
                                   const EntityType this_type,
                                   std::vector<EntityHandle> &L2hloc,
                                   std::vector<EntityHandle> &L2hrem,
                                   std::vector<unsigned int> &L2p,
                                   EntityHandle &new_h);
  
  ErrorCode build_sharedhps_list(const EntityHandle entity,
                                   const unsigned char pstatus,
                                   const int sharedp, 
                                   const std::set<unsigned int> &entprocs,
                                   unsigned int &num_ents,
                                   int *tmp_procs,
                                   EntityHandle *tmp_handles);
  
  /**\brief Get list of tags for which to exchange data
   *
   * Get tags and entities for which to exchange tag data.  This function
   * was originally part of 'pack_tags' requested with the 
   * 'all_possible_tags' parameter.
   *
   *\param all_entities  Input.  The set of entities for which data is to 
   *                      be communicated.
   *\param all_tags      Output.  Populated with the handles of tags to be
   *                      sent.
   *\param tag_ranges    Output.  For each corresponding tag in all_tags, the
   *                      subset of 'all_entities' for which a tag value has
   *                      been set.
   */
  ErrorCode get_tag_send_list( const Range& all_entities,
                                 std::vector<Tag>& all_tags,
                                 std::vector<Range>& tag_ranges );

  /**\brief Serialize entity tag data
   *
   * This function operates in two passes.  The first phase,
   * specified by 'just_count == true' calculates the necesary
   * buffer size for the serialized data.  The second phase
   * writes the actual binary serialized representation of the
   * data to the passed buffer.
   *
   *\NOTE First two arguments are not used.  (Legacy interface?)
   *
   *\param entities      NOT USED
   *\param start_rit     NOT USED
   *\param whole_range   Should be the union of the sets of entities for 
   *                     which tag values are to be serialized.  Also
   *                     specifies ordering for indexes for tag values and
   *                     serves as the superset from which to compose entity
   *                     lists from individual tags if just_count and
   *                     all_possible_tags are both true.
   *\param buff_ptr      Buffer into which to write binary serailzed data
   *\param count         Output:  The size of the serialized data is added
   *                     to this parameter.  NOTE: Should probalby initialize
   *                     to zero before calling.
   *\param just_count    If true, just calculate the buffer size required to
   *                     hold the serialized data.  Will also append to
   *                     'all_tags' and 'tag_ranges' if all_possible_tags
   *                     == true.
   *\param store_handles The data for each tag is preceeded by a list of 
   *                     EntityHandles designating the entity each of
   *                     the subsequent tag values corresponds to.  This value
   *                     may be one of:
   *                     1) If store_handles == false:
   *                        An invalid handle composed of {MBMAXTYPE,idx}, where
   *                        idx is the position of the entity in "whole_range".
   *                     2) If store_hanldes == true and a valid remote
   *                        handle exists, the remote handle.
   *                     3) If store_hanldes == true and no valid remote 
   *                        handle is defined for the entity, the same as 1).
   *\param to_proc       If 'store_handles' is true, the processor rank for
   *                     which to store the corresponding remote entity 
   *                     handles.
   *\param all_tags      List of tags to write
   *\param tag_ranges    List of entities to serialize tag data, one
   *                            for each corresponding tag handle in 'all_tags.
   */
  ErrorCode pack_tags(Range &entities,
                        const std::vector<Tag> &src_tags,
                        const std::vector<Tag> &dst_tags,
                        const std::vector<Range> &tag_ranges,
                        Buffer *buff,
                        const bool store_handles,
                        const int to_proc);

    /**\brief Calculate buffer size required to packtag data
     *\param source_tag The tag for which data will be serialized
     *\param entites    The entities for which tag values will be serialized
     *\param count_out  Output: The required buffer size, in bytes.
     */
  ErrorCode packed_tag_size( Tag source_tag, 
                               const Range& entities, 
                               int& count_out );
  
  /**\brief Serialize tag data
   *\param source_tag    The tag for which data will be serialized
   *\param destination_tag Tag in which to store unpacked tag data.  Typically
   *                     the same as source_tag.
   *\param entites       The entities for which tag values will be serialized
   *\param whole_range   Calculate entity indices as location in this range
   *\param buff_ptr      Input/Output: As input, pointer to the start of the
   *                     buffer in which to serialize data.  As output, the
   *                     position just passed the serialized data.
   *\param count_out     Output: The required buffer size, in bytes.
   *\param store_handles The data for each tag is preceeded by a list of 
   *                     EntityHandles designating the entity each of
   *                     the subsequent tag values corresponds to.  This value
   *                     may be one of:
   *                     1) If store_handles == false:
   *                        An invalid handle composed of {MBMAXTYPE,idx}, where
   *                        idx is the position of the entity in "whole_range".
   *                     2) If store_hanldes == true and a valid remote
   *                        handle exists, the remote handle.
   *                     3) If store_hanldes == true and no valid remote 
   *                        handle is defined for the entity, the same as 1).
   *\param to_proc       If 'store_handles' is true, the processor rank for
   *                     which to store the corresponding remote entity 
   *                     handles.
   */
  ErrorCode pack_tag( Tag source_tag,
                        Tag destination_tag,
                        const Range &entities,
                        const Range &whole_range,
                        Buffer *buff,
                        const bool store_remote_handles,
                        const int to_proc );

  ErrorCode unpack_tags(Buffer *buff,
                          Range &entities,
                          const bool store_handles,
                          const int to_proc);
  
  ErrorCode tag_shared_verts(tuple_list &shared_verts,
                               Range *skin_ents,
                               std::map<std::vector<int>, Range> &proc_nranges,
                               Range &proc_verts);
  
  ErrorCode tag_shared_ents(int resolve_dim,
                              int shared_dim,
                              Range *skin_ents,
                              std::map<std::vector<int>, Range> &proc_nranges);

    // each iterate in proc_nranges contains a set of procs and the entities *possibly*
    // on the interface between those procs; this function makes sets for each,
    // and tags the set with the procs sharing it; interface sets are optionally
    // returned; NOTE: a subsequent step is used to verify entities on the interface
    // and remove them if they're not shared
  ErrorCode create_interface_sets(std::map<std::vector<int>, Range> &proc_nranges,
                                    int resolve_dim, int shared_dim);

    // do the same but working straight from sharedEnts
  ErrorCode create_interface_sets(int resolve_dim, int shared_dim);

    // after verifying shared entities, now parent/child links between sets can be established
  ErrorCode create_iface_pc_links();
  
    //! pack a range map with keys in this_range and values a contiguous series
    //! of handles starting at actual_start
  ErrorCode pack_range_map(Range &this_range, EntityHandle actual_start,
                             HandleMap &handle_map);

    //! returns true if the set is an interface shared with to_proc
  bool is_iface_proc(EntityHandle this_set, int to_proc);
  
    //! for any remote_handles set to zero, remove corresponding sent_ents from
    //! iface_sets corresponding to from_proc
  ErrorCode update_iface_sets(Range &sent_ents,
                                std::vector<EntityHandle> &remote_handles, 
                                int from_proc);
  
    //! for specified bridge/ghost dimension, to_proc, and number
    //! of layers, get the entities to be ghosted, and info on additional procs
    //! needing to communicate with to_proc
  ErrorCode get_ghosted_entities(int bridge_dim,
                                   int ghost_dim,
                                   int to_proc, 
                                   int num_layers,
                                   Range &ghosted_ents);
  
    //! add vertices adjacent to entities in this list
  ErrorCode add_verts(Range &sent_ents);
  
  //! Every processor sends shared entity handle data to every other processor
  //! that it shares entities with.  Passed back map is all received data,
  //! indexed by processor ID. This function is intended to be used for 
  //! debugging.
  ErrorCode exchange_all_shared_handles(  
      std::vector<std::vector<SharedEntityData> > &send_data, 
      std::vector<std::vector<SharedEntityData> > &result);
  
    //! replace handles in from_vec with corresponding handles on
    //! to_proc (by checking shared[p/h]_tag and shared[p/h]s_tag;
    //! if no remote handle and new_ents is non-null, substitute
    //! instead CREATE_HANDLE(MBMAXTYPE, index) where index is handle's
    //! position in new_ents
  ErrorCode get_remote_handles(const bool store_remote_handles,
                                 EntityHandle *from_vec, 
                                 EntityHandle *to_vec_tmp,
                                 int num_ents, int to_proc,
                                 const Range &new_ents);
  
    //! same as other version, except from_range and to_range should be
    //! different here
  ErrorCode get_remote_handles(const bool store_remote_handles,
                                 const Range &from_range, 
                                 Range &to_range,
                                 int to_proc,
                                 const Range &new_ents);
  
    //! same as other version, except packs range into vector
  ErrorCode get_remote_handles(const bool store_remote_handles,
                                 const Range &from_range, 
                                 EntityHandle *to_vec,
                                 int to_proc,
                                 const Range &new_ents);

  std::vector<unsigned int> &buff_procs();

    //! goes through from_vec, and for any with type MBMAXTYPE, replaces with
    //! new_ents value at index corresponding to id of entity in from_vec
  ErrorCode get_local_handles(EntityHandle *from_vec, 
                                int num_ents,
                                const Range &new_ents);

    //! same as above except puts results in range
  ErrorCode get_local_handles(const Range &remote_handles,
                                Range &local_handles,
                                const Range &new_ents);
  
    //! same as above except gets new_ents from vector
  ErrorCode get_local_handles(EntityHandle *from_vec,
                                int num_ents,
                                const std::vector<EntityHandle> &new_ents);
  
  ErrorCode update_remote_data(Range &local_range,
                                 Range &remote_range,
                                 int other_proc,
                                 const unsigned char add_pstat);
  
  ErrorCode update_remote_data(const EntityHandle new_h,
                                 const int *ps,
                                 const EntityHandle *hs,
                                 const int num_ps,
                                 const unsigned char add_pstat);
  
    /** \brief Set pstatus tag interface bit on entities in sets passed in
     */
  ErrorCode tag_iface_entities();

    //! add a pc to the iface instance tag PARALLEL_COMM
  int add_pcomm(ParallelComm *pc);
  
    //! remove a pc from the iface instance tag PARALLEL_COMM
  void remove_pcomm(ParallelComm *pc);
  
    //! check entities to make sure there are no zero-valued remote handles
    //! where they shouldn't be
  ErrorCode check_sent_ents(Range &allsent);

    //! MB interface associated with this writer
  Interface *mbImpl;

    //! Proc config object, keeps info on parallel stuff
  ProcConfig procConfig;
  
    //! Tag server, so we can get more info about tags
  TagServer *tagServer;
  
    //! Sequence manager, to get more efficient access to entities
  SequenceManager *sequenceManager;
  
    //! more data buffers, proc-specific
  std::vector<Buffer*> localOwnedBuffs, remoteOwnedBuffs;

    //! reset message buffers to their initial state
  void reset_all_buffers();

    //! delete all buffers, freeing up any memory held by them
  void delete_all_buffers();

    //! request objects, may be used if store_remote_handles is used
  std::vector<MPI_Request> sendReqs;

    //! processor rank for each buffer index
  std::vector<unsigned int> buffProcs;

    //! the partition, interface sets for this comm'n instance
  Range partitionSets, interfaceSets;

    //! all local entities shared with others, whether ghost or ghosted
  Range sharedEnts;
  
    //! tags used to save sharing procs and handles
  Tag sharedpTag, sharedpsTag, sharedhTag, sharedhsTag, pstatusTag, 
      ifaceSetsTag, partitionTag;
    
  int globalPartCount; //!< Cache of global part count
  
  EntityHandle partitioningSet; //!< entity set containing all parts

  std::ofstream myFile;
  
  int pcommID;

};

inline ParallelComm::Buffer::Buffer(const Buffer &other_buff) 
{
  alloc_size = other_buff.alloc_size;
  mem_ptr = (unsigned char *)malloc(alloc_size);
//  MALLOC(mem_ptr, alloc_size, unsigned char*);
  
  memcpy(mem_ptr, other_buff.mem_ptr, alloc_size);
  buff_ptr = mem_ptr + (other_buff.buff_ptr - other_buff.mem_ptr);
}

inline ParallelComm::Buffer::Buffer(unsigned int new_size) 
        : mem_ptr(NULL), buff_ptr(NULL), alloc_size(0)
{
  if (new_size) this->reserve(new_size);
}

inline ParallelComm::Buffer::~Buffer() 
{
  if (mem_ptr) {
    free(mem_ptr);
    mem_ptr = NULL;
  }
}

#define DEBUG_BUFFER 1

inline void ParallelComm::Buffer::reserve(unsigned int new_size) {
  
#ifdef DEBUG_BUFFER
  int tmp_pos = 0;
  if (mem_ptr) {
    tmp_pos = buff_ptr - mem_ptr;
  }
  MALLOC(buff_ptr, new_size, unsigned char *);
  assert(0 <= tmp_pos && tmp_pos <= (int)alloc_size);  
  if (tmp_pos) memcpy(buff_ptr, mem_ptr, tmp_pos);
  if (mem_ptr) free(mem_ptr);
  mem_ptr = buff_ptr;
  alloc_size = new_size;
  buff_ptr = mem_ptr + tmp_pos;
#else    
  if (mem_ptr && alloc_size < new_size) {
    size_t tmp_pos = mem_ptr ? buff_ptr - mem_ptr : 0;
    REALLOC(mem_ptr, mem_ptr, alloc_size, new_size, unsigned char*);
    alloc_size = new_size;
    buff_ptr = mem_ptr + tmp_pos;
  }
  else if (!mem_ptr) {
    MALLOC(mem_ptr, new_size, unsigned char *);
    alloc_size = new_size;
    buff_ptr = mem_ptr;
  } 
#endif
}

inline void ParallelComm::Buffer::check_space(unsigned int addl_space,
                                              unsigned int num_aligns)
{
#ifdef MEM_ALIGN
  addl_space += num_aligns*MEM_ALIGN;
#endif
  assert(buff_ptr >= mem_ptr && buff_ptr <= mem_ptr+alloc_size);
  unsigned int new_size = buff_ptr - mem_ptr + addl_space;
  if (new_size > alloc_size) 
    reserve(1.5*new_size);
}

inline void ParallelComm::reset_all_buffers() 
{
  std::vector<Buffer*>::iterator vit;
  for (vit = localOwnedBuffs.begin(); vit != localOwnedBuffs.end(); vit++)
    (*vit)->reset_buffer();
  for (vit = remoteOwnedBuffs.begin(); vit != remoteOwnedBuffs.end(); vit++)
    (*vit)->reset_buffer();
}

inline void ParallelComm::delete_all_buffers() 
{
  std::vector<Buffer*>::iterator vit;
  for (vit = localOwnedBuffs.begin(); vit != localOwnedBuffs.end(); vit++)
    delete (*vit);
  localOwnedBuffs.clear();
  
  for (vit = remoteOwnedBuffs.begin(); vit != remoteOwnedBuffs.end(); vit++)
    delete (*vit);
  remoteOwnedBuffs.clear();
}

inline std::vector<unsigned int> &ParallelComm::buff_procs() 
{
  return buffProcs;
}

inline ErrorCode ParallelComm::get_shared_proc_tags(Tag &sharedp,
                                                        Tag &sharedps,
                                                        Tag &sharedh,
                                                        Tag &sharedhs,
                                                        Tag &pstatus) 
{
  sharedp = sharedp_tag();
  sharedps = sharedps_tag();
  sharedh = sharedh_tag();
  sharedhs = sharedhs_tag();
  pstatus = pstatus_tag();
  
  return MB_SUCCESS;
}

inline ErrorCode ParallelComm::exchange_tags( const char *tag_name,
                                              const Range &entities)
{
    // get the tag handle
  std::vector<Tag> tags(1);
  ErrorCode result = mbImpl->tag_get_handle(tag_name, tags[0]);
  if (MB_SUCCESS != result) return result;
  else if (!tags[0]) return MB_TAG_NOT_FOUND;
  
  return exchange_tags(tags, tags, entities);
}
  
inline ErrorCode ParallelComm::exchange_tags( Tag tagh,
                                              const Range &entities)
{
    // get the tag handle
  std::vector<Tag> tags;
  tags.push_back(tagh);
  
  return exchange_tags(tags, tags, entities);
}
  
inline ErrorCode ParallelComm::get_comm_procs(std::set<unsigned int> &procs) 
{
  ErrorCode result = get_interface_procs(procs);
  if (MB_SUCCESS != result) return result;

  std::copy(buffProcs.begin(), buffProcs.end(), std::inserter(procs, procs.begin()));
    
  return MB_SUCCESS;
}

inline ErrorCode ParallelComm::get_owner(EntityHandle entity,
                                             int &owner) 
{
  EntityHandle tmp_handle;
  return get_owner_handle(entity, owner, tmp_handle);
}

    /* \brief Unpack message with remote handles (const pointer to buffer)
     */
inline ErrorCode ParallelComm::unpack_remote_handles(unsigned int from_proc,
                                                     const Buffer *buff_ptr,
                                                     std::vector<EntityHandle> &L2hloc,
                                                     std::vector<EntityHandle> &L2hrem,
                                                     std::vector<unsigned int> &L2p) 
{
    // cast away const-ness, we won't be passing back a modified ptr
  Buffer buff(*buff_ptr);
  
  return unpack_remote_handles(from_proc, &buff, 
                               L2hloc, L2hrem, L2p);
}

inline void ParallelComm::set_rank(unsigned int r) 
{
  procConfig.proc_rank(r);
  if (procConfig.proc_size() < r) procConfig.proc_size(r+1);
}

inline void ParallelComm::set_size(unsigned int s) 
{
  procConfig.proc_size(s);
}

inline ErrorCode ParallelComm::get_sharing_data(const EntityHandle *entities,
                                                    int num_entities,
                                                    std::set<int> &procs,
                                                    int op) 
{
  Range dum_range;
    // cast away constness 'cuz the range is passed as const
  EntityHandle *ents_cast = const_cast<EntityHandle*>(entities);
  std::copy(ents_cast, ents_cast+num_entities, range_inserter(dum_range));
  return get_sharing_data(dum_range, procs, op);
}

inline ErrorCode ParallelComm::get_sharing_data(const EntityHandle entity,
                                                    int *ps, 
                                                    EntityHandle *hs,
                                                    unsigned char &pstat,
                                                    int &num_ps) 
{
  unsigned int dum_ps;
  ErrorCode result = get_sharing_data(entity, ps, hs, pstat, dum_ps);
  if (MB_SUCCESS == result)
    num_ps = dum_ps;
  return result;
}

} // namespace moab
  
#endif
