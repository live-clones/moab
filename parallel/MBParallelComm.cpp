#include "MBInterface.hpp"
#include "MBParallelComm.hpp"
#include "MBWriteUtilIface.hpp"
#include "MBReadUtilIface.hpp"
#include "SequenceManager.hpp"
#include "EntitySequence.hpp"
#include "TagServer.hpp"
#include "MBTagConventions.hpp"
#include "MBSkinner.hpp"
#include "MBParallelConventions.h"
#include "MBCore.hpp"
#include "MBError.hpp"
#include "ElementSequence.hpp"
#include "MBCN.hpp"
#include "RangeMap.hpp"
#include "MeshTopoUtil.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>

#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)
const bool debug = false;

#include <math.h>
#include <assert.h>


extern "C" 
{
#include "minmax.h"
#include "gs.h"
#include "errmem.h"
#include "types.h"
#include "sort.h"
#include "tuple_list.h"
}

#ifdef USE_MPI
#include "mpi.h"
#endif

#define INITIAL_BUFF_SIZE 1024

#undef DEBUG_PACKING
#ifdef DEBUG_PACKING
unsigned int __PACK_num = 0, __UNPACK_num = 0, __PACK_count = 0, __UNPACK_count = 0;
std::string __PACK_string, __UNPACK_string;

#define PC(n, m) {\
          if (__PACK_num == (unsigned int)n && __PACK_string == m) __PACK_count++;\
          else {\
            if (__PACK_count > 1) std::cerr << " (" << __PACK_count << "x)";\
            __PACK_count = 1; __PACK_string = m; __PACK_num = n;\
            std::cerr << std::endl << "PACK: " << n << m;\
          }}
#define UPC(n, m) {\
          if (__UNPACK_num == (unsigned int)n && __UNPACK_string == m) __UNPACK_count++;\
          else {\
            if (__UNPACK_count > 1) std::cerr << "(" << __UNPACK_count << "x)";\
            __UNPACK_count = 1; __UNPACK_string = m; __UNPACK_num = n;\
            std::cerr << std::endl << "UNPACK: " << n << m;\
          }}
#else
#define PC(n, m)
#define UPC(n, m)
#endif

#define PACK_INT(buff, int_val) {int tmp_val = int_val; PACK_INTS(buff, &tmp_val, 1);}

#define PACK_INTS(buff, int_val, num) {memcpy(buff, int_val, (num)*sizeof(int)); buff += (num)*sizeof(int); PC(num, " ints");}

#define PACK_DBL(buff, dbl_val, num) {memcpy(buff, dbl_val, (num)*sizeof(double)); buff += (num)*sizeof(double); PC(num, " doubles");}

#define PACK_EH(buff, eh_val, num) {memcpy(buff, eh_val, (num)*sizeof(MBEntityHandle)); buff += (num)*sizeof(MBEntityHandle); PC(num, " handles");}

#define PACK_CHAR_64(buff, char_val) {strcpy((char*)buff, char_val); buff += 64; PC(64, " chars");}

#define PACK_VOID(buff, val, num) {memcpy(buff, val, num); buff += num; PC(num, " void");}

#define PACK_BYTES(buff, val, num) PACK_INT(buff, num) PACK_VOID(buff, val, num)

#define PACK_RANGE(buff, rng) {int num_subs = num_subranges(rng); PACK_INTS(buff, &num_subs, 1); PC(num_subs, "-subranged range"); \
          for (MBRange::const_pair_iterator cit = rng.const_pair_begin(); cit != rng.const_pair_end(); cit++) { \
            MBEntityHandle eh = (*cit).first; PACK_EH(buff, &eh, 1); \
            eh = (*cit).second; PACK_EH(buff, &eh, 1);}; }

#define UNPACK_INT(buff, int_val) {UNPACK_INTS(buff, &int_val, 1);}

#define UNPACK_INTS(buff, int_val, num) {memcpy(int_val, buff, (num)*sizeof(int)); buff += (num)*sizeof(int); UPC(num, " ints");}

#define UNPACK_DBL(buff, dbl_val, num) {memcpy(dbl_val, buff, (num)*sizeof(double)); buff += (num)*sizeof(double); UPC(num, " doubles");}

#define UNPACK_EH(buff, eh_val, num) {memcpy(eh_val, buff, (num)*sizeof(MBEntityHandle)); buff += (num)*sizeof(MBEntityHandle); UPC(num, " handles");}

#define UNPACK_CHAR_64(buff, char_val) {strcpy(char_val, (char*)buff); buff += 64; UPC(64, " chars");}

#define UNPACK_VOID(buff, val, num) {memcpy(val, buff, num); buff += num; UPC(num, " void");}

#define UNPACK_RANGE(buff, rng) {int num_subs; UNPACK_INTS(buff, &num_subs, 1); UPC(num_subs, "-subranged range"); MBEntityHandle _eh[2]; \
          for (int i = 0; i < num_subs; i++) { UNPACK_EH(buff, _eh, 2); rng.insert(_eh[0], _eh[1]);}}
#define CHECK_BUFF_SPACE(buff_vec, buff_ptr, addl_space) { \
      unsigned int _old_size = buff_ptr - &buff_vec[0],            \
          _new_size = _old_size + (addl_space);    \
      if (_new_size > buff_vec.size()) {               \
        buff_vec.resize(1.5*_new_size);            \
        buff_ptr = &buff_vec[_new_size-(addl_space)];} }
    

#define RANGE_SIZE(rng) (2*sizeof(MBEntityHandle)*num_subranges(rng)+sizeof(int))
#define RR(a) if (MB_SUCCESS != result) {\
          dynamic_cast<MBCore*>(mbImpl)->get_error_handler()->set_last_error(a);\
          return result;}

#define RRA(a) if (MB_SUCCESS != result) {\
      std::string tmp_str; mbImpl->get_last_error(tmp_str);\
      tmp_str.append("\n"); tmp_str.append(a);\
      dynamic_cast<MBCore*>(mbImpl)->get_error_handler()->set_last_error(tmp_str.c_str()); \
      return result;}

#define RRAI(i, a) if (MB_SUCCESS != result) {                \
      std::string tmp_str; i->get_last_error(tmp_str);\
      tmp_str.append("\n"); tmp_str.append(a);\
      dynamic_cast<MBCore*>(i)->get_error_handler()->set_last_error(tmp_str.c_str()); \
      return result;}

/** Name of tag used to store MBParallelComm Index on mesh paritioning sets */
const char* PARTITIONING_PCOMM_TAG_NAME = "__PRTN_PCOMM";
 
/** \brief Tag storing parallel communication objects
 *
 * This tag stores pointers to MBParallelComm communication
 * objects; one of these is allocated for each different
 * communicator used to read mesh.  MBParallelComm stores
 * partition and interface sets corresponding to its parallel mesh.
 * By default, a parallel read uses the first MBParallelComm object
 * on the interface instance; if instantiated with one, ReadParallel
 * adds this object to the interface instance too.
 *
 * Tag type: opaque
 * Tag size: MAX_SHARING_PROCS*sizeof(MBParallelComm*)
 */
#define PARALLEL_COMM_TAG_NAME "__PARALLEL_COMM"


enum MBMessageTag {MB_MESG_ANY=MPI_ANY_TAG, 
                   MB_MESG_SIZE,
                   MB_MESG_ENTS,
                   MB_MESG_REMOTE_HANDLES,
                   MB_MESG_SHAREDHPS,
                   MB_MESG_TAGS };
    
MBParallelComm::MBParallelComm(MBInterface *impl, MPI_Comm comm, int* id ) 
        : mbImpl(impl), procConfig(comm),
          sharedpTag(0), sharedpsTag(0), 
          sharedhTag(0), sharedhsTag(0), pstatusTag(0), ifaceSetsTag(0),
          partitionTag(0), globalPartCount(-1), partitioningSet(0)
{
  myBuffer.resize(INITIAL_BUFF_SIZE);

  tagServer = dynamic_cast<MBCore*>(mbImpl)->tag_server();
  sequenceManager = dynamic_cast<MBCore*>(mbImpl)->sequence_manager();

  int flag = 1;
  int retval = MPI_Initialized(&flag);
  if (MPI_SUCCESS != retval || !flag) {
    int argc = 0;
    char **argv = NULL;
    
      // mpi not initialized yet - initialize here
    retval = MPI_Init(&argc, &argv);
  }

  pcommID = add_pcomm(this);
  if (id)
    *id = pcommID;
}

MBParallelComm::MBParallelComm(MBInterface *impl,
                               std::vector<unsigned char> &tmp_buff, 
                               MPI_Comm comm,
                               int* id) 
    : mbImpl(impl), procConfig(comm),
      sharedpTag(0), sharedpsTag(0), 
      sharedhTag(0), sharedhsTag(0), pstatusTag(0), ifaceSetsTag(0),
      partitionTag(0), globalPartCount(-1), partitioningSet(0)
{
  myBuffer.swap(tmp_buff);
  int flag = 1;
  int retval = MPI_Initialized(&flag);
  if (MPI_SUCCESS != retval || !flag) {
    int argc = 0;
    char **argv = NULL;
    
      // mpi not initialized yet - initialize here
    retval = MPI_Init(&argc, &argv);
  }

  pcommID = add_pcomm(this);
  if (id)
    *id = pcommID;
}

MBParallelComm::~MBParallelComm() 
{
  remove_pcomm(this);
}

int MBParallelComm::add_pcomm(MBParallelComm *pc) 
{
    // add this pcomm to instance tag
  std::vector<MBParallelComm *> pc_array(MAX_SHARING_PROCS, 
                                         (MBParallelComm*)NULL);
  MBTag pc_tag = pcomm_tag(mbImpl, true);
  assert(0 != pc_tag);
  
  MBErrorCode result = mbImpl->tag_get_data(pc_tag, 0, 0, (void*)&pc_array[0]);
  if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) 
    return -1;
  int index = 0;
  while (index < MAX_SHARING_PROCS && pc_array[index]) index++;
  if (index == MAX_SHARING_PROCS) {
    index = -1;
    assert(false);
  }
  else {
    pc_array[index] = pc;
    mbImpl->tag_set_data(pc_tag, 0, 0, (void*)&pc_array[0]);
  }
  return index;
}

void MBParallelComm::remove_pcomm(MBParallelComm *pc) 
{
    // remove this pcomm from instance tag
  std::vector<MBParallelComm *> pc_array(MAX_SHARING_PROCS);
  MBTag pc_tag = pcomm_tag(mbImpl, true);
  
  MBErrorCode result = mbImpl->tag_get_data(pc_tag, 0, 0, (void*)&pc_array[0]);
  std::vector<MBParallelComm*>::iterator pc_it = 
    std::find(pc_array.begin(), pc_array.end(), pc);
  assert(MB_SUCCESS == result && 
         pc_it != pc_array.end());
  *pc_it = NULL;
  mbImpl->tag_set_data(pc_tag, 0, 0, (void*)&pc_array[0]);
}

//! assign a global id space, for largest-dimension or all entities (and
//! in either case for vertices too)
MBErrorCode MBParallelComm::assign_global_ids(MBEntityHandle this_set,
                                              const int dimension, 
                                              const int start_id,
                                              const bool largest_dim_only,
                                              const bool parallel) 
{
  MBRange entities[4];
  int local_num_elements[4];
  MBErrorCode result;
  std::vector<unsigned char> pstatus;
  for (int dim = 0; dim <= dimension; dim++) {
    if (dim == 0 || !largest_dim_only || dim == dimension) {
      result = mbImpl->get_entities_by_dimension(this_set, dim, entities[dim]); 
      RRA("Failed to get vertices in assign_global_ids.");
    }

      // need to filter out non-locally-owned entities!!!
    pstatus.resize(entities[dim].size());
    result = mbImpl->tag_get_data(pstatus_tag(), entities[dim], &pstatus[0]);
    RRA("Failed to get pstatus in assign_global_ids.");
    
    MBRange dum_range;
    MBRange::iterator rit;
    unsigned int i;
    for (rit = entities[dim].begin(), i = 0; rit != entities[dim].end(); rit++, i++)
      if (pstatus[i] & PSTATUS_NOT_OWNED)
        dum_range.insert(*rit);
    entities[dim] = entities[dim].subtract(dum_range);
    
    local_num_elements[dim] = entities[dim].size();
  }
  
    // communicate numbers
  std::vector<int> num_elements(procConfig.proc_size()*4);
#ifdef USE_MPI
  if (procConfig.proc_size() > 1 && parallel) {
    int retval = MPI_Alltoall(local_num_elements, 4, MPI_INT,
                              &num_elements[0], procConfig.proc_size()*4, 
                              MPI_INT, procConfig.proc_comm());
    if (0 != retval) return MB_FAILURE;
  }
  else
#endif
    for (int dim = 0; dim < 4; dim++) num_elements[dim] = local_num_elements[dim];
  
    // my entities start at one greater than total_elems[d]
  int total_elems[4] = {start_id, start_id, start_id, start_id};
  
  for (unsigned int proc = 0; proc < procConfig.proc_rank(); proc++) {
    for (int dim = 0; dim < 4; dim++) total_elems[dim] += num_elements[4*proc + dim];
  }
  
    //assign global ids now
  MBTag gid_tag;
  int zero = 0;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), 
                              MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                              &zero, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) return result;
  
  for (int dim = 0; dim < 4; dim++) {
    if (entities[dim].empty()) continue;
    num_elements.resize(entities[dim].size());
    int i = 0;
    for (MBRange::iterator rit = entities[dim].begin(); rit != entities[dim].end(); rit++)
      num_elements[i++] = total_elems[dim]++;
    
    result = mbImpl->tag_set_data(gid_tag, entities[dim], &num_elements[0]); 
    RRA("Failed to set global id tag in assign_global_ids.");
  }
  
  return MB_SUCCESS;
}

int MBParallelComm::get_buffers(int to_proc, bool *is_new) 
{
  int ind = -1;
  std::vector<unsigned int>::iterator vit = 
    std::find(buffProcs.begin(), buffProcs.end(), to_proc);
  if (vit == buffProcs.end()) {
    ind = buffProcs.size();
    buffProcs.push_back((unsigned int)to_proc);
    ownerSBuffs.push_back(std::vector<unsigned char>());
    ghostRBuffs.push_back(std::vector<unsigned char>(INITIAL_BUFF_SIZE));
      // allocate these other buffs in case we're storing remote handles
    ownerRBuffs.push_back(std::vector<unsigned char>(INITIAL_BUFF_SIZE));
    ghostSBuffs.push_back(std::vector<unsigned char>());
    if (is_new) *is_new = true;
  }
  else {
    ind = vit - buffProcs.begin();
    if (is_new) *is_new = false;
  }
  assert(ind < MAX_SHARING_PROCS);
  return ind;
}

MBErrorCode MBParallelComm::send_buffer(const unsigned int to_proc,
                                        const unsigned char *send_buff,
                                        const unsigned int buff_size,
                                        const int msg_type,
                                        MPI_Request &send_req) 
{
#ifndef USE_MPI
  return MB_FAILURE;
#else

  MBErrorCode result = MB_SUCCESS;
  int success;

    // if the message is large, send a first message to tell how large
  if (INITIAL_BUFF_SIZE < buff_size) {
    int tmp_buff_size = -buff_size;
    int success = MPI_Isend(&tmp_buff_size, sizeof(int), MPI_UNSIGNED_CHAR, 
                            to_proc, MB_MESG_SIZE, procConfig.proc_comm(),
                            &send_req);
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }
    
    // send the buffer
  success = MPI_Isend(const_cast<unsigned char*>(send_buff), buff_size, MPI_UNSIGNED_CHAR, to_proc, 
                      msg_type, procConfig.proc_comm(), &send_req);
  if (success != MPI_SUCCESS) return MB_FAILURE;

  return result;
#endif
}

MBErrorCode MBParallelComm::recv_size_buff(const int from_proc,
                                           std::vector<unsigned char> &recv_buff,
                                           MPI_Request &recv_req,
                                           int mesg_tag) 
{
    // use the received size to resize buffer, then post another irecv
  recv_buff.resize(-(*((int*)&recv_buff[0])));
  int success = MPI_Irecv(&recv_buff[0], recv_buff.size(), MPI_UNSIGNED_CHAR, from_proc, 
                          mesg_tag, procConfig.proc_comm(), &recv_req);
  if (MPI_SUCCESS != success) {
    MBErrorCode result = MB_FAILURE;
    RRA("Failed call to Irecv in recv_size_buff.");
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::broadcast_entities( const int from_proc,
                                                MBRange &entities,
                                                const bool adjacencies,
                                                const bool tags) 
{
#ifndef USE_MPI
  return MB_FAILURE;
#else
  
  MBErrorCode result = MB_SUCCESS;
  int success;
  int buff_size;
  
  std::vector<unsigned char> buff;
  std::vector<int> addl_procs;
  if ((int)procConfig.proc_rank() == from_proc) {
    result = add_verts(entities);
    RRA("Failed to add adj vertices.");

    result = pack_buffer( entities, adjacencies, tags, 
                          false, -1, buff, buff_size); 
    RRA("Failed to compute buffer size in broadcast_entities.");
  }

  success = MPI_Bcast( &buff_size, 1, MPI_INT, from_proc, procConfig.proc_comm() );
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("MPI_Bcast of buffer size failed.");
  }
  
  if (!buff_size) // no data
    return MB_SUCCESS;

  if ((int)procConfig.proc_rank() != from_proc) 
    buff.resize(buff_size);

  success = MPI_Bcast( &buff[0], buff_size, MPI_UNSIGNED_CHAR, from_proc, procConfig.proc_comm() );
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("MPI_Bcast of buffer failed.");
  }

  if ((int)procConfig.proc_rank() != from_proc) {
    std::vector<std::vector<MBEntityHandle> > dum1a, dum1b;
    std::vector<std::vector<int> > dum1p;
    std::vector<MBEntityHandle> dum2;
    std::vector<unsigned int> dum3;
    result = unpack_buffer(&buff[0], false, from_proc, -1, 
                           dum1a, dum1b, dum1p, dum2, dum2, dum3, entities);
    RRA("Failed to unpack buffer in broadcast_entities.");
  }

  return MB_SUCCESS;
#endif
}

MBErrorCode MBParallelComm::pack_buffer(MBRange &orig_ents, 
                                        const bool adjacencies,
                                        const bool tags,
                                        const bool store_remote_handles,
                                        const int to_proc,
                                        std::vector<unsigned char> &buff,
                                        int &buff_size) 
{
    // pack the buffer with the entity ranges, adjacencies, and tags sections
    // 
    // Note: new entities used in subsequent connectivity lists, sets, or tags, 
    //   are referred to as (MBMAXTYPE + index), where index is into vector 
    //   of new entities, 0-based
  MBErrorCode result;

  MBRange set_range;
  std::vector<MBRange> set_ranges;
  std::vector<MBTag> all_tags;
  std::vector<MBRange> tag_ranges;
  std::vector<int> set_sizes;
  std::vector<unsigned int> options_vec;

  MBRange::const_iterator rit;

    // get an estimate of the buffer size
  buff_size = estimate_ents_buffer_size(orig_ents, store_remote_handles);
  buff.clear();
  buff.reserve(1);
  unsigned char *buff_ptr = &buff[0];
  CHECK_BUFF_SPACE(buff, buff_ptr, buff_size);
  
    // entities
  result = pack_entities(orig_ents, buff, buff_ptr,
                         store_remote_handles, to_proc, false); 
  RRA("Packing entities failed.");
  
    // sets
  result = pack_sets(orig_ents, buff, buff_ptr,
                     store_remote_handles, to_proc); 
  RRA("Packing sets (count) failed.");

    // tags
  MBRange final_ents;
  if (tags) {
    result = get_tag_send_list(orig_ents, all_tags, tag_ranges );
    RRA("Failed to get tagged entities.");
    result = pack_tags(orig_ents, all_tags, all_tags, tag_ranges, 
                       buff, buff_ptr, store_remote_handles, to_proc);
    RRA("Packing tags (count) failed.");
  }

  return result;
}
 
MBErrorCode MBParallelComm::unpack_buffer(unsigned char *buff_ptr,
                                          const bool store_remote_handles,
                                          const int from_proc,
                                          const int ind,
                                          std::vector<std::vector<MBEntityHandle> > &L1hloc,
                                          std::vector<std::vector<MBEntityHandle> > &L1hrem,
                                          std::vector<std::vector<int> > &L1p,
                                          std::vector<MBEntityHandle> &L2hloc, 
                                          std::vector<MBEntityHandle> &L2hrem,
                                          std::vector<unsigned int> &L2p,
                                          MBRange &new_ents) 
{
  if (myBuffer.capacity() == 0) return MB_FAILURE;
  
#ifdef DEBUG_PACKING
    unsigned char *tmp_buff = buff_ptr;
#endif  
    MBErrorCode result;
    result = unpack_entities(buff_ptr, store_remote_handles,
                             ind, false, L1hloc, L1hrem, L1p, L2hloc, L2hrem, L2p, new_ents);
  RRA("Unpacking entities failed.");
#ifdef DEBUG_PACKING
    std::cerr << "unpack_entities buffer space: " << buff_ptr - tmp_buff << " bytes." << std::endl;
    tmp_buff = buff_ptr;
#endif  
  result = unpack_sets(buff_ptr, new_ents, store_remote_handles, from_proc);
  RRA("Unpacking sets failed.");
#ifdef DEBUG_PACKING
    std::cerr << "unpack_sets buffer space: " << buff_ptr - tmp_buff << " bytes." << std::endl;
    tmp_buff = buff_ptr;
#endif  
  result = unpack_tags(buff_ptr, new_ents, store_remote_handles, from_proc);
  RRA("Unpacking tags failed.");
#ifdef DEBUG_PACKING
    std::cerr << "unpack_tags buffer space: " << buff_ptr - tmp_buff << " bytes." << std::endl;
    tmp_buff = buff_ptr;
#endif  

#ifdef DEBUG_PACKING
  std::cerr << std::endl;
#endif
  
  return MB_SUCCESS;
}

int MBParallelComm::num_subranges(const MBRange &this_range)
{
    // ok, have all the ranges we'll pack; count the subranges
  int num_sub_ranges = 0;
  for (MBRange::const_pair_iterator pit = this_range.const_pair_begin(); 
       pit != this_range.const_pair_end(); pit++)
    num_sub_ranges++;

  return num_sub_ranges;
}

int MBParallelComm::estimate_ents_buffer_size(MBRange &entities,
                                              const bool store_remote_handles) 
{
  int buff_size = 0;
  std::vector<MBEntityHandle> dum_connect_vec;
  const MBEntityHandle *connect;
  int num_connect;

  int num_verts = entities.num_of_type(MBVERTEX);
    // # verts + coords + handles
  buff_size += 2*sizeof(int) + 3*sizeof(double)*num_verts;
  if (store_remote_handles) buff_size += sizeof(MBEntityHandle)*num_verts;

    // do a rough count by looking at first entity of each type
  for (MBEntityType t = MBEDGE; t < MBENTITYSET; t++) {
    const MBRange::iterator rit = entities.lower_bound(t);
    if (TYPE_FROM_HANDLE(*rit) != t) continue;
    
    MBErrorCode result = mbImpl->get_connectivity(*rit, connect, num_connect, 
                                                  true, &dum_connect_vec);
    RRA("Failed to get connectivity to estimate buffer size.");

      // number, type, nodes per entity
    buff_size += 3*sizeof(int);
    int num_ents = entities.num_of_type(t);
      // connectivity, handle for each ent
    buff_size += (num_connect+1)*sizeof(MBEntityHandle)*num_ents;
  }

      // extra entity type at end, passed as int
  buff_size += sizeof(int);

  return buff_size;
}

int MBParallelComm::estimate_sets_buffer_size(MBRange &entities,
                                              const bool store_remote_handles) 
{
    // number of sets
  int buff_size = sizeof(int);
  
    // do a rough count by looking at first entity of each type
  MBRange::iterator rit = entities.lower_bound(MBENTITYSET);
  MBErrorCode result;
  
  for (; rit != entities.end(); rit++) {
    unsigned int options;
    result = mbImpl->get_meshset_options(*rit, options);
    RRA("Failed to get meshset options.");

    buff_size += sizeof(int);
    
    MBRange set_range;
    if (options & MESHSET_SET) {
        // range-based set; count the subranges
      result = mbImpl->get_entities_by_handle(*rit, set_range);
      RRA("Failed to get set entities.");

        // set range
      buff_size += RANGE_SIZE(set_range);
    }
    else if (options & MESHSET_ORDERED) {
        // just get the number of entities in the set
      int num_ents;
      result = mbImpl->get_number_entities_by_handle(*rit, num_ents);
      RRA("Failed to get number entities in ordered set.");

        // set vec
      buff_size += sizeof(MBEntityHandle) * num_ents + sizeof(int);
    }

      // get numbers of parents/children
    int num_par, num_ch;
    result = mbImpl->num_child_meshsets(*rit, &num_ch);
    RRA("Failed to get num children.");

    result = mbImpl->num_parent_meshsets(*rit, &num_par);
    RRA("Failed to get num parents.");

    buff_size += (num_ch + num_par) * sizeof(MBEntityHandle) + 2*sizeof(int);
  }

  return buff_size;
}

MBErrorCode MBParallelComm::pack_entities(MBRange &entities,
                                          std::vector<unsigned char> &buff,
                                          unsigned char *&buff_ptr,
                                          const bool store_remote_handles,
                                          const int to_proc,
                                          const bool is_iface,
                                          std::vector<std::set<unsigned int> > *entprocs,
                                          MBRange *allsent) 
{
    // packed information:
    // 1. # entities = E
    // 2. for e in E
    //   a. # procs sharing e, incl. sender and receiver = P
    //   b. for p in P (procs sharing e)
    //   c. for p in P (handle for e on p) (Note1)
    // 3. vertex/entity info

    // get an estimate of the buffer size & pre-allocate buffer size
  unsigned int buff_size = estimate_ents_buffer_size(entities, 
                                                     store_remote_handles);
  CHECK_BUFF_SPACE(buff, buff_ptr, buff_size);
  
  MBWriteUtilIface *wu;
  MBErrorCode result = mbImpl->query_interface(std::string("MBWriteUtilIface"), 
                                               reinterpret_cast<void**>(&wu));
  RRA("Couldn't get MBWriteUtilIface.");

  unsigned int num_ents;

    // first pack procs/handles sharing this ent, not including this dest but including
    // others (with zero handles)
  if (store_remote_handles) {

      // buff space is at least proc+handle for each entity; use avg of 4 other procs
      // to estimate buff size, but check later
    CHECK_BUFF_SPACE(buff, buff_ptr,
                     sizeof(int) + (5*sizeof(int) + sizeof(MBEntityHandle))*entities.size());

      // 1. # entities = E
    PACK_INT(buff_ptr, entities.size());
  
    MBRange::iterator rit;
  
      // pre-fetch sharedp and pstatus
    std::vector<int> sharedp_vals(entities.size());
    result = mbImpl->tag_get_data(sharedp_tag(), entities, &sharedp_vals[0]);
    RRA("Failed to get sharedp_tag.");
    std::vector<char> pstatus_vals(entities.size());
    result = mbImpl->tag_get_data(pstatus_tag(), entities, &pstatus_vals[0]);
    RRA("Failed to get sharedp_tag.");
  
    unsigned int i;
    std::vector<std::set<unsigned int> >::iterator vit;
    int tmp_procs[MAX_SHARING_PROCS];
    MBEntityHandle tmp_handles[MAX_SHARING_PROCS];
    std::set<unsigned int> dumprocs;

      // 2. for e in E
    for (rit = entities.begin(), i = 0; 
         rit != entities.end(); rit++, i++) {
      result = build_sharedhps_list(*rit, pstatus_vals[i], sharedp_vals[i],
                                    (entprocs ? (*entprocs)[allsent->index(*rit)] : dumprocs),
                                    num_ents, tmp_procs, tmp_handles);
      RRA("Failed to build sharedhps.");

        // now pack them
      CHECK_BUFF_SPACE(buff, buff_ptr, (num_ents+1)*sizeof(int) + sizeof(MBEntityHandle));
      PACK_INT(buff_ptr, num_ents);
      PACK_INTS(buff_ptr, tmp_procs, num_ents);
      PACK_EH(buff_ptr, tmp_handles, num_ents);
    }
  }
  
    // pack vertices
  MBRange these_ents = entities.subset_by_type(MBVERTEX);
  num_ents = these_ents.size();

  if (num_ents) {
    buff_size = 2*sizeof(int) + 3*num_ents*sizeof(double);
    CHECK_BUFF_SPACE(buff, buff_ptr, buff_size);

    // type, # ents
    PACK_INT(buff_ptr, ((int) MBVERTEX));
    PACK_INT(buff_ptr, ((int) num_ents));

    result = mbImpl->get_coords(these_ents, (double*)buff_ptr);
    RRA("Couldn't get vertex coordinates.");
    PC(3*num_ents, " doubles");

    buff_ptr += 3 * num_ents * sizeof(double);

#ifdef DEBUG_PACKING
  std::cerr << "Packed " << these_ents.size() << " ents of type " 
            << MBCN::EntityTypeName(TYPE_FROM_HANDLE(*these_ents.begin())) << std::endl;
#endif      
  }

    // now entities; go through range, packing by type and equal # verts per element
  MBRange::iterator start_rit = entities.find(*these_ents.rbegin());
  start_rit++;
  int last_nodes = -1;
  MBEntityType last_type = MBMAXTYPE;
  these_ents.clear();
  MBRange::iterator end_rit = start_rit;
  EntitySequence *seq;
  ElementSequence *eseq;
  
  while (start_rit != entities.end() || !these_ents.empty()) {
      // cases:
      // A: !end, last_type == MAXTYPE, seq: save contig sequence in these_ents
      // B: !end, last type & nodes same, seq: save contig sequence in these_ents
      // C: !end, last type & nodes different: pack these_ents, then save contig sequence in these_ents
      // D: end: pack these_ents

      // find the sequence holding current start entity, if we're not at end
    eseq = NULL;
    if (start_rit != entities.end()) {
      result = sequenceManager->find(*start_rit, seq);
      RRA("Couldn't find entity sequence.");
      if (NULL == seq) return MB_FAILURE;
      eseq = dynamic_cast<ElementSequence*>(seq);
    }

      // pack the last batch if at end or next one is different
    if (!these_ents.empty() &&
        (!eseq || eseq->type() != last_type ||
         last_nodes != (int) eseq->nodes_per_element())) {
      result = pack_entity_seq(last_nodes, store_remote_handles,
                               to_proc, these_ents, entities, buff, buff_ptr);
      RRA("Failed to pack entities from a sequence.");
      these_ents.clear();
    }

    if (eseq) {
        // continuation of current range, just save these entities
        // get position in entities list one past end of this sequence
      end_rit = entities.lower_bound(start_rit, entities.end(), eseq->end_handle()+1);

        // put these entities in the range
      std::copy(start_rit, end_rit, mb_range_inserter(these_ents));

      last_type = eseq->type();
      last_nodes = eseq->nodes_per_element();
    }
    else if (start_rit != entities.end() &&
             TYPE_FROM_HANDLE(*start_rit) == MBENTITYSET)
      break;

    start_rit = end_rit;
  }

    // pack MBMAXTYPE to indicate end of ranges
  CHECK_BUFF_SPACE(buff, buff_ptr, sizeof(int));
  PACK_INT(buff_ptr, ((int)MBMAXTYPE));

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::build_sharedhps_list(const MBEntityHandle entity,
                                                 const unsigned char pstatus,
                                                 const int sharedp, 
                                                 const std::set<unsigned int> &entprocs,
                                                 unsigned int &num_ents,
                                                 int *tmp_procs,
                                                 MBEntityHandle *tmp_handles)
{
  num_ents = 0;
  unsigned char pstat;
  MBErrorCode result = get_sharing_data(entity, tmp_procs, tmp_handles,
                                        pstat, num_ents);
  RRA("Failed in get_sharing_data.");
  assert(pstat == pstatus);
  
    // build shared proc/handle lists
    // start with multi-shared, since if it is the owner will be first
  if (pstatus & PSTATUS_MULTISHARED) {
  }
  else if (pstatus & PSTATUS_NOT_OWNED) {
      // if not multishared and not owned, other sharing proc is owner, put that
      // one first
    assert("If not owned, I should be shared too" &&
           pstatus & PSTATUS_SHARED &&
           num_ents == 1);
    tmp_procs[1] = procConfig.proc_rank();
    tmp_handles[1] = entity;
    num_ents = 2;
  }
  else if (pstatus & PSTATUS_SHARED) {
      // if not multishared and owned, I'm owner
    assert("shared and owned, should be only 1 sharing proc" &&
           1 == num_ents);
    tmp_procs[1] = tmp_procs[0];
    tmp_procs[0] = procConfig.proc_rank();
    tmp_handles[1] = tmp_handles[0];
    tmp_handles[0] = entity;
    num_ents = 2;
  }
  else {
      // not shared yet, just add owner (me)
    tmp_procs[0] = procConfig.proc_rank();
    tmp_handles[0] = entity;
    num_ents = 1;
  }

  int tmp_ps = num_ents;
  
    // now add others, with zero handle for now
  for (std::set<unsigned int>::iterator sit = entprocs.begin();
       sit != entprocs.end(); sit++) {
    assert("these procs shouldn't already be in the shared list" &&
           std::find(tmp_procs, tmp_procs+tmp_ps, *sit) == tmp_procs+tmp_ps);
    tmp_procs[num_ents] = *sit;
    tmp_handles[num_ents] = 0;
    num_ents++;
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_entity_seq(const int nodes_per_entity,
                                            const bool store_remote_handles,
                                            const int to_proc,
                                            MBRange &these_ents,
                                            MBRange &entities,
                                            std::vector<unsigned char> &buff,
                                            unsigned char *&buff_ptr) 
{
  int tmp_space = 3*sizeof(int) + nodes_per_entity*these_ents.size()*sizeof(MBEntityHandle);
  CHECK_BUFF_SPACE(buff, buff_ptr, tmp_space);
  
    // pack the entity type
  PACK_INT(buff_ptr, ((int)TYPE_FROM_HANDLE(*these_ents.begin())));

    // pack # ents
  PACK_INT(buff_ptr, these_ents.size());
      
    // pack the nodes per entity
  PACK_INT(buff_ptr, nodes_per_entity);
      
    // pack the connectivity
  const MBEntityHandle *connect;
  int num_connect;
  std::vector<MBEntityHandle> dum_connect;
  MBEntityHandle *start_vec = (MBEntityHandle*)buff_ptr;
  MBErrorCode result = MB_SUCCESS;
  for (MBRange::const_iterator rit = these_ents.begin(); rit != these_ents.end(); rit++) {
    result = mbImpl->get_connectivity(*rit, connect, num_connect, false,
                                      &dum_connect);
    RRA("Failed to get connectivity.");
    assert(num_connect == nodes_per_entity);
    PACK_EH(buff_ptr, connect, num_connect);
  }

    // substitute destination handles
  result = get_remote_handles(store_remote_handles, start_vec, start_vec,
                              nodes_per_entity*these_ents.size(), to_proc,
                              entities);
  RRA("Trouble getting remote handles when packing entities.");

#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Packed " << these_ents.size() << " ents of type " 
            << MBCN::EntityTypeName(TYPE_FROM_HANDLE(*these_ents.begin())) << std::endl;
#endif      

  return result;
}


MBErrorCode MBParallelComm::get_remote_handles(const bool store_remote_handles,
                                               MBEntityHandle *from_vec, 
                                               MBEntityHandle *to_vec_tmp,
                                               int num_ents, int to_proc,
                                               const MBRange &new_ents) 
{
    // NOTE: THIS IMPLEMENTATION IS JUST LIKE THE RANGE-BASED VERSION, NO REUSE
    // AT THIS TIME, SO IF YOU FIX A BUG IN THIS VERSION, IT MAY BE IN THE
    // OTHER VERSION TOO!!!
  if (0 == num_ents) return MB_SUCCESS;
  
    // use a local destination ptr in case we're doing an in-place copy
  std::vector<MBEntityHandle> tmp_vector;
  MBEntityHandle *to_vec = to_vec_tmp;
  if (to_vec == from_vec) {
    tmp_vector.resize(num_ents);
    to_vec = &tmp_vector[0];
  }

  if (!store_remote_handles) {
    int err;
      // in this case, substitute position in new_ents list
    for (int i = 0; i < num_ents; i++) {
      int ind = new_ents.index(from_vec[i]);
      to_vec[i] = CREATE_HANDLE(MBMAXTYPE, ind, err);
      assert(to_vec[i] != 0 && !err && -1 != ind);
    }
  }
  else {
    MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
    MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                              sharedh_tag, sharedhs_tag, pstatus_tag);
  
      // get single-proc destination handles and shared procs
    std::vector<int> sharing_procs(num_ents);
    result = mbImpl->tag_get_data(sharedh_tag, from_vec, num_ents,
                                  to_vec);
    RRA("Failed to get shared handle tag for remote_handles.");
    result = mbImpl->tag_get_data(sharedp_tag, from_vec, num_ents, &sharing_procs[0]);
    RRA("Failed to get sharing proc tag in remote_handles.");
    for (int j = 0; j < num_ents; j++) {
      if (to_vec[j] && sharing_procs[j] != to_proc)
        to_vec[j] = 0;
    }
    
    MBEntityHandle tmp_handles[MAX_SHARING_PROCS];
    int tmp_procs[MAX_SHARING_PROCS];
    int i;
      // go through results, and for 0-valued ones, look for multiple shared proc
    MBEntityHandle *tmp_eh;
    for (tmp_eh = to_vec, i = 0; i < num_ents; i++) {
      if (!to_vec[i]) {
        result = mbImpl->tag_get_data(sharedps_tag, from_vec+i, 1, tmp_procs);
        if (MB_SUCCESS == result) {
          for (int j = 0; j < MAX_SHARING_PROCS; j++) {
            if (-1 == tmp_procs[j]) break;
            else if (tmp_procs[j] == to_proc) {
              result = mbImpl->tag_get_data(sharedhs_tag, from_vec+i, 1, tmp_handles);
              RRA("Trouble getting sharedhs tag.");
              to_vec[i] = tmp_handles[j];
              assert(to_vec[i]);
              break;
            }
          }
        }
        if (!to_vec[i]) {
          int j = new_ents.index(from_vec[i]);
          if (-1 == j) {
            result = MB_FAILURE;
            RRA("Failed to find new entity in send list.");
          }
          int err;
          to_vec[i] = CREATE_HANDLE(MBMAXTYPE, j, err);
          if (err) {
            result = MB_FAILURE;
            RRA("Failed to create handle in remote_handles.");
          }
        }
      }
    }
  }
  
    // memcpy over results if from_vec and to_vec are the same
  if (to_vec_tmp == from_vec) 
    memcpy(from_vec, to_vec, num_ents * sizeof(MBEntityHandle));
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_remote_handles(const bool store_remote_handles,
                                               const MBRange &from_range, 
                                               MBEntityHandle *to_vec,
                                               int to_proc,
                                               const MBRange &new_ents) 
{
    // NOTE: THIS IMPLEMENTATION IS JUST LIKE THE VECTOR-BASED VERSION, NO REUSE
    // AT THIS TIME, SO IF YOU FIX A BUG IN THIS VERSION, IT MAY BE IN THE
    // OTHER VERSION TOO!!!
  if (from_range.empty()) return MB_SUCCESS;
  
  if (!store_remote_handles) {
    int err;
      // in this case, substitute position in new_ents list
    MBRange::iterator rit;
    unsigned int i;
    for (rit = from_range.begin(), i = 0; rit != from_range.end(); rit++, i++) {
      int ind = new_ents.index(*rit);
      to_vec[i] = CREATE_HANDLE(MBMAXTYPE, ind, err);
      assert(to_vec[i] != 0 && !err && -1 != ind);
    }
  }
  else {
    MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
    MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                              sharedh_tag, sharedhs_tag, pstatus_tag);
  
      // get single-proc destination handles and shared procs
    std::vector<int> sharing_procs(from_range.size());
    result = mbImpl->tag_get_data(sharedh_tag, from_range, to_vec);
    RRA("Failed to get shared handle tag for remote_handles.");
    result = mbImpl->tag_get_data(sharedp_tag, from_range, &sharing_procs[0]);
    RRA("Failed to get sharing proc tag in remote_handles.");
    for (unsigned int j = 0; j < from_range.size(); j++) {
      if (to_vec[j] && sharing_procs[j] != to_proc)
        to_vec[j] = 0;
    }
    
    MBEntityHandle tmp_handles[MAX_SHARING_PROCS];
    int tmp_procs[MAX_SHARING_PROCS];
      // go through results, and for 0-valued ones, look for multiple shared proc
    MBRange::iterator rit;
    unsigned int i;
    for (rit = from_range.begin(), i = 0; rit != from_range.end(); rit++, i++) {
      if (!to_vec[i]) {
        result = mbImpl->tag_get_data(sharedhs_tag, &(*rit), 1, tmp_handles);
        if (MB_SUCCESS == result) {
          result = mbImpl->tag_get_data(sharedps_tag, &(*rit), 1, tmp_procs);
          RRA("Trouble getting sharedps tag.");
          for (int j = 0; j < MAX_SHARING_PROCS; j++)
            if (tmp_procs[j] == to_proc) {
              to_vec[i] = tmp_handles[j];
              break;
            }
        }
      
        if (!to_vec[i]) {
          int j = new_ents.index(*rit);
          if (-1 == j) {
            result = MB_FAILURE;
            RRA("Failed to find new entity in send list.");
          }
          int err;
          to_vec[i] = CREATE_HANDLE(MBMAXTYPE, j, err);
          if (err) {
            result = MB_FAILURE;
            RRA("Failed to create handle in remote_handles.");
          }
        }
      }
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_remote_handles(const bool store_remote_handles,
                                               const MBRange &from_range, 
                                               MBRange &to_range,
                                               int to_proc,
                                               const MBRange &new_ents) 
{
  std::vector<MBEntityHandle> to_vector(from_range.size());

  MBErrorCode result =
    get_remote_handles(store_remote_handles, from_range, &to_vector[0],
                       to_proc, new_ents);
  RRA("Trouble getting remote handles.");
  std::copy(to_vector.begin(), to_vector.end(), mb_range_inserter(to_range));
  return result;
}

MBErrorCode MBParallelComm::unpack_entities(unsigned char *&buff_ptr,
                                            const bool store_remote_handles,
                                            const int from_ind,
                                            const bool is_iface,
                                            std::vector<std::vector<MBEntityHandle> > &L1hloc,
                                            std::vector<std::vector<MBEntityHandle> > &L1hrem,
                                            std::vector<std::vector<int> > &L1p,
                                            std::vector<MBEntityHandle> &L2hloc, 
                                            std::vector<MBEntityHandle> &L2hrem,
                                            std::vector<unsigned int> &L2p,
                                            MBRange &new_ents) 
{
    // general algorithm:
    // - unpack # entities
    // - save start of remote handle info, then scan forward to entity definition data
    // - for all vertices or entities w/ same # verts:
    //   . get entity type, num ents, and (if !vert) # verts 
    //   . for each ent:
    //      o get # procs/handles in remote handle info
    //      o if # procs/handles > 2, check for already-created entity:
    //        x get index of owner proc (1st in proc list), resize L1 list if nec
    //        x look for already-arrived entity in L2 by owner handle
    //      o if no existing entity:
    //        x if iface, look for existing entity with same connect & type
    //        x if none found, create vertex or element
    //        x if !iface & multi-shared, save on L2
    //        x if !iface, put new entity on new_ents list
    //      o update proc/handle, pstatus tags, adjusting to put owner first if iface
    //      o if !iface, save new handle on L1 for all sharing procs

    // lists of handles/procs to return to sending/other procs
    // L1hloc[p], L1hrem[p]: handle pairs [h, h'], where h is the local proc handle
    //         and h' is either the remote proc handle (if that is known) or
    //         the owner proc handle (otherwise);
    // L1p[p]: indicates whether h is remote handle (= -1) or owner (rank of owner)
    // L2hloc, L2hrem: local/remote handles for entities shared by > 2 procs;
    //         remote handles are on owning proc
    // L2p: owning procs for handles in L2hrem

  MBErrorCode result;
  bool done = false;
  MBReadUtilIface *ru = NULL;

  result = mbImpl->query_interface(std::string("MBReadUtilIface"), 
                                   reinterpret_cast<void**>(&ru));
  RRA("Failed to get MBReadUtilIface.");

    // procs the sending proc is telling me I'll be receiving from
  std::set<unsigned int> comm_procs;
  
    // 1. # entities = E
  int num_ents;
  unsigned char *buff_save = buff_ptr;
  int i, j;

  if (store_remote_handles) {
    UNPACK_INT(buff_ptr, num_ents);

    buff_save = buff_ptr;
    
      // save place where remote handle info starts, then scan forward to ents
    for (i = 0; i < num_ents; i++) {
      UNPACK_INT(buff_ptr, j);
      buff_ptr += j * (sizeof(int)+sizeof(MBEntityHandle));
    }
  }

  std::vector<MBEntityHandle> msg_ents;
  
  while (!done) {
    MBEntityType this_type = MBMAXTYPE;
    UNPACK_INT(buff_ptr, this_type);
    assert(this_type >= MBVERTEX && 
           (this_type == MBMAXTYPE || this_type < MBENTITYSET));

      // MBMAXTYPE signifies end of entities data
    if (MBMAXTYPE == this_type) break;

    assert(!is_iface || this_type != MBVERTEX);
    
      // get the number of ents
    int num_ents2, verts_per_entity;
    UNPACK_INT(buff_ptr, num_ents2);

      // unpack the nodes per entity
    if (MBVERTEX != this_type) {
      UNPACK_INT(buff_ptr, verts_per_entity);
    }
      
    std::vector<int> ps(MAX_SHARING_PROCS, -1);
    std::vector<MBEntityHandle> hs(MAX_SHARING_PROCS, 0);
    for (int e = 0; e < num_ents2; e++) {
        // check for existing entity, otherwise make new one
      MBEntityHandle new_h = 0;

      MBEntityHandle *connect;
      double *coords;
      int num_ps = -1;

        //=======================================
        // unpack all the data at once, to make sure the buffer pointers
        // are tracked correctly
        //=======================================
      if (store_remote_handles) {
          // pointers to other procs/handles
        UNPACK_INT(buff_save, num_ps);
        UNPACK_INTS(buff_save, &ps[0], num_ps);
        UNPACK_EH(buff_save, &hs[0], num_ps);
      }

      if (MBVERTEX == this_type) {
        coords = (double*) buff_ptr;
        buff_ptr += 3*sizeof(double);
      }
      else {
        connect = (MBEntityHandle*) buff_ptr;
        buff_ptr += verts_per_entity * sizeof(MBEntityHandle);

          // update connectivity to local handles
        result = get_local_handles(connect, verts_per_entity, msg_ents);
        RRA("Couldn't get local handles.");
      }

        //=======================================
        // now, process that data; begin by finding an identical 
        // entity, if there is one
        //=======================================
      if (store_remote_handles) {
        result = find_existing_entity(is_iface, ps[0], hs[0], num_ps, 
                                      connect, verts_per_entity,
                                      this_type,
                                      L2hloc, L2hrem, L2p,
                                      new_h);
        RRA("Trouble getting existing entity.");
      }

        //=======================================
        // if we didn't find one, we'll have to create one
        //=======================================

      bool created_here = false;
      if (!new_h) {
        
        if (MBVERTEX == this_type) {
            // create a vertex
          result = mbImpl->create_vertex(coords, new_h);
          RRA("Couldn't make new vertex.");
        }
        else {
            // create the element
          result = mbImpl->create_element(this_type, connect, verts_per_entity, new_h);
          RRA("Couldn't make new vertex.");

            // update adjacencies
          result = ru->update_adjacencies(new_h, 1, 
                                          verts_per_entity, connect);
          RRA("Failed to update adjacencies.");
        }

          // should have a new handle now
        assert(new_h);
        
          // if a new multi-shared entity, save owner for subsequent lookup in L2 lists
        if (store_remote_handles && !is_iface && num_ps > 2) {
          L2hrem.push_back(hs[0]);
          L2hloc.push_back(new_h);
          L2p.push_back(ps[0]);
        }

        created_here = true;
      }

        //=======================================
        // take care of sharing data
        //=======================================

        // need to save entities found in order, for interpretation of
        // later parts of this message
      if (!is_iface) msg_ents.push_back(new_h);

      if (created_here) new_ents.insert(new_h);

      if (store_remote_handles) {
        
          // update sharing data and pstatus, adjusting order if iface
        result = update_remote_data(new_h, &ps[0], &hs[0], num_ps, 
                                    (is_iface ? PSTATUS_INTERFACE :
                                     (created_here ? (PSTATUS_GHOST | PSTATUS_NOT_OWNED) : 0)));
        RRA("");
      
          // need to send this new handle to all sharing procs
        if (!is_iface) {
          for (j = 0; j < num_ps; j++) {
            if (ps[j] == (int)procConfig.proc_rank()) continue;
            int idx = get_buffers(ps[j]);
            if (idx == (int)L1hloc.size()) {
              L1hloc.resize(idx+1);
              L1hrem.resize(idx+1);
              L1p.resize(idx+1);
            }
            
              // don't bother adding if it's already in the list
            std::vector<MBEntityHandle>::iterator vit = 
                std::find(L1hloc[idx].begin(), L1hloc[idx].end(), new_h);
            if (vit != L1hloc[idx].end()) {
                // if it's in the list but remote handle isn't known but we know
                // it, replace in the list
              if (L1p[idx][vit-L1hloc[idx].begin()] != -1 && hs[j]) {
                L1hrem[idx][vit-L1hloc[idx].begin()] = hs[j];
                L1p[idx][vit-L1hloc[idx].begin()] = -1;
              }
              else continue;
            }
            else {
              if (!hs[j]) {
                assert(-1 != ps[0] && num_ps > 2);
                L1p[idx].push_back(ps[0]);
                L1hrem[idx].push_back(hs[0]);
              }
              else {
                assert("either this remote handle isn't in the remote list, or it's for another proc" &&
                       (std::find(L1hrem[idx].begin(), L1hrem[idx].end(), hs[j]) == 
                        L1hrem[idx].end() ||
                        L1p[idx][std::find(L1hrem[idx].begin(), L1hrem[idx].end(), hs[j]) - 
                                 L1hrem[idx].begin()] != -1));
                L1p[idx].push_back(-1);
                L1hrem[idx].push_back(hs[j]);
              }
              L1hloc[idx].push_back(new_h);
            }
          }
        }

        assert("Shouldn't be here for non-shared entities" &&
               -1 != num_ps);
        std::fill(&ps[0], &ps[num_ps], -1);
        std::fill(&hs[0], &hs[num_ps], 0);
      }
    }
    
    
#ifdef DEBUG_PACKING
      std::cerr << "Unpacked " << num_ents << " ents of type " 
                << MBCN::EntityTypeName(TYPE_FROM_HANDLE(this_type)) << std::endl;
#endif      

  }

#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done unpacking entities." << std::endl;
#endif

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::print_buffer(unsigned char *buff_ptr) 
{
    // 1. # entities = E
  int num_ents;
  int i, j, k;
  std::vector<int> ps;
  std::vector<MBEntityHandle> hs;

  UNPACK_INT(buff_ptr, num_ents);
  std::cout << num_ents << " entities..." << std::endl;

    // save place where remote handle info starts, then scan forward to ents
  for (i = 0; i < num_ents; i++) {
    UNPACK_INT(buff_ptr, j);
    ps.resize(j);
    hs.resize(j);
    std::cout << "Entity " << i << ": # procs = " << j << std::endl;
    UNPACK_INTS(buff_ptr, &ps[0], j);
    UNPACK_EH(buff_ptr, &hs[0], j);
    std::cout << "   Procs: ";
    for (k = 0; k < j; k++) std::cout << ps[k] << " ";
    std::cout << std::endl;
    std::cout << "   Handles: ";
    for (k = 0; k < j; k++) std::cout << hs[k] << " ";
    std::cout << std::endl;
  }
  
  while (true) {
    MBEntityType this_type = MBMAXTYPE;
    UNPACK_INT(buff_ptr, this_type);
    assert(this_type >= MBVERTEX && 
           (this_type == MBMAXTYPE || this_type < MBENTITYSET));

      // MBMAXTYPE signifies end of entities data
    if (MBMAXTYPE == this_type) break;

      // get the number of ents
    int num_ents2, verts_per_entity;
    UNPACK_INT(buff_ptr, num_ents2);

      // unpack the nodes per entity
    if (MBVERTEX != this_type) {
      UNPACK_INT(buff_ptr, verts_per_entity);
    }

    std::cout << "Type: " << MBCN::EntityTypeName(this_type)
              << "; num_ents = " << num_ents2;
    if (MBVERTEX != this_type) std::cout << "; verts_per_ent = " << verts_per_entity;
    std::cout << std::endl;
    
    for (int e = 0; e < num_ents2; e++) {
        // check for existing entity, otherwise make new one
      MBEntityHandle *connect;
      double *coords;

      if (MBVERTEX == this_type) {
        coords = (double*) buff_ptr;
        buff_ptr += 3*sizeof(double);
        std::cout << "xyz = " << coords[0] << ", " << coords[1] << ", " 
                  << coords[2] << std::endl;
      }
      else {
        connect = (MBEntityHandle*) buff_ptr;
        buff_ptr += verts_per_entity * sizeof(MBEntityHandle);

          // update connectivity to local handles
        std::cout << "Connectivity: ";
        for (k = 0; k < verts_per_entity; k++) std::cout << connect[k] << " ";
        std::cout << std::endl;
      }
    }
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::list_entities(const MBEntityHandle *ents, int num_ents) 
{
  if (NULL == ents && 0 == num_ents) {
    sharedEnts.print("Shared entities:\n");
    return MB_SUCCESS;
  }
  
  else if (NULL == ents || 0 == num_ents) {
    return list_entities(sharedEnts);
  }
    
  unsigned char pstat;
  MBEntityHandle tmp_handles[MAX_SHARING_PROCS];
  int tmp_procs[MAX_SHARING_PROCS];
  unsigned int num_ps;
  MBErrorCode result;

  for (int i = 0; i < num_ents; i++) {
    result = get_sharing_data(ents[i], tmp_procs, tmp_handles, pstat, num_ps);
    RRA("Failed to get sharing data.");

    result = mbImpl->list_entities(ents+i, 1);
    std::cout << "Pstatus: ";
    if (!num_ps)
      std::cout << "local " << std::endl;
    else {
      if (pstat & PSTATUS_NOT_OWNED) std::cout << "NOT_OWNED; ";
      if (pstat & PSTATUS_SHARED) std::cout << "SHARED; ";
      if (pstat & PSTATUS_MULTISHARED) std::cout << "MULTISHARED; ";
      if (pstat & PSTATUS_INTERFACE) std::cout << "INTERFACE; ";
      if (pstat & PSTATUS_GHOST) std::cout << "GHOST; ";
      std::cout << std::endl;
      for (unsigned int j = 0; j < num_ps; j++) {
        std::cout << "  proc " << tmp_procs[j] << " id (handle) " 
                  << mbImpl->id_from_handle(tmp_handles[j]) 
                  << "(" << tmp_handles[j] << ")" << std::endl;
      }
    }
    std::cout << std::endl;
  }

  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::list_entities(const MBRange &ents) 
{
  for (MBRange::iterator rit = ents.begin(); rit != ents.end(); rit++)
    list_entities(&(*rit), 1);
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::update_remote_data(MBRange &local_range,
                                               MBRange &remote_range,
                                               int other_proc,
                                               const unsigned char add_pstat) 
{
  MBRange::iterator rit, rit2;
  MBErrorCode result = MB_SUCCESS;

    // for each pair of local/remote handles:
  for (rit = local_range.begin(), rit2 = remote_range.begin(); 
       rit != local_range.end(); rit++, rit2++) {

    result = update_remote_data(*rit, &other_proc, &(*rit2), 1, add_pstat);
    RRA(" ");
  }

  return result;
}
  
MBErrorCode MBParallelComm::update_remote_data(const MBEntityHandle new_h,
                                               const int *ps,
                                               const MBEntityHandle *hs,
                                               const int num_ps,
                                               const unsigned char add_pstat) 
{
  MBEntityHandle tag_hs[MAX_SHARING_PROCS];
  int tag_ps[MAX_SHARING_PROCS];
  unsigned char pstat;
    // get initial sharing data; tag_ps and tag_hs get terminated with -1 and 0
    // in this function, so no need to initialize
  unsigned int num_exist;
  MBErrorCode result = get_sharing_data(new_h, tag_ps, tag_hs, pstat, num_exist);
  RRA("");
  
    // add any new sharing data
  bool changed = false;
  int idx;
  if (!num_exist) {
      // just take what caller passed
    memcpy(tag_ps, ps, num_ps*sizeof(int));
    memcpy(tag_hs, hs, num_ps*sizeof(MBEntityHandle));
    num_exist = num_ps;
      // if it's only one, hopefully I'm not there yet...
    assert("I shouldn't be the only proc there." &&
           (1 != num_exist || ps[0] != (int)procConfig.proc_rank()));
    changed = true;
  }
  else {
    for (int i = 0; i < num_ps; i++) {
      idx = std::find(tag_ps, tag_ps+num_exist, ps[i]) - tag_ps;
      if (idx == (int) num_exist) {
        tag_ps[num_exist] = ps[i];
        tag_hs[num_exist] = hs[i];
        num_exist++;
        changed = true;
      }
      else if (0 == tag_hs[idx]) {
        tag_hs[idx] = hs[i];
        changed = true;
      }
      else if (0 != hs[i]) {
        assert(hs[i] == tag_hs[idx]);
      }
    }
  }
  
    // adjust for interface layer if necessary
  if (add_pstat & PSTATUS_INTERFACE) {
    idx = std::min_element(tag_ps, tag_ps+num_exist) - tag_ps;
    if (idx) {
      int tag_proc = tag_ps[idx];
      tag_ps[idx] = tag_ps[0];
      tag_ps[0] = tag_proc;
      MBEntityHandle tag_h = tag_hs[idx];
      tag_hs[idx] = tag_hs[0];
      tag_hs[0] = tag_h;
      changed = true;
      if (tag_ps[0] != (int)procConfig.proc_rank()) pstat |= PSTATUS_NOT_OWNED;
    }
  }
    
  if (!changed) return MB_SUCCESS;
  
  assert("interface entities should have > 1 proc" &&
         (!(add_pstat & PSTATUS_INTERFACE) || num_exist > 1));
  assert("ghost entities should have > 1 proc" &&
         (!(add_pstat & PSTATUS_GHOST) || num_exist > 1));
  
    // if it's multi-shared and we created the entity in this unpack,
    // local handle probably isn't in handle list yet
  if (num_exist > 2) {
    idx = std::find(tag_ps, tag_ps+num_exist, procConfig.proc_rank()) - tag_ps;
    assert(idx < (int) num_exist);
    if (!tag_hs[idx])
      tag_hs[idx] = new_h;
  }
      
  int tag_p;
  MBEntityHandle tag_h;

    // reset single shared proc/handle if was shared and moving to multi-shared
  if (num_exist > 2 && !(pstat & PSTATUS_MULTISHARED) &&
      (pstat & PSTATUS_SHARED)) {
      // must remove sharedp/h first, which really means set to default value
    tag_p = -1;
    result = mbImpl->tag_set_data(sharedp_tag(), &new_h, 1, &tag_p);
    RRA("Couldn't set sharedp tag.");
    tag_h = 0;
    result = mbImpl->tag_set_data(sharedh_tag(), &new_h, 1, &tag_h);
    RRA("Couldn't set sharedh tag.");
  }

    // update pstat
  pstat |= add_pstat;
  
    // set sharing tags
  if (num_exist > 2) {
    std::fill(tag_ps+num_exist, tag_ps+MAX_SHARING_PROCS, -1);
    std::fill(tag_hs+num_exist, tag_hs+MAX_SHARING_PROCS, 0);
    result = mbImpl->tag_set_data(sharedps_tag(), &new_h, 1, tag_ps);
    RRA("Couldn't set sharedps tag.");
    result = mbImpl->tag_set_data(sharedhs_tag(), &new_h, 1, tag_hs);
    RRA("Couldn't set sharedhs tag.");
    pstat |= (PSTATUS_MULTISHARED | PSTATUS_SHARED);
  }
  else if (num_exist == 2 || num_exist == 1) {
    if (tag_ps[0] == (int) procConfig.proc_rank()) {
      assert(2 == num_exist);
      tag_ps[0] = tag_ps[1];
      tag_hs[0] = tag_hs[1];
    }
    assert(tag_ps[0] != -1 && tag_hs[0] != 0);
    result = mbImpl->tag_set_data(sharedp_tag(), &new_h, 1, tag_ps);
    RRA("Couldn't set sharedp tag.");
    result = mbImpl->tag_set_data(sharedh_tag(), &new_h, 1, tag_hs);
    RRA("Couldn't set sharedh tag.");
    pstat |= PSTATUS_SHARED;
  }

    // now set new pstatus
  result = mbImpl->tag_set_data(pstatus_tag(), &new_h, 1, &pstat);
  RRA("Couldn't set pstatus tag.");

  if (pstat & PSTATUS_SHARED) sharedEnts.insert(new_h);
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_sharing_data(MBEntityHandle entity,
                                             int *ps, 
                                             MBEntityHandle *hs,
                                             unsigned char &pstat,
                                             unsigned int &num_ps)
{
  MBErrorCode result = mbImpl->tag_get_data(pstatus_tag(), &entity, 1, &pstat);
  RRA("Couldn't get pstatus tag.");
  if (pstat & PSTATUS_MULTISHARED) {
    result = mbImpl->tag_get_data(sharedps_tag(), &entity, 1, ps);
    RRA("Couldn't get sharedps tag.");
    if (hs) {
      result = mbImpl->tag_get_data(sharedhs_tag(), &entity, 1, hs);
      RRA("Couldn't get sharedhs tag.");
    }
    num_ps = std::find(ps, ps+MAX_SHARING_PROCS, -1) - ps;
  }
  else if (pstat & PSTATUS_SHARED) {
    result = mbImpl->tag_get_data(sharedp_tag(), &entity, 1, ps);
    RRA("Couldn't get sharedp tag.");
    if (hs) {
      result = mbImpl->tag_get_data(sharedh_tag(), &entity, 1, hs);
      RRA("Couldn't get sharedh tag.");
      hs[1] = 0;
    }
      // initialize past end of data
    ps[1] = -1;
    num_ps = 1;
  }
  else {
    ps[0] = -1;
    if (hs) hs[0] = 0;
    num_ps = 0;
  }

  assert(0 <= num_ps && MAX_SHARING_PROCS >= num_ps);
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::find_existing_entity(const bool is_iface,
                                                 const int owner_p,
                                                 const MBEntityHandle owner_h,
                                                 const int num_ps,
                                                 const MBEntityHandle *connect,
                                                 const int num_connect,
                                                 const MBEntityType this_type,
                                                 std::vector<MBEntityHandle> &L2hloc,
                                                 std::vector<MBEntityHandle> &L2hrem,
                                                 std::vector<unsigned int> &L2p,
                                                 MBEntityHandle &new_h) 
{
  new_h = 0;
  if (!is_iface && num_ps > 2) {
    for (unsigned int i = 0; i < L2hrem.size(); i++) {
      if (L2hrem[i] == owner_h && owner_p == (int) L2p[i]) {
        new_h = L2hloc[i];
        return MB_SUCCESS;
      }
    }        
  }

    // if we got here and it's a vertex, we don't need to look further
  if (MBVERTEX == this_type || !connect || !num_connect) return MB_SUCCESS;
  
  MBRange tmp_range;
  MBErrorCode result = mbImpl->get_adjacencies(connect, num_connect, 
                                               MBCN::Dimension(this_type), false, 
                                               tmp_range);
  RRA("Problem getting existing entity.");
  if (!tmp_range.empty()) {
      // found a corresponding entity - return target
    new_h = *tmp_range.begin();
  }  
  else {
    new_h = 0;
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_local_handles(const MBRange &remote_handles,
                                              MBRange &local_handles,
                                              const MBRange &new_ents) 
{
  std::vector<MBEntityHandle> rh_vec;
  rh_vec.reserve(remote_handles.size());
  std::copy(remote_handles.begin(), remote_handles.end(), std::back_inserter(rh_vec));
  MBErrorCode result = get_local_handles(&rh_vec[0], remote_handles.size(), new_ents);
  std::copy(rh_vec.begin(), rh_vec.end(), mb_range_inserter(local_handles));
  return result;
}
  
MBErrorCode MBParallelComm::get_local_handles(MBEntityHandle *from_vec, 
                                              int num_ents,
                                              const MBRange &new_ents) 
{
  std::vector<MBEntityHandle> tmp_ents;
  std::copy(new_ents.begin(), new_ents.end(), std::back_inserter(tmp_ents));
  return get_local_handles(from_vec, num_ents, tmp_ents);
}

MBErrorCode MBParallelComm::get_local_handles(MBEntityHandle *from_vec,
                                              int num_ents,
                                              const std::vector<MBEntityHandle> &new_ents) 
{
  for (int i = 0; i < num_ents; i++) {
    if (TYPE_FROM_HANDLE(from_vec[i]) == MBMAXTYPE) {
      assert(ID_FROM_HANDLE(from_vec[i]) < (int) new_ents.size());
      from_vec[i] = new_ents[ID_FROM_HANDLE(from_vec[i])];
    }
  }
  
  return MB_SUCCESS;
}

template <typename T> void
insert_in_array( T* array, size_t array_size, size_t location, T value )
{
  assert( location+1 < array_size );
  for (size_t i = array_size-1; i > location; --i)
    array[i] = array[i-1];
  array[location] = value;
}

MBErrorCode MBParallelComm::pack_range_map(MBRange &key_range, MBEntityHandle val_start,
                                           HandleMap &handle_map) 
{
  for (MBRange::const_pair_iterator key_it = key_range.const_pair_begin(); 
       key_it != key_range.const_pair_end(); key_it++) {
    int tmp_num = (*key_it).second - (*key_it).first + 1;
    handle_map.insert((*key_it).first, val_start, tmp_num);
    val_start += tmp_num;
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_sets(MBRange &entities,
                                      std::vector<unsigned char> &buff,
                                      unsigned char *&buff_ptr,
                                      const bool store_remote_handles,
                                      const int to_proc)
{
    // SETS:
    // . #sets
    // . for each set:
    //   - options[#sets] (unsigned int)
    //   - if (unordered) set range 
    //   - else if ordered
    //     . #ents in set
    //     . handles[#ents]
    //   - #parents
    //   - if (#parents) handles[#parents]
    //   - #children
    //   - if (#children) handles[#children]
  
    // now the sets; assume any sets the application wants to pass are in the entities list
  MBErrorCode result;
  MBRange all_sets = entities.subset_by_type(MBENTITYSET);

  int buff_size = estimate_sets_buffer_size(all_sets, store_remote_handles);
  CHECK_BUFF_SPACE(buff, buff_ptr, buff_size);

    // number of sets
  PACK_INT(buff_ptr, all_sets.size());

    // options for all sets
  std::vector<unsigned int> options(all_sets.size());
  MBRange::iterator rit;
  std::vector<MBEntityHandle> members;
  int i;
  for (rit = all_sets.begin(), i = 0; rit != all_sets.end(); rit++, i++) {
      result = mbImpl->get_meshset_options(*rit, options[i]);
      RRA("Failed to get meshset options.");
  }
  CHECK_BUFF_SPACE(buff, buff_ptr, all_sets.size()*sizeof(unsigned int));
  PACK_VOID(buff_ptr, &options[0], all_sets.size()*sizeof(unsigned int));
  
    // vectors/ranges
  for (rit = all_sets.begin(), i = 0; rit != all_sets.end(); rit++, i++) {
      MBRange set_range;
      if (options[i] & MESHSET_SET) {
        MBRange set_range;
        result = mbImpl->get_entities_by_handle(*rit, set_range);
        RRA("Failed to get set entities.");

        buff_size = RANGE_SIZE(set_range);
        CHECK_BUFF_SPACE(buff, buff_ptr, buff_size);
        PACK_RANGE(buff_ptr, set_range);
      }
      else if (options[i] & MESHSET_ORDERED) {
        members.clear();
        result = mbImpl->get_entities_by_handle(*rit, members);
        RRA("Failed to get entities in ordered set.");
        
        CHECK_BUFF_SPACE(buff, buff_ptr,
                         members.size()*sizeof(MBEntityHandle)+sizeof(int));
        PACK_INT(buff_ptr, members.size());
        PACK_EH(buff_ptr, &members[0], members.size());
      }
  }
    // pack numbers of parents/children
  unsigned int tot_pch = 0;
  int num_pch;
  CHECK_BUFF_SPACE(buff, buff_ptr, 2*all_sets.size()*sizeof(int));
  for (rit = all_sets.begin(), i = 0; rit != all_sets.end(); rit++, i++) {
      // pack parents
    result = mbImpl->num_parent_meshsets(*rit, &num_pch);
    RRA("Failed to get num parents.");
    PACK_INT(buff_ptr, num_pch);
    tot_pch += num_pch;
    result = mbImpl->num_child_meshsets(*rit, &num_pch);
    RRA("Failed to get num children.");
    PACK_INT(buff_ptr, num_pch);
    tot_pch += num_pch;
  }

    // now pack actual parents/children
  members.clear();
  members.reserve(tot_pch);
  std::vector<MBEntityHandle> tmp_pch;
  for (rit = all_sets.begin(), i = 0; rit != all_sets.end(); rit++, i++) {
    result = mbImpl->get_parent_meshsets(*rit, tmp_pch);
    RRA("Failed to get parents.");
    std::copy(tmp_pch.begin(), tmp_pch.end(), std::back_inserter(members));
    tmp_pch.clear();
    result = mbImpl->get_child_meshsets(*rit, tmp_pch);
    RRA("Failed to get children.");
    std::copy(tmp_pch.begin(), tmp_pch.end(), std::back_inserter(members));
    tmp_pch.clear();
  }
  assert(members.size() == tot_pch);
  if (!members.empty()) {
    result = get_remote_handles(store_remote_handles,
                                &members[0], &members[0], 
                                members.size(), to_proc,
                                entities);
    RRA("Trouble getting remote handles for set parent/child sets.");
#ifndef NDEBUG
      // check that all handles are either sets or maxtype
    for (unsigned int __j = 0; __j < members.size(); __j++)
      assert((TYPE_FROM_HANDLE(members[__j]) == MBMAXTYPE &&
              ID_FROM_HANDLE(members[__j]) < (int)entities.size()) ||
             TYPE_FROM_HANDLE(members[__j]) == MBENTITYSET);
#endif        
    CHECK_BUFF_SPACE(buff, buff_ptr, members.size()*sizeof(MBEntityHandle));
    PACK_EH(buff_ptr, &members[0], members.size());
  }
    
    // pack the handles
  if (store_remote_handles && !all_sets.empty()) {
    buff_size = RANGE_SIZE(all_sets);
    CHECK_BUFF_SPACE(buff, buff_ptr, buff_size);
    PACK_RANGE(buff_ptr, all_sets);
  }
  
#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done packing sets." << std::endl;
#endif

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::unpack_sets(unsigned char *&buff_ptr,
                                        MBRange &entities,
                                        const bool store_remote_handles,
                                        const int from_proc)
{
  
    // now the sets; assume any sets the application wants to pass are in the entities list
  MBErrorCode result;

  MBRange new_sets;
  int num_sets;
  UNPACK_INT(buff_ptr, num_sets);

  if (!num_sets) return MB_SUCCESS;
         
  std::vector<MBEntityHandle> members;
  int num_ents;
  std::vector<unsigned int> options_vec(num_sets);
      // option value
  if (num_sets)
    UNPACK_VOID(buff_ptr, &options_vec[0], num_sets*sizeof(unsigned int));

    // create sets
  int i;
  MBRange::const_iterator rit;
  for (i = 0; i < num_sets; i++) {
    
      // create the set
    MBEntityHandle set_handle;
    result = mbImpl->create_meshset(options_vec[i], set_handle);
    RRA("Failed to create set in unpack.");

    // make sure new sets handles are monotonically increasing
    assert(set_handle > *new_sets.rbegin());

    new_sets.insert(set_handle);
  }

  entities.merge(new_sets);
  
  for (rit = new_sets.begin(), i = 0; rit != new_sets.end(); rit++, i++) {
    if (options_vec[i] & MESHSET_SET) {
        // unpack entities as a range
      MBRange set_range, tmp_range;
      UNPACK_RANGE(buff_ptr, tmp_range);
      result = get_local_handles(tmp_range, set_range, entities);      
      RRA("Failed to get local handles for unordered set contents.");
      result = mbImpl->add_entities(*rit, set_range);
      RRA("Failed to add ents to unordered set in unpack.");
    }
    else if (options_vec[i] & MESHSET_ORDERED) {
        // unpack entities as vector, with length
      UNPACK_INT(buff_ptr, num_ents);
      members.resize(num_ents);
      if (num_ents) UNPACK_EH(buff_ptr, &members[0], num_ents);
      result = get_local_handles(&members[0], num_ents, entities);
      RRA("Failed to get local handles for ordered set contents.");
      result = mbImpl->add_entities(*rit, &members[0], num_ents);
      RRA("Failed to add ents to ordered set in unpack.");
    }
  }

  std::vector<int> num_pch(2*new_sets.size());
  std::vector<int>::iterator vit;
  int tot_pch = 0;
  for (vit = num_pch.begin(); vit != num_pch.end(); vit++) {
    UNPACK_INT(buff_ptr, *vit);
    tot_pch += *vit;
  }
  
  members.resize(tot_pch);
  UNPACK_EH(buff_ptr, &members[0], tot_pch);
  result = get_local_handles(&members[0], tot_pch, entities);
  RRA("Couldn't get local handle for parent/child sets.");

  int num = 0;
  MBEntityHandle *mem_ptr = &members[0];
  for (rit = new_sets.begin(); rit != new_sets.end(); rit++) {
      // unpack parents/children
    int num_par = num_pch[num++], num_child = num_pch[num++];
    if (num_par+num_child) {
      for (i = 0; i < num_par; i++) {
        assert(0 != mem_ptr[i]);
        result = mbImpl->add_parent_meshset(*rit, mem_ptr[i]);
        RRA("Failed to add parent to set in unpack.");
      }
      mem_ptr += num_par;
      for (i = 0; i < num_child; i++) {
        assert(0 != mem_ptr[i]);
        result = mbImpl->add_child_meshset(*rit, mem_ptr[i]);
        RRA("Failed to add child to set in unpack.");
      }
      mem_ptr += num_child;
    }
  }

    // unpack source handles
  MBRange dum_range;
  if (store_remote_handles && !new_sets.empty()) {
    UNPACK_RANGE(buff_ptr, dum_range);
    result = update_remote_data(new_sets, dum_range, from_proc, 0);
    RRA("Couldn't set sharing data for sets");
  }

#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done unpacking sets." << std::endl;
#endif

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_adjacencies(MBRange &entities,
                                             MBRange::const_iterator &start_rit,
                                             MBRange &whole_range,
                                             unsigned char *&buff_ptr,
                                             int &count,
                                             const bool just_count,
                                             const bool store_handles,
                                             const int to_proc)
{
  return MB_FAILURE;
}

MBErrorCode MBParallelComm::unpack_adjacencies(unsigned char *&buff_ptr,
                                               MBRange &entities,
                                               const bool store_handles,
                                               const int from_proc)
{
  return MB_FAILURE;
}

MBErrorCode MBParallelComm::pack_tags(MBRange &entities,
                                      const std::vector<MBTag> &src_tags,
                                      const std::vector<MBTag> &dst_tags,
                                      const std::vector<MBRange> &tag_ranges,
                                      std::vector<unsigned char> &buff,
                                      unsigned char *&buff_ptr,
                                      const bool store_remote_handles,
                                      const int to_proc)
{
  

  MBErrorCode result;
  std::vector<MBTag>::const_iterator tag_it, dst_it;
  std::vector<MBRange>::const_iterator rit;
  int count = 0;
  
  for (tag_it = src_tags.begin(), rit = tag_ranges.begin(); 
       tag_it != src_tags.end(); tag_it++, rit++) {

    result = packed_tag_size( *tag_it, *rit, count );
    if (MB_SUCCESS != result)
      return result;
  }
    
    // number of tags
  count += sizeof(int);

  CHECK_BUFF_SPACE(buff, buff_ptr, count);
  
  PACK_INT(buff_ptr, src_tags.size());
    
  for (tag_it = src_tags.begin(), dst_it = dst_tags.begin(), rit = tag_ranges.begin(); 
       tag_it != src_tags.end(); tag_it++, dst_it++, rit++) {
    
    result = pack_tag( *tag_it, *dst_it, *rit, entities, buff, buff_ptr, 
                       store_remote_handles, to_proc );
    if (MB_SUCCESS != result)
      return result;
  }
  
#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done packing tags." << std::endl;
#endif

  return MB_SUCCESS;
}
         

MBErrorCode MBParallelComm::packed_tag_size( MBTag tag,
                                             const MBRange &tagged_entities,
                                             int &count )
{
    // for dense tags, compute size assuming all entities have that tag
    // for sparse tags, get number of entities w/ that tag to compute size

  std::vector<int> var_len_sizes;
  std::vector<const void*> var_len_values;
    
  const TagInfo *tinfo = tagServer->get_tag_info(tag);
    // default value
  count += sizeof(int);
  if (NULL != tinfo->default_value()) 
    count += tinfo->default_value_size();

    // size, type, data type
  count += 3*sizeof(int);

    // name
  count += sizeof(int);
  count += tinfo->get_name().size();

    // range of tag
  count += sizeof(int) + tagged_entities.size() * sizeof(MBEntityHandle);

  if (tinfo->get_size() == MB_VARIABLE_LENGTH) {
    const int num_ent = tagged_entities.size();
      // send a tag size for each entity
    count += num_ent * sizeof(int);
      // send tag data for each entity
    var_len_sizes.resize( num_ent );
    var_len_values.resize( num_ent );
    MBErrorCode result = tagServer->get_data( tag,
                                              tagged_entities, 
                                              &var_len_values[0], 
                                              &var_len_sizes[0] );
    RRA("Failed to get lenghts of variable-length tag values.");
    count += std::accumulate( var_len_sizes.begin(), var_len_sizes.end(), 0 );
  }
  else {
      // tag data values for range or vector
    count += tagged_entities.size() * tinfo->get_size();
  }
  
  return MB_SUCCESS;
}


MBErrorCode MBParallelComm::pack_tag( MBTag src_tag,
                                      MBTag dst_tag,
                                      const MBRange &tagged_entities,
                                      const MBRange &whole_range,
                                      std::vector<unsigned char> &buff,
                                      unsigned char *&buff_ptr,
                                      const bool store_remote_handles,
                                      const int to_proc )
{
  MBErrorCode result;
  std::vector<int> var_len_sizes;
  std::vector<const void*> var_len_values;

  const TagInfo* tinfo = tagServer->get_tag_info(src_tag);
  if (!tinfo)
    return MB_TAG_NOT_FOUND;
    
  const TagInfo* dst_tinfo;
  if (src_tag == dst_tag) {
    dst_tinfo = tinfo;
  }
  else {
    dst_tinfo = tagServer->get_tag_info(dst_tag);
    if (!dst_tinfo)
      return MB_TAG_NOT_FOUND;
    if (dst_tinfo->get_size() != tinfo->get_size())
      return MB_TYPE_OUT_OF_RANGE;
    if (dst_tinfo->get_data_type() != tinfo->get_data_type() && 
        dst_tinfo->get_data_type() != MB_TYPE_OPAQUE &&
            tinfo->get_data_type() != MB_TYPE_OPAQUE)
      return MB_TYPE_OUT_OF_RANGE;
  }
    
    // size, type, data type
  CHECK_BUFF_SPACE(buff, buff_ptr, 3*sizeof(int));
  PACK_INT(buff_ptr, tinfo->get_size());
  MBTagType this_type;
  result = mbImpl->tag_get_type(dst_tag, this_type);
  PACK_INT(buff_ptr, (int)this_type);
  PACK_INT(buff_ptr, (int)(tinfo->get_data_type()));

    // default value
  if (NULL == tinfo->default_value()) {
    CHECK_BUFF_SPACE(buff, buff_ptr, sizeof(int));
    PACK_INT(buff_ptr, 0);
  }
  else {
    CHECK_BUFF_SPACE(buff, buff_ptr, tinfo->default_value_size());
    PACK_BYTES(buff_ptr, tinfo->default_value(), tinfo->default_value_size());
  }

    // name
  CHECK_BUFF_SPACE(buff, buff_ptr, tinfo->get_name().size());
  PACK_BYTES(buff_ptr, dst_tinfo->get_name().c_str(), dst_tinfo->get_name().size());

#ifdef DEBUG_PACKING
  std::cerr << "Packing tag \"" << tinfo->get_name() << "\"";
  if (tinfo != dst_tinfo)
    std::cerr << " (as tag \"" << dst_tinfo->get_name() << "\")";
  std::cerr << std::endl;
#endif    
    // pack entities
  CHECK_BUFF_SPACE(buff, buff_ptr, tagged_entities.size()*sizeof(MBEntityHandle)+sizeof(int));
  PACK_INT(buff_ptr, tagged_entities.size());
  result = get_remote_handles(store_remote_handles,
                              tagged_entities, (MBEntityHandle*)buff_ptr, to_proc,
                              whole_range);
#ifdef DEBUG_PACKING
  if (MB_SUCCESS != result) {
    std::cerr << "Trouble getting remote handles for tagged entities:" << std::endl;
    tagged_entities.print("  ");
  }
#else
  RRA("Trouble getting remote handles for tagged entities.");
#endif

  buff_ptr += tagged_entities.size() * sizeof(MBEntityHandle);

  const size_t num_ent = tagged_entities.size();
  if (tinfo->get_size() == MB_VARIABLE_LENGTH) {
    var_len_sizes.resize( num_ent, 0 );
    var_len_values.resize( num_ent, 0 );
    result = mbImpl->tag_get_data(src_tag, tagged_entities, &var_len_values[0], 
                                  &var_len_sizes[0] );
    RRA("Failed to get variable-length tag data in pack_tags.");
    CHECK_BUFF_SPACE(buff, buff_ptr, num_ent*sizeof(int));
    PACK_INTS(buff_ptr, &var_len_sizes[0], num_ent);
    for (unsigned int i = 0; i < num_ent; ++i) {
      CHECK_BUFF_SPACE(buff, buff_ptr, var_len_sizes[i]);
      PACK_VOID(buff_ptr, var_len_values[i], var_len_sizes[i]);
    }
  }
  else {
    CHECK_BUFF_SPACE(buff, buff_ptr, num_ent * tinfo->get_size());
    result = mbImpl->tag_get_data(src_tag, tagged_entities, buff_ptr);
    RRA("Failed to get tag data in pack_tags.");
    buff_ptr += num_ent * tinfo->get_size();
    PC(num_ent*tinfo->get_size(), " void");
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_tag_send_list( const MBRange& whole_range,
                                               std::vector<MBTag>& all_tags,
                                               std::vector<MBRange>& tag_ranges )
{
  std::vector<MBTag> tmp_tags;
  MBErrorCode result = tagServer->get_tags(tmp_tags);
  RRA("Failed to get tags in pack_tags.");

  std::vector<MBTag>::iterator tag_it;
  for (tag_it = tmp_tags.begin(); tag_it != tmp_tags.end(); tag_it++) {
    std::string tag_name;
    result = mbImpl->tag_get_name(*tag_it, tag_name);
    if (tag_name.c_str()[0] == '_' && tag_name.c_str()[1] == '_')
      continue;

    MBRange tmp_range;
    result = tagServer->get_entities(*tag_it, tmp_range);
    RRA("Failed to get entities for tag in pack_tags.");
    tmp_range = tmp_range.intersect(whole_range);

    if (tmp_range.empty()) continue;
        
      // ok, we'll be sending this tag
    all_tags.push_back( *tag_it );
    tag_ranges.push_back( MBRange() );
    tag_ranges.back().swap( tmp_range );
  }
  
  return MB_SUCCESS;
}



MBErrorCode MBParallelComm::unpack_tags(unsigned char *&buff_ptr,
                                        MBRange &entities,
                                        const bool store_remote_handles,
                                        const int from_proc)
{
    // tags
    // get all the tags
    // for dense tags, compute size assuming all entities have that tag
    // for sparse tags, get number of entities w/ that tag to compute size

  MBErrorCode result;
  
  int num_tags;
  UNPACK_INT(buff_ptr, num_tags);
  std::vector<MBEntityHandle> tag_ents;
  std::vector<const void*> var_len_vals;
  std::vector<int> var_lengths;

  for (int i = 0; i < num_tags; i++) {
    
        // tag handle
    MBTag tag_handle;

      // size, data type
    int tag_size, tag_data_type, tag_type;
    UNPACK_INT(buff_ptr, tag_size);
    UNPACK_INT(buff_ptr, tag_type);
    UNPACK_INT(buff_ptr, tag_data_type);
      
      // default value
    int def_val_size;
    UNPACK_INT(buff_ptr, def_val_size);
    void *def_val_ptr = NULL;
    if (def_val_size) {
      def_val_ptr = buff_ptr;
      buff_ptr += def_val_size;
      UPC(tag_size, " void");
    }
    
      // name
    int name_len;
    UNPACK_INT(buff_ptr, name_len);
    std::string tag_name( reinterpret_cast<char*>(buff_ptr), name_len );
    buff_ptr += name_len;
    UPC(64, " chars");
#ifdef DEBUG_PACKING
    std::cerr << "Unpacking tag " << tag_name << std::endl;
#endif    

      // create the tag
    if (tag_size == MB_VARIABLE_LENGTH) 
      result = mbImpl->tag_create_variable_length( tag_name.c_str(), (MBTagType)tag_type,
                                                   (MBDataType)tag_data_type, tag_handle,
                                                   def_val_ptr, def_val_size );
    else
      result = mbImpl->tag_create(tag_name.c_str(), tag_size, (MBTagType) tag_type, 
                                  (MBDataType) tag_data_type, tag_handle,
                                  def_val_ptr);
    if (MB_ALREADY_ALLOCATED == result) {
        // already allocated tag, check to make sure it's the same size, type, etc.
      const TagInfo *tag_info = tagServer->get_tag_info(tag_name.c_str());
      MBTagType this_type;
      result = mbImpl->tag_get_type(tag_handle, this_type);
      if (tag_size != tag_info->get_size() ||
          tag_type != this_type ||
          tag_data_type != tag_info->get_data_type() ||
          (def_val_ptr && !tag_info->default_value()) ||
          (!def_val_ptr && tag_info->default_value())) {
        RRA("Didn't get correct tag info when unpacking tag.");
      }
    }
    else if (MB_SUCCESS != result) return result;

      // go through handle vec (in buffer) and convert to local handles in-place
    int num_ents;
    UNPACK_INT(buff_ptr, num_ents);
    MBEntityHandle *handle_vec = (MBEntityHandle*)buff_ptr;
    result = get_local_handles(handle_vec, num_ents, entities);
    RRA("Failed to get local handles for tagged entities.");
    buff_ptr += num_ents * sizeof(MBEntityHandle);

      // if it's a handle type, also convert tag vals in-place in buffer
    if (MB_TYPE_HANDLE == tag_type) {
      MBEntityHandle *val_vec = (MBEntityHandle*)buff_ptr;
      result = get_local_handles(val_vec, num_ents, entities);
      RRA("Failed to get local handles for tag vals.");
    }

    if (tag_size == MB_VARIABLE_LENGTH) {
        // Be careful of alignment here.  If the integers are aligned
        // in the buffer, we can use them directly.  Otherwise we must
        // copy them.
      const int* size_arr;
      if (((size_t)buff_ptr)%4) {
        var_lengths.resize( num_ents );
        memcpy( &var_lengths[0], buff_ptr, num_ents*sizeof(int) );
        size_arr = &var_lengths[0];
      }
      else {
        size_arr = reinterpret_cast<const int*>(buff_ptr);
      }
      buff_ptr += sizeof(int) * num_ents;
      UPC(sizeof(int) * num_ents, " void");
      
        // get pointers into buffer for each tag value
      var_len_vals.resize(num_ents);
      for (std::vector<MBEntityHandle>::size_type i = 0; 
           i < (std::vector<MBEntityHandle>::size_type) num_ents; ++i) {
        var_len_vals[i] = buff_ptr;
        buff_ptr += size_arr[i];
        UPC(size_arr[i], " void");
      }
      result = mbImpl->tag_set_data( tag_handle, handle_vec, num_ents,
                                     &var_len_vals[0], size_arr );
      RRA("Trouble setting tag data when unpacking variable-length tag.");
    }
    else {
      result = mbImpl->tag_set_data(tag_handle, handle_vec,
                                    num_ents, buff_ptr);
      RRA("Trouble setting range-based tag data when unpacking tag.");
      buff_ptr += num_ents * tag_size;
      UPC(num_ents * tag_size, " void");
    }
  }
  
#ifdef DEBUG_PACKING
  std::cerr << std::endl << "Done unpacking tags." << std::endl;
#endif

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::resolve_shared_ents(MBEntityHandle this_set,
                                                int resolve_dim,
                                                int shared_dim) 
{
  MBErrorCode result;
  MBRange proc_ents;
      // get the entities in the partition sets
  for (MBRange::iterator rit = partitionSets.begin(); rit != partitionSets.end(); rit++) {
    MBRange tmp_ents;
    result = mbImpl->get_entities_by_handle(*rit, tmp_ents, true);
    if (MB_SUCCESS != result) return result;
    proc_ents.merge(tmp_ents);
  }

    // resolve dim is maximal dim of entities in proc_ents
  if (-1 == resolve_dim) {
    resolve_dim = mbImpl->dimension_from_handle(*proc_ents.rbegin()); 
    RRA("Couldn't get dimension.");
    
  }

    // proc_ents should all be of same dimension
  if (resolve_dim > shared_dim &&
      mbImpl->dimension_from_handle(*proc_ents.rbegin()) !=
      mbImpl->dimension_from_handle(*proc_ents.begin())) {
    MBRange::iterator lower = proc_ents.lower_bound(MBCN::TypeDimensionMap[0].first),
      upper = proc_ents.upper_bound(MBCN::TypeDimensionMap[resolve_dim-1].second);
    proc_ents.erase(lower, upper);
  }
  
    // must call even if we don't have any entities, to make sure
    // collective comm'n works
  return resolve_shared_ents(this_set, proc_ents, resolve_dim, shared_dim);
}
  
MBErrorCode MBParallelComm::resolve_shared_ents(MBEntityHandle this_set,
                                                MBRange &proc_ents,
                                                int resolve_dim,
                                                int shared_dim) 
{
  MBErrorCode result;
  if (debug) std::cerr << "Resolving shared entities." << std::endl;

  if (-1 == shared_dim) {
    if (0 == resolve_dim) {
      result = mbImpl->get_dimension(shared_dim); 
      RRA("Couldn't get dimension.");
    }
    else shared_dim = mbImpl->dimension_from_handle(*proc_ents.begin())-1;
  }
  assert(shared_dim >= 0 && resolve_dim >= 0);
  
    // get the skin entities by dimension
  MBRange skin_ents[4];
  std::vector<int> gid_data;
  std::vector<MBEntityHandle> handle_vec;
  int skin_dim;

    // get the entities to be skinned
  if (resolve_dim < shared_dim) {
      // for vertex-based partition, it's the elements adj to the vertices
    result = mbImpl->get_adjacencies(proc_ents, shared_dim,
                                     false, skin_ents[resolve_dim],
                                     MBInterface::UNION);
    RRA("Failed getting skinned entities.");
    skin_dim = shared_dim-1;
  }
  else {
      // for element-based partition, it's just the elements
    skin_ents[resolve_dim] = proc_ents;
    skin_dim = resolve_dim-1;
  }

    // find the skin
  MBSkinner skinner(mbImpl);
  result = skinner.find_skin(skin_ents[skin_dim+1], skin_ents[skin_dim],
                             skin_ents[skin_dim], true);
  RRA("Failed to find skin.");
  if (debug) std::cerr << "Found skin, now resolving." << std::endl;

    // get entities adjacent to skin ents from shared_dim down to
    // zero; don't create them if they don't exist already
  for (int this_dim = skin_dim-1; this_dim >= 0; this_dim--) {
    result = mbImpl->get_adjacencies(skin_ents[skin_dim], this_dim,
                                     false, skin_ents[this_dim],
                                     MBInterface::UNION);
    RRA("Failed getting skin adjacencies.");
  }

    // resolve shared vertices first

    // global id tag
  MBTag gid_tag; int def_val = -1;
  result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int),
                              MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                              &def_val, true);
  if (MB_FAILURE == result) return result;

  else if (MB_ALREADY_ALLOCATED != result) {
      // just created it, so we need global ids
    result = assign_global_ids(0, skin_dim+1);
    RRA("Failed assigning global ids.");
  }

    // store index in temp tag; reuse gid_data 
  gid_data.resize(2*skin_ents[0].size());
  int idx = 0;
  for (MBRange::iterator rit = skin_ents[0].begin(); 
       rit != skin_ents[0].end(); rit++) 
    gid_data[idx] = idx, idx++;
  MBTag idx_tag;
  result = mbImpl->tag_create("__idx_tag", sizeof(int), MB_TAG_DENSE,
                              MB_TYPE_INTEGER, idx_tag, &def_val, true);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) return result;
  result = mbImpl->tag_set_data(idx_tag, skin_ents[0], &gid_data[0]);
  RRA("Couldn't assign index tag.");

    // get gids for skin ents in a vector, to pass to gs
  result = mbImpl->tag_get_data(gid_tag, skin_ents[0], &gid_data[0]);
  RRA("Couldn't get gid tag for skin vertices.");

    // put handles in vector for passing to gs setup
  std::copy(skin_ents[0].begin(), skin_ents[0].end(), 
            std::back_inserter(handle_vec));
  
    // get a crystal router
  crystal_data *cd = procConfig.crystal_router();

/*  
    // get total number of entities; will overshoot highest global id, but
    // that's ok
  int num_total[2] = {0, 0}, num_local[2] = {0, 0};
  result = mbImpl->get_number_entities_by_dimension(0, 0, num_local);
  if (MB_SUCCESS != result) return result;
  int failure = MPI_Allreduce(num_local, num_total, 1,
                              MPI_INTEGER, MPI_SUM, procConfig.proc_comm());
  if (failure) {
    result = MB_FAILURE;
    RRA("Allreduce for total number of shared ents failed.");
  }
  
*/
    // call gather-scatter to get shared ids & procs
  gs_data *gsd;
  assert(sizeof(ulong_) == sizeof(MBEntityHandle));
  if (sizeof(int) != sizeof(ulong_)) {
    std::vector<long> lgid_data(gid_data.size());
    std::copy(gid_data.begin(), gid_data.end(), lgid_data.begin());
    gsd = gs_data_setup(skin_ents[0].size(), &lgid_data[0], 
                        (ulong_*)&handle_vec[0], 2, 1, 1, cd);
  }
  else {
    gsd = gs_data_setup(skin_ents[0].size(), (long*)&gid_data[0], 
                        (ulong_*)&handle_vec[0], 2, 1, 1, cd);
  }
  
  if (NULL == gsd) {
    result = MB_FAILURE;
    RRA("Couldn't create gs data.");
  }

    // get shared proc tags
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                sharedh_tag, sharedhs_tag, pstatus_tag);
  RRA("Couldn't get shared proc tags.");
  
    // load shared verts into a tuple, then sort by index
  tuple_list shared_verts;
  tuple_list_init_max(&shared_verts, 2, 0, 1, 0, 
                      skin_ents[0].size()*(MAX_SHARING_PROCS+1));
  unsigned int i = 0, j = 0;
  for (unsigned int p = 0; p < gsd->nlinfo->np; p++) 
    for (unsigned int np = 0; np < gsd->nlinfo->nshared[p]; np++) {
      shared_verts.vi[i++] = gsd->nlinfo->sh_ind[j];
      shared_verts.vi[i++] = gsd->nlinfo->target[p];
      shared_verts.vul[j] = gsd->nlinfo->ulabels[j];
      j++;
      shared_verts.n++;
    }
  
  int max_size = skin_ents[0].size()*(MAX_SHARING_PROCS+1);
  buffer sort_buffer;
  buffer_init(&sort_buffer, max_size);
  tuple_list_sort(&shared_verts, 0, &sort_buffer);
  buffer_free(&sort_buffer);

    // set sharing procs and handles tags on skin ents
  int maxp = -1;
  std::vector<int> sharing_procs(MAX_SHARING_PROCS);
  std::fill(sharing_procs.begin(), sharing_procs.end(), maxp);
  j = 0; i = 0;

    // get ents shared by 1 or n procs
  std::map<std::vector<int>, MBRange> proc_nranges;
  MBRange proc_verts;
  result = mbImpl->get_adjacencies(proc_ents, 0, false, proc_verts,
                                   MBInterface::UNION);
  RRA("Couldn't get proc_verts.");
  
  result = tag_shared_verts(shared_verts, skin_ents,
                            proc_nranges, proc_verts);
  RRA("Trouble tagging shared verts.");

    // get entities shared by 1 or n procs
  result = tag_shared_ents(resolve_dim, shared_dim, skin_ents,
                           proc_nranges);
  RRA("Trouble tagging shared entities.");

  tuple_list_free(&shared_verts);
  
  if (debug) {
    for (std::map<std::vector<int>, MBRange>::const_iterator mit = proc_nranges.begin();
         mit != proc_nranges.end(); mit++) {
      std::cout << "Iface: ";
      for (std::vector<int>::const_iterator vit = (mit->first).begin();
           vit != (mit->first).end(); vit++) std::cout << " " << *vit;
      std::cout << std::endl;
    }
  }
  
    // create the sets for each interface; store them as tags on
    // the interface instance
  MBRange iface_sets;
  result = create_interface_sets(proc_nranges, resolve_dim, shared_dim);
  RRA("Trouble creating iface sets.");

    // establish comm procs and buffers for them
  std::set<unsigned int> procs;
  result = get_interface_procs(procs, true);
  RRA("Trouble getting iface procs.");

    // resolve shared entity remote handles; implemented in ghost cell exchange
    // code because it's so similar
  result = exchange_ghost_cells(-1, -1, 0, true, true);
  RRA("Trouble resolving shared entity remote handles.");

    // now set the shared/interface tag on non-vertex entities on interface
  result = tag_iface_entities();
  RRA("Failed to tag iface entities.");

    // now build parent/child links for interface sets
  result = create_iface_pc_links();
  RRA("Trouble creating interface parent/child links.");

  gs_data_free(gsd);

    // done
  return result;
}

MBErrorCode MBParallelComm::resolve_shared_ents(MBParallelComm **pc, 
                                                const unsigned int np, 
                                                const int part_dim) 
{
  std::vector<MBRange> verts(np);
  int tot_verts = 0;
  unsigned int p, i, j, v, vtot;
  MBErrorCode rval;
  for (p = 0; p < np; p++) {
    MBSkinner skinner(pc[p]->get_moab());
    MBRange part_ents, skin_ents;
    rval = pc[p]->get_moab()->get_entities_by_dimension(0, part_dim, part_ents);
    if (MB_SUCCESS != rval) return rval;
    rval = skinner.find_skin(part_ents, skin_ents, skin_ents, true);
    if (MB_SUCCESS != rval) return rval;
    rval = pc[p]->get_moab()->get_adjacencies(skin_ents, 0, true, verts[p],
                                              MBInterface::UNION);
    if (MB_SUCCESS != rval) return rval;
    tot_verts += verts[p].size();
  }
  
  tuple_list shared_ents;
  tuple_list_init_max(&shared_ents, 2, 0, 1, 0, tot_verts);

  i = 0; j = 0;
  std::vector<int> gids;
  MBRange::iterator rit;
  MBTag gid_tag;
  int dum_default = -1;
  for (p = 0; p < np; p++) {
    rval = pc[p]->get_moab()->tag_create(GLOBAL_ID_TAG_NAME, 
                                         sizeof(int), MB_TAG_DENSE,
                                         MB_TYPE_INTEGER, gid_tag, 
                                         &dum_default, true);
    gids.resize(verts[p].size());
    rval = pc[p]->get_moab()->tag_get_data(gid_tag, verts[p], &gids[0]);
    if (MB_SUCCESS != rval) return rval;
    
    for (v = 0, rit = verts[p].begin(); v < gids.size(); v++, rit++) {
      shared_ents.vi[i++] = gids[v];
      shared_ents.vi[i++] = p;
      shared_ents.vul[j] = *rit;
      j++;
      shared_ents.n++;
    }
  }
  
  buffer sort_buffer;
  buffer_init(&sort_buffer, vtot);
  tuple_list_sort(&shared_ents, 0, &sort_buffer);
  buffer_free(&sort_buffer);

  j = 0; i = 0;
  std::vector<MBEntityHandle> handles;
  std::vector<int> procs;
  
  while (i < shared_ents.n) {
    handles.clear();
    procs.clear();
    
      // count & accumulate sharing procs
    int this_gid = shared_ents.vi[j];
    while (i < shared_ents.n && shared_ents.vi[j] == this_gid) {
      j++;
      procs.push_back( shared_ents.vi[j++] );
      handles.push_back( shared_ents.vul[i++] );
    }
    if (1 == procs.size()) continue;
    
    for (v = 0; v < procs.size(); v++) {
      rval = pc[procs[v]]->update_remote_data(handles[v], 
                                              &procs[0], &handles[0], procs.size(),
                                              PSTATUS_INTERFACE);
      if (MB_SUCCESS != rval) return rval;
    }
  }

  std::set<unsigned int> psets;
  for (p = 0; p < np; p++) {
    rval = pc[p]->create_interface_sets(part_dim, part_dim-1);
    if (MB_SUCCESS != rval) return rval;
      // establish comm procs and buffers for them
    psets.clear();
    rval = pc[p]->get_interface_procs(psets, true);
    if (MB_SUCCESS != rval) return rval;
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::tag_iface_entities() 
{
  MBErrorCode result = MB_SUCCESS;
  MBRange iface_ents, tmp_ents, rmv_ents;
  std::vector<unsigned char> pstat;
  MBRange::iterator rit2;
  unsigned int i;
  
  for (MBRange::iterator rit = interfaceSets.begin(); rit != interfaceSets.end(); rit++) {
    iface_ents.clear();
    
    result = mbImpl->get_entities_by_handle(*rit, iface_ents);
    RRA("Couldn't get iface set contents.");
    pstat.resize(iface_ents.size());
    result = mbImpl->tag_get_data(pstatus_tag(), iface_ents, &pstat[0]);
    RRA("Couldn't get pstatus values.");
    rmv_ents.clear();
    for (rit2 = iface_ents.begin(), i = 0; rit2 != iface_ents.end(); rit2++, i++) {
      if (!(pstat[i] & PSTATUS_INTERFACE)) rmv_ents.insert(*rit2);
    }
    result = mbImpl->remove_entities(*rit, rmv_ents);
    RRA("Couldn't remove entities from set.");
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::set_pstatus_entities(MBRange &pstatus_ents,
                                                 unsigned char pstatus_val,
                                                 bool lower_dim_ents,
                                                 bool verts_too,
                                                 int operation) 
{
  std::vector<unsigned char> pstatus_vals(pstatus_ents.size());
  MBRange all_ents, *range_ptr = &pstatus_ents;
  MBErrorCode result;
  if (lower_dim_ents || verts_too) {
    all_ents = pstatus_ents;
    range_ptr = &all_ents;
    int start_dim = (lower_dim_ents ? mbImpl->dimension_from_handle(*pstatus_ents.rbegin())-1 : 0);
    for (; start_dim >= 0; start_dim--) {
      result = mbImpl->get_adjacencies(all_ents, start_dim, true, all_ents,
                                       MBInterface::UNION);
      RRA(" ");
    }
  }
  if (MBInterface::UNION == operation) {
    result = mbImpl->tag_get_data(pstatus_tag(), *range_ptr, &pstatus_vals[0]);
    RRA("Couldn't get pstatus tag value.");
    for (unsigned int i = 0; i < pstatus_vals.size(); i++)
      pstatus_vals[i] |= pstatus_val;
  }
  else {
    for (unsigned int i = 0; i < pstatus_vals.size(); i++)
      pstatus_vals[i] = pstatus_val;
  }
  result = mbImpl->tag_set_data(pstatus_tag(), *range_ptr, &pstatus_vals[0]);
  RRA("Couldn't set pstatus tag value.");
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::set_pstatus_entities(MBEntityHandle *pstatus_ents,
                                                 int num_ents,
                                                 unsigned char pstatus_val,
                                                 bool lower_dim_ents,
                                                 bool verts_too,
                                                 int operation) 
{
  std::vector<unsigned char> pstatus_vals(num_ents);
  MBErrorCode result;
  if (lower_dim_ents || verts_too) {
      // in this case, call the range-based version
    MBRange tmp_range;
    std::copy(pstatus_ents, pstatus_ents+num_ents, mb_range_inserter(tmp_range));
    return set_pstatus_entities(tmp_range, pstatus_val, lower_dim_ents, 
                                verts_too, operation);
  }

  if (MBInterface::UNION == operation) {
    result = mbImpl->tag_get_data(pstatus_tag(), pstatus_ents, num_ents, &pstatus_vals[0]);
    RRA("Couldn't get pstatus tag value.");
    for (unsigned int i = 0; i < (unsigned int) num_ents; i++)
      pstatus_vals[i] |= pstatus_val;
  }
  else {
    for (unsigned int i = 0; i < (unsigned int) num_ents; i++)
      pstatus_vals[i] = pstatus_val;
  }
  result = mbImpl->tag_set_data(pstatus_tag(), pstatus_ents, num_ents, &pstatus_vals[0]);
  RRA("Couldn't set pstatus tag value.");
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::create_interface_sets(int resolve_dim, int shared_dim) 
{
  std::map<std::vector<int>, MBRange> proc_nranges;
  
    // build up the list of shared entities
  int procs[MAX_SHARING_PROCS];
  MBEntityHandle handles[MAX_SHARING_PROCS];
  MBErrorCode result;
  int nprocs;
  unsigned char pstat;
  for (MBRange::iterator rit = sharedEnts.begin(); rit != sharedEnts.end(); rit++) {
    if (shared_dim != -1 && mbImpl->dimension_from_handle(*rit) > shared_dim)
      continue;
    result = get_sharing_data(*rit, procs, handles, pstat, nprocs);
    RRA("");
    std::sort(procs, procs+nprocs);
    std::vector<int> tmp_procs(procs, procs + nprocs);
    proc_nranges[tmp_procs].insert(*rit);
  }
                                                  
  MBSkinner skinner(mbImpl);
  MBRange skin_ents[4];
  result = mbImpl->get_entities_by_dimension(0, resolve_dim, skin_ents[resolve_dim]);
  RRA("");
  result = skinner.find_skin(skin_ents[resolve_dim], skin_ents[resolve_dim-1], 
                             skin_ents[resolve_dim-1], true);
  RRA("Failed to find skin.");
  if (shared_dim > 1) {
    result = mbImpl->get_adjacencies(skin_ents[resolve_dim-1], resolve_dim-2, true,
                                     skin_ents[resolve_dim-2], MBInterface::UNION);
    RRA("");
  }

  result = tag_shared_ents(resolve_dim, shared_dim, skin_ents,
                           proc_nranges);
    
  return create_interface_sets(proc_nranges, resolve_dim, shared_dim);
}
  
MBErrorCode MBParallelComm::create_interface_sets(std::map<std::vector<int>, MBRange> &proc_nranges,
                                                  int resolve_dim, int shared_dim) 
{
  if (proc_nranges.empty()) return MB_SUCCESS;
  
  int proc_ids[MAX_SHARING_PROCS];
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag,
                                            pstatus_tag);
  RRA("Trouble getting shared proc tags in create_interface_sets.");
  MBRange::iterator rit;

    // create interface sets, tag them, and tag their contents with iface set tag
  std::vector<MBEntityHandle> tag_vals;
  std::vector<unsigned char> pstatus;
  for (std::map<std::vector<int>,MBRange>::iterator mit = proc_nranges.begin();
       mit != proc_nranges.end(); mit++) {
      // create the set
    MBEntityHandle new_set;
    result = mbImpl->create_meshset(MESHSET_SET, new_set); 
    RRA("Failed to create interface set.");
    interfaceSets.insert(new_set);

      // add entities
    result = mbImpl->add_entities(new_set, mit->second); 
    RRA("Failed to add entities to interface set.");
      // tag set with the proc rank(s)
    if (mit->first.size() == 1)
      result = mbImpl->tag_set_data(sharedp_tag, &new_set, 1, 
                                    &(mit->first)[0]); 
    else {
      // pad tag data out to MAX_SHARING_PROCS with -1
      assert( mit->first.size() <= MAX_SHARING_PROCS );
      std::copy( mit->first.begin(), mit->first.end(), proc_ids );
      std::fill( proc_ids + mit->first.size(), proc_ids + MAX_SHARING_PROCS, -1 );
      result = mbImpl->tag_set_data(sharedps_tag, &new_set, 1, proc_ids );
    }
    RRA("Failed to tag interface set with procs.");
    
      // get the owning proc, then set the pstatus tag on iface set
    int min_proc = (mit->first)[0];
    unsigned char pval = (PSTATUS_SHARED | PSTATUS_INTERFACE);
    if (min_proc < (int) procConfig.proc_rank()) pval |= PSTATUS_NOT_OWNED;
    if (mit->first.size() > 1) pval |= PSTATUS_MULTISHARED;
    result = mbImpl->tag_set_data(pstatus_tag, &new_set, 1, &pval); 
    RRA("Failed to tag interface set with pstatus.");

      // tag the vertices with the same thing
    pstatus.clear();
    MBRange verts = (mit->second).subset_by_type(MBVERTEX);
    pstatus.resize(verts.size(), pval);
    result = mbImpl->tag_set_data(pstatus_tag, verts, &pstatus[0]); 
    RRA("Failed to tag interface set vertices with pstatus.");
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::create_iface_pc_links() 
{
    // now that we've resolved the entities in the iface sets, 
    // set parent/child links between the iface sets

    // first tag all entities in the iface sets
  MBTag tmp_iface_tag;
  MBEntityHandle tmp_iface_set = 0;
  MBErrorCode result = mbImpl->tag_create("__tmp_iface", sizeof(MBEntityHandle),
                                          MB_TAG_DENSE, MB_TYPE_HANDLE,
                                          tmp_iface_tag, &tmp_iface_set);
  if (MB_ALREADY_ALLOCATED != result && MB_SUCCESS != result) 
    RRA("Failed to create temporary iface set tag.");

  MBRange iface_ents;
  std::vector<MBEntityHandle> tag_vals;
  MBRange::iterator rit;
  
  for (rit = interfaceSets.begin(); rit != interfaceSets.end(); rit++) {
      // tag entities with interface set
    iface_ents.clear();
    result = mbImpl->get_entities_by_handle(*rit, iface_ents);
    RRA("Couldn't get entities in iface set.");
    
    if (iface_ents.empty()) continue;
    
    tag_vals.resize(iface_ents.size());
    std::fill(tag_vals.begin(), tag_vals.end(), *rit);
    result = mbImpl->tag_set_data(tmp_iface_tag, iface_ents, &tag_vals[0]); 
    RRA("Failed to tag iface entities with interface set.");
  }
  
    // now go back through interface sets and add parent/child links
  MBRange tmp_ents2;
  for (int d = 2; d >= 0; d--) {
    for (rit = interfaceSets.begin(); rit != interfaceSets.end(); rit++) {
        // get entities on this interface
      iface_ents.clear();
      result = mbImpl->get_entities_by_handle(*rit, iface_ents, true);
      RRA("Couldn't get entities by dimension.");
      if (iface_ents.empty() ||
          mbImpl->dimension_from_handle(*iface_ents.rbegin()) != d) continue;

        // get higher-dimensional entities and their interface sets
      result = mbImpl->get_adjacencies(&(*iface_ents.begin()), 1, d+1,
                                       false, tmp_ents2);
      RRA("Couldn't get adjacencies for interface sets.");
      tag_vals.resize(tmp_ents2.size());
      result = mbImpl->tag_get_data(tmp_iface_tag, tmp_ents2, &tag_vals[0]);
      RRA("Couldn't get iface set tag for interface sets.");
      
        // go through and for any on interface make it a parent
      MBEntityHandle last_set = 0;
      for (unsigned int i = 0; i < tag_vals.size(); i++) {
        if (tag_vals[i] && tag_vals[i] != last_set) {
          result = mbImpl->add_parent_child(tag_vals[i], *rit);
          RRA("Couldn't add parent/child link for interface set.");
          last_set = tag_vals[i];
        }
      }
    }
  }
  
    // delete the temporary tag
  result = mbImpl->tag_delete(tmp_iface_tag);
  RRA("Couldn't delete tmp iface tag.");

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::tag_shared_ents(int resolve_dim,
                                            int shared_dim,
                                            MBRange *skin_ents,
                                            std::map<std::vector<int>, MBRange> &proc_nranges) 
{
    // set sharing procs tags on other skin ents
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag, pstatus_tag);
  RRA("Trouble getting shared proc tags in tag_shared_ents.");
  const MBEntityHandle *connect; int num_connect;
  std::vector<int> sharing_procs, sharing_procs1, sharing_procs2;
  std::vector<int>::iterator vii;
  std::vector<unsigned char> pstatus_flags;
  std::vector<MBEntityHandle> dum_connect;

  for (int d = 3; d > 0; d--) {
    if (resolve_dim == d) continue;
    
    for (MBRange::iterator rit = skin_ents[d].begin();
         rit != skin_ents[d].end(); rit++) {
        // get connectivity
      result = mbImpl->get_connectivity(*rit, connect, num_connect, true,
                                        &dum_connect);
      RRA("Failed to get connectivity on non-vertex skin entities.");
 
        // if any vertices not shared, this entity isn't
      bool is_shared = true;
      pstatus_flags.resize( num_connect );
      result = mbImpl->tag_get_data(pstatus_tag, connect, num_connect,
                                    &pstatus_flags[0]);
      RRA("Couldn't get pstatus flag.");
      for (int nc = 0; nc < num_connect; nc++) {
        if (!(pstatus_flags[nc] & PSTATUS_SHARED)) {
          is_shared = false;
          break;
        }
      }
      if (!is_shared) continue;

      for (int nc = 0; nc < num_connect; nc++) {
        sharing_procs2.clear();
        
          // get sharing procs
        sharing_procs2.resize(1);
        result = mbImpl->tag_get_data(sharedp_tag, connect+nc, 1, &sharing_procs2[0]);
        RRA("Couldn't get sharedp_tag on skin vertices in entity.");
        if (sharing_procs2[0] == -1) {
          sharing_procs2.resize(MAX_SHARING_PROCS);
          result = mbImpl->tag_get_data(sharedps_tag, connect+nc, 1, &sharing_procs2[0]);
          RRA("Couldn't get sharedps_tag on skin vertices in entity.");
        }
        assert(-1 != sharing_procs2[0]);
          // remove any unnecessary entries
        vii = std::find( sharing_procs2.begin(), sharing_procs2.end(), -1 );
        sharing_procs2.erase( vii, sharing_procs2.end() );
        
          // build range of sharing procs for this vertex
          // intersect with range for this skin ent
        if (0 == nc) {
          sharing_procs.swap( sharing_procs2 );
        }
        else if (resolve_dim < shared_dim) {
          sharing_procs1.clear();
          set_union( sharing_procs.begin(), sharing_procs.end(), 
                     sharing_procs2.begin(), sharing_procs2.end(),
                     std::back_inserter( sharing_procs1 ) );
          sharing_procs.swap( sharing_procs1 );
        }
        else {
          sharing_procs1.clear();
          set_intersection( sharing_procs.begin(), sharing_procs.end(), 
                            sharing_procs2.begin(), sharing_procs2.end(),
                            std::back_inserter( sharing_procs1 ) );
          sharing_procs.swap( sharing_procs1 );
        }
      }

      if (sharing_procs.empty() && resolve_dim < shared_dim) continue;

        // intersection is the owning proc(s) for this skin ent
      if (sharing_procs.empty()) continue;

      proc_nranges[sharing_procs].insert(*rit);

        // reset sharing proc(s) tags
      sharing_procs.clear();
    }
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::tag_shared_verts(tuple_list &shared_ents,
                                             MBRange *skin_ents,
                                             std::map<std::vector<int>, MBRange> &proc_nranges,
                                             MBRange &proc_verts) 
{
  MBTag sharedp_tag, sharedps_tag, sharedh_tag, sharedhs_tag, pstatus_tag;
  MBErrorCode result = get_shared_proc_tags(sharedp_tag, sharedps_tag, 
                                            sharedh_tag, sharedhs_tag, pstatus_tag);
  RRA("Trouble getting shared proc tags in tag_shared_verts.");
  
  unsigned int j = 0, i = 0;
  std::vector<int> sharing_procs, sharing_procs2;
  std::vector<MBEntityHandle> sharing_handles, sharing_handles2;
  
  while (j < 2*shared_ents.n) {
      // count & accumulate sharing procs
    int this_idx = shared_ents.vi[j];
    MBEntityHandle this_ent = skin_ents[0][this_idx];
    while (j < 2*shared_ents.n && shared_ents.vi[j] == this_idx) {
      j++;
      sharing_procs.push_back( shared_ents.vi[j++] );
      sharing_handles.push_back( shared_ents.vul[i++] );
    }

    if (sharing_procs.size() > 1) {
        // add current proc/handle to list
      sharing_procs.push_back(procConfig.proc_rank());
      sharing_handles.push_back(this_ent);
    }
      
      // sort sharing_procs and sharing_handles such that
      // sharing_procs is in ascending order.  Use temporary
      // lists and binary search to re-order sharing_handles.
    sharing_procs2 = sharing_procs;
    std::sort( sharing_procs2.begin(), sharing_procs2.end() );
    sharing_handles2.resize( sharing_handles.size() );
    for (size_t k = 0; k < sharing_handles.size(); ++k) {
      size_t idx = std::lower_bound( sharing_procs2.begin(), 
                                     sharing_procs2.end(), 
                                     sharing_procs[k] ) - sharing_procs2.begin();
      sharing_handles2[idx] = sharing_handles[k];
    }
    sharing_procs.swap( sharing_procs2 );
    sharing_handles.swap( sharing_handles2 );
    
    
    proc_nranges[sharing_procs].insert(this_ent);

    unsigned char share_flag = PSTATUS_SHARED, 
        ms_flag = (PSTATUS_SHARED | PSTATUS_MULTISHARED);
    if (sharing_procs.size() == 1) {
      result = mbImpl->tag_set_data(sharedp_tag, &this_ent, 1,
                                    &sharing_procs[0]);
      result = mbImpl->tag_set_data(sharedh_tag, &this_ent, 1,
                                    &sharing_handles[0]);
      result = mbImpl->tag_set_data(pstatus_tag, &this_ent, 1, &share_flag);
      RRA("Couldn't set shared tag on shared vertex.");
      sharedEnts.insert(this_ent);
    }
    else {
        // pad lists 
      assert( sharing_procs.size() <= MAX_SHARING_PROCS );
      sharing_procs.resize( MAX_SHARING_PROCS, -1 );
      sharing_handles.resize( MAX_SHARING_PROCS, 0 );
      result = mbImpl->tag_set_data(sharedps_tag, &this_ent, 1,
                                    &sharing_procs[0]);
      result = mbImpl->tag_set_data(sharedhs_tag, &this_ent, 1,
                                    &sharing_handles[0]);
      result = mbImpl->tag_set_data(pstatus_tag, &this_ent, 1, &ms_flag);
      RRA("Couldn't set multi-shared tag on shared vertex.");
      sharedEnts.insert(this_ent);
    }
    RRA("Failed setting shared_procs tag on skin vertices.");

      // reset sharing proc(s) tags
    sharing_procs.clear();
    sharing_handles.clear();
  }

  return MB_SUCCESS;
}
  
  //! get processors with which this processor communicates; sets are sorted by processor
MBErrorCode MBParallelComm::get_interface_procs(std::set<unsigned int> &procs_set,
                                                bool get_buffs)
{
    // make sure the sharing procs vector is empty
  procs_set.clear();

    // pre-load vector of single-proc tag values
  unsigned int i, j;
  std::vector<int> iface_proc(interfaceSets.size());
  MBErrorCode result = mbImpl->tag_get_data(sharedp_tag(), interfaceSets, &iface_proc[0]);
  RRA("Failed to get iface_proc for iface sets.");

    // get sharing procs either from single-proc vector or by getting
    // multi-proc tag value
  int tmp_iface_procs[MAX_SHARING_PROCS];
  std::fill(tmp_iface_procs, tmp_iface_procs+MAX_SHARING_PROCS, -1);
  MBRange::iterator rit;
  for (rit = interfaceSets.begin(), i = 0; rit != interfaceSets.end(); rit++, i++) {
    if (-1 != iface_proc[i]) procs_set.insert((unsigned int) iface_proc[i]);
    else {
        // get the sharing_procs tag
      result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1,
                                    tmp_iface_procs);
      RRA("Failed to get iface_procs for iface set.");
      for (j = 0; j < MAX_SHARING_PROCS; j++) {
        if (-1 != tmp_iface_procs[j] && tmp_iface_procs[j] != (int)procConfig.proc_rank()) 
          procs_set.insert((unsigned int) tmp_iface_procs[j]);
        else if (-1 == tmp_iface_procs[j]) {
          std::fill(tmp_iface_procs, tmp_iface_procs+j, -1);
          break;
        }
      }
    }
  }

  if (get_buffs) {
    for (std::set<unsigned int>::iterator sit = procs_set.begin(); sit != procs_set.end(); sit++)
      get_buffers(*sit);
  }
  
  return MB_SUCCESS;
}
  
MBErrorCode MBParallelComm::get_pstatus_entities(int dim,
                                                 unsigned char pstatus_val,
                                                 MBRange &pstatus_ents)
{
  MBRange ents;
  MBErrorCode result;
  
  if (-1 == dim) result = mbImpl->get_entities_by_handle(0, ents);
  else result = mbImpl->get_entities_by_dimension(0, dim, ents);
  RRA(" ");
  
  std::vector<unsigned char> pstatus(ents.size());
  result = mbImpl->tag_get_data(pstatus_tag(), ents, &pstatus[0]);
  RRA("Couldn't get pastatus tag.");
  MBRange::iterator rit = ents.begin();
  int i = 0;
  if (pstatus_val) {
    for (; rit != ents.end(); i++, rit++)
      if (pstatus[i]&pstatus_val &&
          (-1 == dim || mbImpl->dimension_from_handle(*rit) == dim)) 
        pstatus_ents.insert(*rit);
  }
  else {
    for (; rit != ents.end(); i++, rit++)
      if (!pstatus[i] &&
          (-1 == dim || mbImpl->dimension_from_handle(*rit) == dim)) 
        pstatus_ents.insert(*rit);
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::check_global_ids(MBEntityHandle this_set,
                                             const int dimension, 
                                             const int start_id,
                                             const bool largest_dim_only,
                                             const bool parallel)
{
    // global id tag
  MBTag gid_tag; int def_val = -1;
  MBErrorCode result = mbImpl->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int),
                                          MB_TAG_DENSE, MB_TYPE_INTEGER, gid_tag,
                                          &def_val, true);
  if (MB_ALREADY_ALLOCATED != result &&
      MB_SUCCESS != result) {
    RRA("Failed to create/get gid tag handle.");
  }

  MBRange dum_range;
  if (MB_ALREADY_ALLOCATED == result) {
    void *tag_ptr = &def_val;
    MBErrorCode tmp_result = mbImpl->get_entities_by_type_and_tag(this_set, MBVERTEX, 
                                                                  &gid_tag, &tag_ptr, 1,
                                                                  dum_range);
    if (MB_SUCCESS != tmp_result) {
      result = tmp_result;
      RRA("Failed to get gid tag.");
    }
  }
  
  if (MB_ALREADY_ALLOCATED != result || !dum_range.empty()) {
      // just created it, so we need global ids
    result = assign_global_ids(this_set, dimension, start_id, largest_dim_only,
                               parallel);
    RRA("Failed assigning global ids.");
  }

  return MB_SUCCESS;
}

bool MBParallelComm::is_iface_proc(MBEntityHandle this_set,
                                   int to_proc) 
{
  int sharing_procs[MAX_SHARING_PROCS];
  std::fill(sharing_procs, sharing_procs+MAX_SHARING_PROCS, -1);
  MBErrorCode result = mbImpl->tag_get_data(sharedp_tag(), &this_set, 1,
                                            sharing_procs);
  if (to_proc == sharing_procs[0]) return true;
  
  result = mbImpl->tag_get_data(sharedps_tag(), &this_set, 1,
                                sharing_procs);
  for (int i = 0; i < MAX_SHARING_PROCS; i++) {
    if (to_proc == sharing_procs[i]) return true;
    else if (-1 == sharing_procs[i]) return false;
  }
  
  return false;
}

MBErrorCode MBParallelComm::filter_pstatus( MBRange &ents,
                                            unsigned char pstat,
                                            unsigned char op,
                                            int to_proc,
                                            MBRange *returned_ents)
{
  MBRange tmp_ents;

  assert(!ents.empty());

    // Put into tmp_ents any entities which are not owned locally or
    // who are already shared with to_proc
  std::vector<unsigned char> shared_flags(ents.size()), shared_flags2;
  MBErrorCode result = mbImpl->tag_get_data(pstatus_tag(), ents,
                                            &shared_flags[0]);
  RRA("Failed to get pstatus flag.");
  MBRange::const_iterator rit;
  int i;
  if (op == PSTATUS_OR) {
    for (rit = ents.begin(), i = 0; rit != ents.end(); rit++, i++) 
      if (((shared_flags[i] & ~pstat)^shared_flags[i]) & pstat) {
        tmp_ents.insert(*rit);
        if (-1 != to_proc) shared_flags2.push_back(shared_flags[i]);
      }
  }
  else if (op == PSTATUS_AND) {
    for (rit = ents.begin(), i = 0; rit != ents.end(); rit++, i++)
      if ((shared_flags[i] & pstat) == pstat) {
        tmp_ents.insert(*rit);
        if (-1 != to_proc) shared_flags2.push_back(shared_flags[i]);
      }
  }
  else if (op == PSTATUS_NOT) {
    for (rit = ents.begin(), i = 0; rit != ents.end(); rit++, i++)
      if (!(shared_flags[i] & pstat)) {
        tmp_ents.insert(*rit);
        if (-1 != to_proc) shared_flags2.push_back(shared_flags[i]);
      }
  }
  else {
    assert(false);
    return MB_FAILURE;
  }

  if (-1 != to_proc) {

    int sharing_procs[MAX_SHARING_PROCS];
    std::fill(sharing_procs, sharing_procs+MAX_SHARING_PROCS, -1);
    MBRange tmp_ents2;

    for (rit = tmp_ents.begin(), i = 0; rit != tmp_ents.end(); rit++, i++) {
        // we need to check sharing procs
      if (shared_flags2[i] & PSTATUS_MULTISHARED) {
        result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1,
                                      sharing_procs);
        assert(-1 != sharing_procs[0]);
        RRA(" ");
        for (unsigned int j = 0; j < MAX_SHARING_PROCS; j++) {
            // if to_proc shares this entity, add it to list
          if (sharing_procs[j] == to_proc) {
            tmp_ents2.insert(*rit);
          }
          else if (sharing_procs[j] == -1) break;

          sharing_procs[j] = -1;
        }
      }
      else if (shared_flags2[i] & PSTATUS_SHARED) {
        result = mbImpl->tag_get_data(sharedp_tag(), &(*rit), 1,
                                      sharing_procs);
        RRA(" ");
        assert(-1 != sharing_procs[0]);
        if (sharing_procs[0] == to_proc) tmp_ents2.insert(*rit);
        sharing_procs[0] = -1;
      }
      else
        assert("should never get here" && false);
    }

    tmp_ents.swap(tmp_ents2);
  }
  
  if (returned_ents)
    returned_ents->swap(tmp_ents);
  else
    ents.swap(tmp_ents);
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::exchange_ghost_cells(int ghost_dim, int bridge_dim,
                                                 int num_layers,
                                                 bool store_remote_handles,
                                                 bool wait_all)
{
    // if we're only finding out about existing ents, we have to be storing
    // remote handles too
  assert(num_layers > 0 || store_remote_handles);
  
  const bool is_iface = !num_layers;
  
    // get the b-dimensional interface(s) with with_proc, where b = bridge_dim
  
  int success;
  unsigned char *buff_ptr;
  MBErrorCode result = MB_SUCCESS;

    // when this function is called, buffProcs should already have any 
    // communicating procs

    //===========================================
    // post ghost irecv's for ghost entities from all communicating procs
    //===========================================
    // index reqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_reqs(buffProcs.size(), MPI_REQUEST_NULL),
      send_reqs(buffProcs.size(), MPI_REQUEST_NULL);
  std::vector<unsigned int>::iterator proc_it;
  int ind;
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, buffProcs[ind],
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    //===========================================
    // get entities to be sent to neighbors
    //===========================================

  MBRange sent_ents[MAX_SHARING_PROCS], allsent, tmp_range;
  std::vector<std::set<unsigned int> > entprocs(allsent.size());
  result = get_sent_ents(is_iface, bridge_dim, ghost_dim, num_layers,
                         sent_ents, allsent, entprocs);
  RRA("get_sent_ents failed.");
  
    //===========================================
    // pack and send ents from this proc to others
    //===========================================
    // initialize sendReqs
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {

      // buff_ptr points to the END (one past last occupied byte) of buffer
    buff_ptr = &ownerSBuffs[ind][0];

      // entities
    result = pack_entities(sent_ents[ind], ownerSBuffs[ind], buff_ptr,
                           store_remote_handles, buffProcs[ind], is_iface,
                           &entprocs, &allsent); 
    RRA("Packing entities failed.");

      // now we're ready to send the buffer
    result = send_buffer(*proc_it, &ownerSBuffs[ind][0], 
                         buff_ptr-&ownerSBuffs[ind][0], MB_MESG_ENTS,
                         send_reqs[ind]);
    RRA("Failed to Isend in ghost exchange.");
  }

    //===========================================
    // receive/unpack new entities
    //===========================================
    // number of incoming messages for ghosts is the number of procs we 
    // communicate with; for iface, it's the number of those with lower rank
  int num_incoming = buffProcs.size();
  std::vector<MPI_Status> status(buffProcs.size());
  std::vector<std::vector<MBEntityHandle> > recd_ents(num_incoming);
  std::vector<std::vector<MBEntityHandle> > L1hloc(buffProcs.size()), L1hrem(buffProcs.size());
  std::vector<std::vector<int> > L1p(buffProcs.size());
  std::vector<MBEntityHandle> L2hloc, L2hrem;
  std::vector<unsigned int> L2p;
  MBRange new_ents;
  
  while (num_incoming) {
      // wait for all recvs of ghost ents before proceeding,
      // b/c some procs may have sent to a 3rd proc ents owned by me;
    success = MPI_Waitany(buffProcs.size(), &recv_reqs[0], &ind, &status[0]);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
      // branch on message type
    if (MB_MESG_SIZE == status[0].MPI_TAG) {
        // incoming message just has size; resize buffer and re-call recv,
        // then re-increment incoming count
      int new_size = *((int*)&ghostRBuffs[ind][0]);
      assert(0 > new_size);
      result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind],
                              MB_MESG_ENTS);
      RRA("Failed to resize recv buffer.");
      num_incoming++;
    }
    else if (MB_MESG_ENTS == status[0].MPI_TAG) {
      
        // incoming ghost entities; unpack; returns entities received
        // both from sending proc and from owning proc (which may be different)
      unsigned char *buff_ptr = &ghostRBuffs[ind][0];
      result = unpack_entities(buff_ptr,
                               store_remote_handles, ind, is_iface,
                               L1hloc, L1hrem, L1p, L2hloc, L2hrem, L2p, new_ents);
      RRA("Failed to unpack entities.");
    }
    else {
      assert(false);
      return MB_FAILURE;
    }
  }

    // add requests for any new addl procs
  if (recv_reqs.size() != buffProcs.size()) {
    recv_reqs.resize(buffProcs.size());
    send_reqs.resize(buffProcs.size());
  }
    
  if (is_iface) {
#ifdef NDEBUG
    result = check_sent_ents(allsent);
    RRA("Failed check on shared entities.");
    result = check_all_shared_handles();
    RRA("Failed check on all shared handles.");
#endif
    return MB_SUCCESS;
  }
  
    //===========================================
    // post recvs for remote handles of my sent ents
    //===========================================
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {
      // skip if iface layer and lower-rank proc
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, buffProcs[ind],
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    //===========================================
    // send local handles for new ghosts to owner, then add
    // those to ghost list for that owner
    //===========================================
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {
      // skip if iface layer and higher-rank proc
    buff_ptr = &ghostSBuffs[ind][0];
    result = pack_remote_handles(L1hloc[ind], L1hrem[ind], L1p[ind], *proc_it,
                                   ghostSBuffs[ind], buff_ptr);
    RRA("Failed to pack remote handles.");
    result = send_buffer(buffProcs[ind], &ghostSBuffs[ind][0], 
                         buff_ptr - &ghostSBuffs[ind][0], 
                         MB_MESG_REMOTE_HANDLES, send_reqs[ind]);
    RRA("Failed to send remote handles.");
  }
  
    //===========================================
    // process remote handles of my ghosteds
    //===========================================
  num_incoming = buffProcs.size();
  while (num_incoming) {
    success = MPI_Waitany(buffProcs.size(), &recv_reqs[0], &ind, &status[0]);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
      // branch on message type
    if (MB_MESG_SIZE == status[0].MPI_TAG) {
        // incoming message just has size; resize buffer and re-call recv,
        // then re-increment incoming count
      int new_size = *((int*)&ghostRBuffs[ind][0]);
      assert(0 > new_size);
      result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind],
                              MB_MESG_REMOTE_HANDLES);
      RRA("Failed to resize recv buffer.");
      num_incoming++;
    }
    else if (MB_MESG_REMOTE_HANDLES == status[0].MPI_TAG) {
        // incoming remote handles
      result = unpack_remote_handles(buffProcs[ind], &ghostRBuffs[ind][0],
                                     L2hloc, L2hrem, L2p);
      RRA("Failed to unpack remote handles.");
    }
    else assert(false);
  }
    
#ifdef NDEBUG
  result = check_sent_ents(allsent);
  RRA("Failed check on shared entities.");
  result = check_all_shared_handles();
  RRA("Failed check on all shared handles.");
#endif

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_sent_ents(const bool is_iface, 
                                          const int bridge_dim, const int ghost_dim,
                                          const int num_layers,
                                          MBRange *sent_ents, MBRange &allsent,
                                          std::vector<std::set<unsigned int> > &entprocs) 
{
  MBErrorCode result;
  int ind;
  std::vector<unsigned int>::iterator proc_it;
  MBRange tmp_range;
  
    // done in a separate loop over procs because sometimes later procs 
    // need to add info to earlier procs' messages
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {
    if (!is_iface) {
      result = get_ghosted_entities(bridge_dim, ghost_dim, buffProcs[ind],
                                    num_layers, sent_ents[ind]);
      RRA("Failed to get ghost layers.");
    }
    else {
      result = get_iface_entities(buffProcs[ind], -1, sent_ents[ind]);
      RRA("Failed to get interface layers.");
    }

      // filter out entities already shared with destination
    tmp_range.clear();
    result = filter_pstatus(sent_ents[ind], PSTATUS_SHARED, PSTATUS_AND,
                            buffProcs[ind], &tmp_range);
    RRA("Couldn't filter on owner.");
    if (!tmp_range.empty()) 
      sent_ents[ind] = sent_ents[ind].subtract(tmp_range);

    allsent.merge(sent_ents[ind]);
  }

    //===========================================
    // need to get procs each entity is sent to
    //===========================================
  MBRange::iterator rit;
  entprocs.resize(allsent.size());
  for (ind = 0, proc_it = buffProcs.begin(); 
       proc_it != buffProcs.end(); proc_it++, ind++) {
    for (rit = sent_ents[ind].begin(); rit != sent_ents[ind].end(); rit++) {
      int rind = allsent.index(*rit);
      assert(rind < (int) allsent.size() && rind >= 0);
      entprocs[rind].insert(*proc_it);
    }
  }

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::exchange_ghost_cells(MBParallelComm **pcs,
                                                 unsigned int num_procs,
                                                 int ghost_dim, int bridge_dim,
                                                 int num_layers,
                                                 bool store_remote_handles)
{
    // static version of function, exchanging info through buffers rather 
    // than through messages

    // if we're only finding out about existing ents, we have to be storing
    // remote handles too
  assert(num_layers > 0 || store_remote_handles);
  
  const bool is_iface = !num_layers;
  
  unsigned int ind;
  unsigned char *buff_ptr;
  MBParallelComm *pc;
  MBErrorCode result = MB_SUCCESS;

    // when this function is called, buffProcs should already have any 
    // communicating procs

    //===========================================
    // get entities to be sent to neighbors
    //===========================================

    // done in a separate loop over procs because sometimes later procs 
    // need to add info to earlier procs' messages
  MBRange sent_ents[MAX_SHARING_PROCS][MAX_SHARING_PROCS], 
      allsent[MAX_SHARING_PROCS];

    //===========================================
    // get entities to be sent to neighbors
    //===========================================

  std::vector<std::set<unsigned int> > entprocs[MAX_SHARING_PROCS];
  for (unsigned int p = 0; p < num_procs; p++) {
    pc = pcs[p];
    result = pc->get_sent_ents(is_iface, bridge_dim, ghost_dim, num_layers,
                               sent_ents[p], allsent[p], entprocs[p]);
    RRAI(pc->get_moab(), "get_sent_ents failed.");
  
    //===========================================
    // pack entities into buffers
    //===========================================

    for (ind = 0; ind < pc->buffProcs.size(); ind++) {
      
        // buff_ptr points to the END (one past last occupied byte) of buffer
      buff_ptr = &pc->ownerSBuffs[ind][0];

        // entities
      result = pc->pack_entities(sent_ents[p][ind], pc->ownerSBuffs[ind], buff_ptr,
                                 store_remote_handles, pc->buffProcs[ind], is_iface,
                                 &entprocs[p], &allsent[p]); 
      RRAI(pc->get_moab(), "Packing entities failed.");
    }
  }

    //===========================================
    // receive/unpack new entities
    //===========================================
    // number of incoming messages for ghosts is the number of procs we 
    // communicate with; for iface, it's the number of those with lower rank
  std::vector<std::vector<MBEntityHandle> > L1hloc[MAX_SHARING_PROCS], L1hrem[MAX_SHARING_PROCS];
  std::vector<std::vector<int> > L1p[MAX_SHARING_PROCS];
  std::vector<MBEntityHandle> L2hloc[MAX_SHARING_PROCS], L2hrem[MAX_SHARING_PROCS];
  std::vector<unsigned int> L2p[MAX_SHARING_PROCS];
  MBRange new_ents[MAX_SHARING_PROCS];
  
  for (unsigned int p = 0; p < num_procs; p++) {
    L1hloc[p].resize(pcs[p]->buffProcs.size());
    L1hrem[p].resize(pcs[p]->buffProcs.size());
    L1p[p].resize(pcs[p]->buffProcs.size());
  }
  
  for (unsigned int p = 0; p < num_procs; p++) {
  
    MBParallelComm *pc = pcs[p];
    
    for (ind = 0; ind < pc->buffProcs.size(); ind++) {
        // incoming ghost entities; unpack; returns entities received
        // both from sending proc and from owning proc (which may be different)

        // buffer could be empty, which means there isn't any message to
        // unpack (due to this comm proc getting added as a result of indirect
        // communication); just skip this unpack
      if (pc->ownerSBuffs[ind].empty()) continue;

      unsigned int to_p = pc->buffProcs[ind];
      unsigned char *buff_ptr = &pc->ownerSBuffs[ind][0];
      result = pcs[to_p]->unpack_entities(buff_ptr,
                                          store_remote_handles, ind, is_iface,
                                          L1hloc[to_p], L1hrem[to_p], L1p[to_p], L2hloc[to_p], 
                                          L2hrem[to_p], L2p[to_p], new_ents[to_p]);
      RRAI(pc->get_moab(), "Failed to unpack entities.");
    }
  }

  if (is_iface) {
#ifdef NDEBUG
    for (unsigned int p = 0; p < num_procs; p++) {
      result = pcs[p]->check_sent_ents(allsent[p]);
      RRAI(pcs[p]->get_moab(), "Failed check on shared entities.");
      result = pcs[p]->check_all_shared_handles();
      RRAI(pcs[p]->get_moab(), "Failed check on all shared handles.");
    }
#endif
    return MB_SUCCESS;
  }
  
      //===========================================
      // send local handles for new ghosts to owner, then add
      // those to ghost list for that owner
      //===========================================
  std::vector<unsigned int>::iterator proc_it;
  for (unsigned int p = 0; p < num_procs; p++) {
    pc = pcs[p];
  
    for (ind = 0, proc_it = pc->buffProcs.begin(); 
         proc_it != pc->buffProcs.end(); proc_it++, ind++) {
        // skip if iface layer and higher-rank proc
      unsigned char *buff_ptr = &pc->ghostSBuffs[ind][0];
      result = pc->pack_remote_handles(L1hloc[p][ind], L1hrem[p][ind], L1p[p][ind], *proc_it,
                                       pc->ghostSBuffs[ind], buff_ptr);
      RRAI(pc->get_moab(), "Failed to pack remote handles.");
    }
  }
  
    //===========================================
    // process remote handles of my ghosteds
    //===========================================
  for (unsigned int p = 0; p < num_procs; p++) {
    pc = pcs[p];
  
    for (ind = 0, proc_it = pc->buffProcs.begin(); 
         proc_it != pc->buffProcs.end(); proc_it++, ind++) {
        // incoming remote handles
      unsigned int to_p = pc->buffProcs[ind];
      result = pcs[to_p]->unpack_remote_handles(p, &pc->ghostSBuffs[ind][0],
                                                L2hloc[to_p], L2hrem[to_p], L2p[to_p]);
      RRAI(pc->get_moab(), "Failed to unpack remote handles.");
    }
  }
    
#ifdef NDEBUG
  for (unsigned int p = 0; p < num_procs; p++) {
    result = pcs[p]->check_sent_ents(allsent[p]);
    RRAI(pcs[p]->get_moab(), "Failed check on shared entities.");
    result = pcs[p]->check_all_shared_handles();
    RRAI(pcs[p]->get_moab(), "Failed check on all shared handles.");
  }
#endif

  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_iface_entities(int other_proc,
                                               int dim,
                                               MBRange &iface_ents) 
{
  MBRange iface_sets;
  MBErrorCode result = MB_SUCCESS;
  
  for (MBRange::iterator rit = interfaceSets.begin(); rit != interfaceSets.end(); rit++) {
    if (-1 != other_proc && !is_iface_proc(*rit, other_proc)) continue;
    
    if (-1 == dim) result = mbImpl->get_entities_by_handle(*rit, iface_ents);
    else result = mbImpl->get_entities_by_dimension(*rit, dim, iface_ents);
    RRA(" Failed to get entities in iface set.");
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::check_sent_ents(MBRange &allsent) 
{
    // check entities to make sure there are no zero-valued remote handles
    // where they shouldn't be
  std::vector<unsigned char> pstat(allsent.size());
  MBErrorCode result = mbImpl->tag_get_data(pstatus_tag(), allsent, &pstat[0]);
  RRA("Trouble getting pstatus.");
  std::vector<MBEntityHandle> handles(allsent.size());
  result = mbImpl->tag_get_data(sharedh_tag(), allsent, &handles[0]);
  RRA("Trouble getting shared handles.");
  std::vector<int> procs(allsent.size());
  result = mbImpl->tag_get_data(sharedp_tag(), allsent, &procs[0]);
  RRA("Trouble getting shared procs.");

  MBRange bad_entities;
  
  MBRange::iterator rit;
  unsigned int i;
  MBEntityHandle dum_hs[MAX_SHARING_PROCS];
  int dum_ps[MAX_SHARING_PROCS];
  
  for (rit = allsent.begin(), i = 0; rit != allsent.end(); rit++, i++) {
    if (-1 != procs[i] && 0 == handles[i]) bad_entities.insert(*rit);
    else {
        // might be multi-shared...
      result = mbImpl->tag_get_data(sharedps_tag(), &(*rit), 1, dum_ps);
      if (MB_TAG_NOT_FOUND == result) continue;
      RRA("Trouble getting sharedps.");
      result = mbImpl->tag_get_data(sharedhs_tag(), &(*rit), 1, dum_hs);
      RRA("Trouble getting sharedhs.");

        // find first non-set proc
      int *ns_proc = std::find(dum_ps, dum_ps+MAX_SHARING_PROCS, -1);
      int num_procs = ns_proc-dum_ps;
      assert(num_procs <= MAX_SHARING_PROCS);
        // now look for zero handles in active part of dum_hs
      MBEntityHandle *ns_handle = std::find(dum_hs, dum_hs+num_procs, 0);
      int num_handles = ns_handle-dum_hs;
      assert(num_handles <= num_procs);
      if (num_handles != num_procs) bad_entities.insert(*rit);
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::pack_remote_handles(std::vector<MBEntityHandle> &L1hloc,
                                                std::vector<MBEntityHandle> &L1hrem,
                                                std::vector<int> &L1p,
                                                unsigned int to_proc,
                                                std::vector<unsigned char> &buff,
                                                unsigned char *&buff_ptr) 
{
    // 2 vectors of handles plus ints
  CHECK_BUFF_SPACE(buff, buff_ptr, ((L1p.size()+1)*sizeof(int) + 
                                    L1hloc.size()*sizeof(MBEntityHandle)));
  
    // should be in pairs of handles
  PACK_INT(buff_ptr, L1hloc.size());
  PACK_INTS(buff_ptr, &L1p[0], L1p.size());
  PACK_EH(buff_ptr, &L1hrem[0], L1hrem.size());
  PACK_EH(buff_ptr, &L1hloc[0], L1hloc.size());
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::unpack_remote_handles(unsigned int from_proc,
                                                  unsigned char *&buff_ptr,
                                                  std::vector<MBEntityHandle> &L2hloc,
                                                  std::vector<MBEntityHandle> &L2hrem,
                                                  std::vector<unsigned int> &L2p)
{
    // incoming remote handles; use to set remote handles
  int num_eh;
  UNPACK_INT(buff_ptr, num_eh);

  unsigned char *buff_proc = buff_ptr;
  buff_ptr += num_eh * sizeof(int);
  unsigned char *buff_rem = buff_ptr + num_eh * sizeof(MBEntityHandle);
  MBErrorCode result;
  MBEntityHandle hpair[2], dum_h;
  int proc;
  for (int i = 0; i < num_eh; i++) {
    UNPACK_INT(buff_proc, proc);
    UNPACK_EH(buff_ptr, hpair, 1);
    UNPACK_EH(buff_rem, hpair+1, 1);

    if (-1 != proc) {
      result = find_existing_entity(false, proc, hpair[0], 3, NULL, 0,
                                    mbImpl->type_from_handle(hpair[1]),
                                    L2hloc, L2hrem, L2p, dum_h);
      RRA("Didn't get existing entity.");
      if (dum_h) hpair[0] = dum_h;
      else hpair[0] = 0;
    }
    assert(hpair[0] && hpair[1]);
    int this_proc = from_proc;
    result = update_remote_data(hpair[0], &this_proc, hpair+1, 1, 0);
    RRA("Trouble setting remote data range on sent entities in ghost exchange.");
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_ghosted_entities(int bridge_dim,
                                                 int ghost_dim,
                                                 int to_proc, 
                                                 int num_layers,
                                                 MBRange &ghosted_ents) 
{
    // get bridge ents on interface(s)
  MBRange from_ents;
  MBErrorCode result = MB_SUCCESS;
  assert(0 < num_layers);
  for (MBRange::iterator rit = interfaceSets.begin(); rit != interfaceSets.end();
       rit++) {
    if (!is_iface_proc(*rit, to_proc)) continue;
      
      // get starting "from" entities
    if (bridge_dim == -1)
      result = mbImpl->get_entities_by_handle(*rit, from_ents);
    else
      result = mbImpl->get_entities_by_dimension(*rit, bridge_dim, from_ents);
    RRA("Couldn't get bridge ents in the set.");

      // need to get layers of bridge-adj entities
    result = MeshTopoUtil(mbImpl).get_bridge_adjacencies(from_ents, bridge_dim,
                                                         ghost_dim, ghosted_ents, 
                                                         num_layers);
    RRA("Couldn't get bridge adjacencies.");
  }
  
  result = add_verts(ghosted_ents);
  RRA("Couldn't add verts.");

  return result;
}

MBErrorCode MBParallelComm::add_verts(MBRange &sent_ents) 
{
      // get the verts adj to these entities, since we'll have to send those too

    // first check sets
  std::pair<MBRange::const_iterator, MBRange::const_iterator>
      set_range = sent_ents.equal_range(MBENTITYSET);
  MBErrorCode result = MB_SUCCESS, tmp_result;
  for (MBRange::const_iterator rit = set_range.first; rit != set_range.second; rit++) {
    tmp_result = mbImpl->get_entities_by_type(*rit, MBVERTEX, sent_ents);
    if (MB_SUCCESS != tmp_result) result = tmp_result;
  }
  RRA("Failed to get contained verts.");
  
    // now non-sets
  MBRange tmp_ents;
  std::copy(sent_ents.begin(), set_range.first, mb_range_inserter(tmp_ents));
  result = mbImpl->get_adjacencies(tmp_ents, 0, false, sent_ents,
                                   MBInterface::UNION);
  RRA("Couldn't get vertices adj to ghosted ents.");

  return result;
}


MBErrorCode MBParallelComm::exchange_tags(std::vector<MBTag> &src_tags,
                                          std::vector<MBTag> &dst_tags,
                                          MBRange &entities)
{
  MBErrorCode result;
  int success;

    // get all procs interfacing to this proc
  std::set<unsigned int> exch_procs;
  result = get_comm_procs(exch_procs);  

    // post ghost irecv's for all interface procs
    // index greqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_reqs(MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  std::vector<MPI_Status> gstatus(MAX_SHARING_PROCS);
  std::vector<unsigned int>::iterator sit;
  int ind;
  for (ind = 0, sit = buffProcs.begin(); sit != buffProcs.end(); sit++, ind++) {
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, *sit,
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    // pack and send tags from this proc to others
    // make sendReqs vector to simplify initialization
  std::fill(sendReqs, sendReqs+MAX_SHARING_PROCS, MPI_REQUEST_NULL);

    // take all shared entities if incoming list is empty
  if (entities.empty()) entities = sharedEnts;
  
  for (ind = 0, sit = buffProcs.begin(); sit != buffProcs.end(); sit++, ind++) {
    
    MBRange tag_ents = entities;
    
      // get ents shared by proc *sit
    result = filter_pstatus(tag_ents, PSTATUS_SHARED, PSTATUS_AND, *sit);
    RRA("Failed pstatus AND check.");
    
      // remote nonowned entities
    if (!tag_ents.empty()) {
      result = filter_pstatus(tag_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT);
      RRA("Failed pstatus NOT check.");
    }
    
      // pack-send; this also posts receives if store_remote_handles is true
    std::vector<MBRange> tag_ranges;
    for (std::vector<MBTag>::iterator vit = src_tags.begin(); vit != src_tags.end(); vit++) {
      const void* ptr;
      int size;
      if (tagServer->get_default_data_ref( *vit, ptr, size ) != MB_SUCCESS) {
        MBRange tagged_ents;
        tagServer->get_entities( *vit, tagged_ents );
        tag_ranges.push_back(tag_ents.intersect(tagged_ents));
      } 
      else {
        tag_ranges.push_back(tag_ents);
      }
    }
    
      // pack the data
    unsigned char *buff_ptr = &ownerSBuffs[ind][0];
    result = pack_tags(tag_ents,
                       src_tags, dst_tags, tag_ranges, 
                       ownerSBuffs[ind], buff_ptr, true, *sit);
    RRA("Failed to count buffer in pack_send_tag.");

      // now send it
    result = send_buffer(*sit, &ownerSBuffs[ind][0], 
                         buff_ptr-&ownerSBuffs[ind][0], 
                         MB_MESG_TAGS, sendReqs[ind]);
    RRA("Failed to send buffer.");
                         
  }
  
    // receive/unpack tags
  int num_incoming = exch_procs.size();
  
  while (num_incoming) {
    int ind;
    MPI_Status status;
    success = MPI_Waitany(MAX_SHARING_PROCS, &recv_reqs[0], &ind, &status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
    int new_size;
    unsigned char *buff_ptr;
    MBRange dum_range;
    
      // branch on message type
    switch (status.MPI_TAG) {
      case MB_MESG_SIZE:
          // incoming message just has size; resize buffer and re-call recv,
          // then re-increment incoming count
        assert(ind < MAX_SHARING_PROCS);
        new_size = *((int*)&ghostRBuffs[ind][0]);
        assert(0 > new_size);
        result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind],
                                MB_MESG_TAGS);
        RRA("Failed to resize recv buffer.");
        num_incoming++;
        break;
      case MB_MESG_TAGS:
          // incoming ghost entities; process
          buff_ptr = &ghostRBuffs[ind][0];
          result = unpack_tags(buff_ptr, dum_range, true,
                               buffProcs[ind]);
        RRA("Failed to recv-unpack-tag message.");
        break;
      default:
        result = MB_FAILURE;
        RRA("Failed to get message of correct type in exch_tags.");
        break;
    }
  }
  
    // ok, now wait
  MPI_Status status[MAX_SHARING_PROCS];
  success = MPI_Waitall(MAX_SHARING_PROCS, &sendReqs[0], status);
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("Failure in waitall in tag exchange.");
  }
  
  return MB_SUCCESS;
}

/*
MBErrorCode MBParallelComm::exchange_tags( MBTag src_tag, 
                                           MBTag dst_tag, 
                                           const MBRange& entities )
{
  MBErrorCode result;
  int success;

    // get all procs interfacing to this proc
  std::set<unsigned int> exch_procs;
  result = get_comm_procs(exch_procs);  

    // post ghost irecv's for all interface procs
    // index greqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_reqs(MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  std::vector<MPI_Status> gstatus(MAX_SHARING_PROCS);
  std::vector<unsigned int>::iterator sit;
  int ind;
  for (ind = 0, sit = buffProcs.begin(); sit != buffProcs.end(); sit++, ind++) {
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, *sit,
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    // figure out which entities are shared with which processors
  std::map<int,MBRange> proc_ents;
  int other_procs[MAX_SHARING_PROCS], num_sharing;
  for (MBRange::const_iterator i = entities.begin(); i != entities.end(); ++i) {
    int owner;
    result = get_owner( *i, owner );
    RRA("Failed to get entity owner.");

      // only send entities that this proc owns
    if ((unsigned)owner != proc_config().proc_rank()) 
      continue;
    
    result = get_sharing_parts( *i, other_procs, num_sharing );
    RRA("Failed to get procs sharing entity.");
    if (num_sharing == 0) // keep track of non-shared entities for later
      proc_ents[proc_config().proc_rank()].insert( *i );
    for (int j = 0; j < num_sharing; ++j)
      proc_ents[other_procs[j]].insert( *i );
  }
  
    // pack and send tags from this proc to others
    // make sendReqs vector to simplify initialization
  std::fill(sendReqs, sendReqs+MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  std::map<unsigned int,MBRange>::const_iterator mit;
  
  for (ind = 0, sit = buffProcs.begin(); sit != buffProcs.end(); sit++, ind++) {
    
      // count first
      // buffer needs to begin with the number of tags (one)
    int buff_size = sizeof(int);
    result = packed_tag_size( src_tag, proc_ents[*sit], buff_size );
    RRA("Failed to count buffer in pack_send_tag.");

    unsigned char *buff_ptr = &ownerSBuffs[ind][0];
    CHECK_BUFF_SPACE(ownerSBuffs[ind], buff_ptr, buff_size);
    PACK_INT( buff_ptr, 1 ); // number of tags
    result = pack_tag( src_tag, dst_tag, proc_ents[*sit], proc_ents[*sit],
                       ownerSBuffs[ind], buff_ptr, true, *sit );
    RRA("Failed to pack buffer in pack_send_tag.");

      // if the message is large, send a first message to tell how large
    if (INITIAL_BUFF_SIZE < buff_size) {
      int tmp_buff_size = -buff_size;
      int success = MPI_Send(&tmp_buff_size, sizeof(int), MPI_UNSIGNED_CHAR, 
                             *sit, MB_MESG_SIZE, procConfig.proc_comm());
      if (success != MPI_SUCCESS) return MB_FAILURE;
    }
    
      // send the buffer
    success = MPI_Isend(&ownerSBuffs[ind][0], buff_size, MPI_UNSIGNED_CHAR, *sit, 
                        MB_MESG_TAGS, procConfig.proc_comm(), &sendReqs[ind]);
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }
  
    // receive/unpack tags
  int num_incoming = exch_procs.size();
  
  while (num_incoming) {
    int ind;
    MPI_Status status;
    success = MPI_Waitany(MAX_SHARING_PROCS, &recv_reqs[0], &ind, &status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
    int new_size;
    unsigned char *buff_ptr;
    MBRange dum_range;
    
      // branch on message type
    switch (status.MPI_TAG) {
      case MB_MESG_SIZE:
          // incoming message just has size; resize buffer and re-call recv,
          // then re-increment incoming count
        assert(ind < MAX_SHARING_PROCS);
        new_size = *((int*)&ghostRBuffs[ind][0]);
        assert(0 > new_size);
        result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind],
                                MB_MESG_TAGS);
        RRA("Failed to resize recv buffer.");
        num_incoming++;
        break;
      case MB_MESG_TAGS:
          // incoming ghost entities; process
          buff_ptr = &ghostRBuffs[ind][0];
          result = unpack_tags(buff_ptr, dum_range, true,
                               buffProcs[ind]);
        RRA("Failed to recv-unpack-tag message.");
        break;
      default:
        result = MB_FAILURE;
        RRA("Failed to get message of correct type in exch_tags.");
        break;
    }
  }
  
    // ok, now wait
  MPI_Status status[MAX_SHARING_PROCS];
  success = MPI_Waitall(MAX_SHARING_PROCS, &sendReqs[0], status);
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("Failure in waitall in tag exchange.");
  }
  
    // if src and destination tags aren't the same, need to copy 
    // values for local entities
  if (src_tag != dst_tag) {
    const MBRange& myents = proc_ents[proc_config().proc_rank()];
    std::vector<const void*> data_ptrs(myents.size());
    std::vector<int> data_sizes(myents.size());
    result = get_moab()->tag_get_data( src_tag, myents, &data_ptrs[0], &data_sizes[0] );
    RRA("Failure to get pointers to local data.");
    result = get_moab()->tag_set_data( dst_tag, myents, &data_ptrs[0], &data_sizes[0] );
    RRA("Failure to get pointers to local data.");
  }  
  
  return MB_SUCCESS;
}
*/

MBErrorCode MBParallelComm::update_shared_mesh()
{
  MBErrorCode result;
  int success;

    // ,,,
    /*

    // get all procs interfacing to this proc
  std::set<unsigned int> iface_procs;
  result = get_interface_procs(iface_procs);
  RRA("Failed to get iface sets, procs");

    // post ghost irecv's for all interface procs
    // index greqs the same as buffer/sharing procs indices
  std::vector<MPI_Request> recv_reqs(2*MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  std::vector<MPI_Status> gstatus(MAX_SHARING_PROCS);
  std::set<unsigned int>::iterator sit;
  for (sit = iface_procs.begin(); sit != iface_procs.end(); sit++) {
    int ind = get_buffers(*sit);
    success = MPI_Irecv(&ghostRBuffs[ind][0], ghostRBuffs[ind].size(), 
                        MPI_UNSIGNED_CHAR, *sit,
                        MB_MESG_ANY, procConfig.proc_comm(), 
                        &recv_reqs[ind]);
    if (success != MPI_SUCCESS) {
      result = MB_FAILURE;
      RRA("Failed to post irecv in ghost exchange.");
    }
  }
  
    // pack and send vertex coordinates from this proc to others
    // make sendReqs vector to simplify initialization
  std::fill(sendReqs, sendReqs+2*MAX_SHARING_PROCS, MPI_REQUEST_NULL);
  MBRange recd_ents[MAX_SHARING_PROCS];
  
  for (sit = iface_procs.begin(); sit != iface_procs.end(); sit++) {
    int ind = get_buffers(*sit);
    
    MBRange vertices;
    for (MBRange::iterator rit = interfaceSets.begin(); rit != interfaceSets.end();
         rit++) {
      if (!is_iface_proc(*rit, *sit)) 
        continue;
      
      result = mbImpl->get_entities_by_type( *rit, MBVERTEX, vertices );
      RRA("Bad interface set.");
    }
    std::map<unsigned int,MBRange>::iterator ghosted = ghostedEnts.find(*sit);
    if (ghosted != ghostedEnts.end()) {
      MBRange::iterator e = ghosted->second.upper_bound(MBVERTEX);
      vertices.merge( ghosted->second.begin(), e );
    }

      // pack-send; this also posts receives if store_remote_handles is true
    MBRange sent;
    result = pack_send_entities(*sit, vertices, false, false, 
                                false, true,
                                ownerSBuffs[ind], ownerRBuffs[MAX_SHARING_PROCS+ind], 
                                sendReqs[ind], recv_reqs[MAX_SHARING_PROCS+ind], 
                                sent);
    RRA("Failed to pack-send in mesh update exchange.");
  }
  
    // receive/unpack entities
    // number of incoming messages depends on whether we're getting back
    // remote handles
  int num_incoming = iface_procs.size();
  
  while (num_incoming) {
    int ind;
    MPI_Status status;
    success = MPI_Waitany(2*MAX_SHARING_PROCS, &recv_reqs[0], &ind, &status);
    if (MPI_SUCCESS != success) {
      result = MB_FAILURE;
      RRA("Failed in waitany in ghost exchange.");
    }
    
      // ok, received something; decrement incoming counter
    num_incoming--;
    
    std::vector<MBEntityHandle> remote_handles_v, sent_ents_tmp;
    MBRange remote_handles_r;
    int new_size;
    
      // branch on message type
    switch (status.MPI_TAG) {
      case MB_MESG_SIZE:
          // incoming message just has size; resize buffer and re-call recv,
          // then re-increment incoming count
        assert(ind < MAX_SHARING_PROCS);
        new_size = *((int*)&ghostRBuffs[ind][0]);
        assert(0 > new_size);
        result = recv_size_buff(buffProcs[ind], ghostRBuffs[ind], recv_reqs[ind],
                                MB_MESG_ENTS);
        RRA("Failed to resize recv buffer.");
        num_incoming++;
        break;
      case MB_MESG_ENTS:
          // incoming ghost entities; process
        result = recv_unpack_entities(buffProcs[ind], true,
                                      false, 
                                      ghostRBuffs[ind], ghostSBuffs[ind], 
                                      sendReqs[ind], recd_ents[ind]);
        RRA("Failed to recv-unpack message.");
        break;
    }
  }
  
    // ok, now wait if requested
  MPI_Status status[2*MAX_SHARING_PROCS];
  success = MPI_Waitall(2*MAX_SHARING_PROCS, &sendReqs[0], status);
  if (MPI_SUCCESS != success) {
    result = MB_FAILURE;
    RRA("Failure in waitall in ghost exchange.");
  }
  
  return MB_SUCCESS;
}
MBErrorCode MBParallelComm::update_iface_sets(MBRange &sent_ents,
                                              std::vector<MBEntityHandle> &remote_handles, 
                                              int from_proc) 
{
  std::vector<MBEntityHandle>::iterator remote_it = remote_handles.begin();
  MBRange::iterator sent_it = sent_ents.begin();
  MBRange ents_to_remove;
  for (; sent_it != sent_ents.end(); sent_it++, remote_it++) {
    if (!*remote_it) ents_to_remove.insert(*sent_it);
  }
  
  for (MBRange::iterator set_it = interfaceSets.begin(); set_it != interfaceSets.end(); set_it++) {
    if (!is_iface_proc(*set_it, from_proc)) continue;
    MBErrorCode result = mbImpl->remove_entities(*set_it, ents_to_remove);
    RRA("Couldn't remove entities from iface set in update_iface_sets.");
  }

*/
  
  return MB_SUCCESS;
}

  //! return sharedp tag
MBTag MBParallelComm::sharedp_tag()
{
  if (!sharedpTag) {
    int def_val = -1;
    MBErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_PROC_TAG_NAME, 
                                            sizeof(int), 
                                            MB_TAG_DENSE,
                                            MB_TYPE_INTEGER, sharedpTag, 
                                            &def_val, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
      return 0;
  }
  
  return sharedpTag;
}

  //! return sharedps tag
MBTag MBParallelComm::sharedps_tag()
{
  if (!sharedpsTag) {
    MBErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_PROCS_TAG_NAME, 
                                            MAX_SHARING_PROCS*sizeof(int), 
                                            MB_TAG_SPARSE,
                                            MB_TYPE_INTEGER, sharedpsTag, NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
      return 0;
  }
  
  return sharedpsTag;
}
  
  //! return sharedh tag
MBTag MBParallelComm::sharedh_tag()
{
  if (!sharedhTag) {
    MBEntityHandle def_val = 0;
    MBErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_HANDLE_TAG_NAME, 
                                            sizeof(MBEntityHandle), 
                                            MB_TAG_DENSE,
                                            MB_TYPE_HANDLE, sharedhTag, 
                                            &def_val, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return sharedhTag;
}
  
  //! return sharedhs tag
MBTag MBParallelComm::sharedhs_tag()
{  
  if (!sharedhsTag) {
    MBErrorCode result = mbImpl->tag_create(PARALLEL_SHARED_HANDLES_TAG_NAME, 
                                            MAX_SHARING_PROCS*sizeof(MBEntityHandle), 
                                            MB_TAG_SPARSE,
                                            MB_TYPE_INTEGER, sharedhsTag, NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) 
      return 0;
  }

  return sharedhsTag;
}
  
  //! return pstatus tag
MBTag MBParallelComm::pstatus_tag()
{  
  if (!pstatusTag) {
    unsigned char tmp_pstatus = 0;
    MBErrorCode result = mbImpl->tag_create(PARALLEL_STATUS_TAG_NAME, 
                                            sizeof(unsigned char),
                                            MB_TAG_DENSE,
                                            MB_TYPE_OPAQUE, pstatusTag, 
                                            &tmp_pstatus, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return pstatusTag;
}
  
  //! return partition set tag
MBTag MBParallelComm::partition_tag()
{  
  if (!partitionTag) {
    MBErrorCode result = mbImpl->tag_create(PARALLEL_PARTITION_TAG_NAME, 
                                            sizeof(int),
                                            MB_TAG_SPARSE,
                                            MB_TYPE_INTEGER, 
                                            partitionTag, 
                                            NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return partitionTag;
}
  
  //! return pcomm tag; passes in impl 'cuz this is a static function
MBTag MBParallelComm::pcomm_tag(MBInterface *impl,
                                bool create_if_missing)
{
  MBTag this_tag = 0;
  MBErrorCode result;
  result = impl->tag_get_handle(PARALLEL_COMM_TAG_NAME, this_tag);
  if ((MB_TAG_NOT_FOUND == result || 0 == this_tag) &&
      create_if_missing) {
    result = impl->tag_create(PARALLEL_COMM_TAG_NAME, 
                              MAX_SHARING_PROCS*sizeof(MBParallelComm*),
                              MB_TAG_SPARSE,
                              MB_TYPE_OPAQUE, this_tag,
                              NULL, true);
    if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
      return 0;
  }
  
  return this_tag;
}

    //! get the indexed pcomm object from the interface
MBParallelComm *MBParallelComm::get_pcomm(MBInterface *impl, const int index) 
{
  MBTag pc_tag = pcomm_tag(impl, false);
  if (0 == pc_tag) return NULL;
  
  MBParallelComm *pc_array[MAX_SHARING_PROCS];
  MBErrorCode result = impl->tag_get_data(pc_tag, 0, 0, (void*)pc_array);
  if (MB_SUCCESS != result) return NULL;
  
  return pc_array[index];
}

MBErrorCode MBParallelComm::get_all_pcomm( MBInterface* impl, std::vector<MBParallelComm*>& list )
{
  MBTag pc_tag = pcomm_tag(impl, false);
  if (0 == pc_tag)
    return MB_TAG_NOT_FOUND;
  
  MBParallelComm *pc_array[MAX_SHARING_PROCS];
  MBErrorCode rval = impl->tag_get_data( pc_tag, 0, 0, pc_array );
  if (MB_SUCCESS != rval)
    return rval;
  
  for (int i = 0; i < MAX_SHARING_PROCS; ++i)
    if (pc_array[i])
      list.push_back( pc_array[i] );
  
  return MB_SUCCESS;
}
  

    //! get the indexed pcomm object from the interface
MBParallelComm *MBParallelComm::get_pcomm( MBInterface *impl, 
                                           MBEntityHandle prtn,
                                           const MPI_Comm* comm ) 
{
  MBErrorCode rval;
  MBParallelComm* result = 0;
  
  MBTag prtn_tag;
  rval = impl->tag_create( PARTITIONING_PCOMM_TAG_NAME, 
                           sizeof(int),
                           MB_TAG_SPARSE,
                           MB_TYPE_INTEGER,
                           prtn_tag,
                           0, true );
  if (MB_SUCCESS != rval)
    return 0;
  
  int pcomm_id;
  rval = impl->tag_get_data( prtn_tag, &prtn, 1, &pcomm_id );
  if (MB_SUCCESS == rval) {
    result= get_pcomm( impl, pcomm_id );
  }
  else if (MB_TAG_NOT_FOUND == rval && comm) {
    result = new MBParallelComm( impl, *comm, &pcomm_id );
    if (!result)
      return 0;
    result->set_partitioning( prtn );
    
    rval = impl->tag_set_data( prtn_tag, &prtn, 1, &pcomm_id );
    if (MB_SUCCESS != rval) {
      delete result;
      result = 0;
    }
  }
  
  return result;
}

MBErrorCode MBParallelComm::set_partitioning( MBEntityHandle set) 
{
  MBErrorCode rval;
  MBTag prtn_tag;
  rval = mbImpl->tag_create( PARTITIONING_PCOMM_TAG_NAME, 
                           sizeof(int),
                           MB_TAG_SPARSE,
                           MB_TYPE_INTEGER,
                           prtn_tag,
                           0, true );
  if (MB_SUCCESS != rval)
    return rval;

    // get my id
  MBParallelComm* pcomm_arr[MAX_SHARING_PROCS];
  MBTag pc_tag = pcomm_tag(mbImpl, false);
  if (0 == pc_tag) 
    return MB_FAILURE;
  MBErrorCode result = mbImpl->tag_get_data(pc_tag, 0, 0, pcomm_arr);
  if (MB_SUCCESS != result) 
    return MB_FAILURE;  
  int id = std::find(pcomm_arr,pcomm_arr+MAX_SHARING_PROCS,this) - pcomm_arr;
  if (id == MAX_SHARING_PROCS)
    return MB_FAILURE;

  MBEntityHandle old = partitioningSet;
  if (old) {
    rval = mbImpl->tag_delete_data( prtn_tag, &old, 1 );
    if (MB_SUCCESS != rval)
      return rval;
    partitioningSet = 0;
  }
  
  if (!set) 
    return MB_SUCCESS;
  
  MBRange contents;
  if (old) {
    rval = mbImpl->get_entities_by_handle( old, contents );
    if (MB_SUCCESS != rval)
      return rval;
  }
  else {
    contents = partition_sets();
  }

  rval = mbImpl->add_entities( set, contents );
  if (MB_SUCCESS != rval)
    return rval;
  
    // store pcomm id on new partition set
  rval = mbImpl->tag_set_data( prtn_tag, &set, 1, &id );
  if (MB_SUCCESS != rval)
    return rval;
  
  partitioningSet = set;
  return MB_SUCCESS;
}
  

  //! return all the entities in parts owned locally
MBErrorCode MBParallelComm::get_part_entities(MBRange &ents, int dim) 
{
  MBErrorCode result;
  
  for (MBRange::iterator rit = partitionSets.begin(); 
       rit != partitionSets.end(); rit++) {
    MBRange tmp_ents;
    if (-1 == dim) 
      result = mbImpl->get_entities_by_handle(*rit, tmp_ents, true);
    else
      result = mbImpl->get_entities_by_dimension(*rit, dim, tmp_ents, true);

    if (MB_SUCCESS != result) return result;
    ents.merge(tmp_ents);
  }
  
  return MB_SUCCESS;
}

  /** \brief Return the rank of the entity owner
   */
MBErrorCode MBParallelComm::get_owner_handle(MBEntityHandle entity,
                                             int &owner,
                                             MBEntityHandle &handle) 
{
  unsigned char pstat;
  int sharing_procs[MAX_SHARING_PROCS];
  MBEntityHandle sharing_handles[MAX_SHARING_PROCS];

  MBErrorCode result = mbImpl->tag_get_data(pstatus_tag(), &entity, 1,
                                            &pstat);
  if (!(pstat & PSTATUS_NOT_OWNED)) {
    owner = proc_config().proc_rank();
    handle = entity;
  }
  
  else if (pstat & PSTATUS_MULTISHARED) {
    result = mbImpl->tag_get_data(sharedps_tag(), &entity, 1,
                                  sharing_procs);
    owner = sharing_procs[0];
    result = mbImpl->tag_get_data(sharedhs_tag(), &entity, 1,
                                  sharing_handles);
    handle = sharing_handles[0];
  }
  else if (pstat & PSTATUS_SHARED) {
    result = mbImpl->tag_get_data(sharedp_tag(), &entity, 1,
                                  sharing_procs);
    RRA(" ");
    owner = sharing_procs[0];
    result = mbImpl->tag_get_data(sharedh_tag(), &entity, 1,
                                  sharing_handles);
    handle = sharing_handles[0];
  }
  else {
    owner = -1;
    handle = 0;
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_global_part_count( int& count_out ) const
{
  count_out = globalPartCount;
  return count_out < 0 ? MB_FAILURE : MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_part_owner( int part_id, int& owner ) const
{
  // FIXME: assumes only 1 local part
  owner = part_id;
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_part_id( MBEntityHandle /*part*/, int& id_out ) const
{
  // FIXME: assumes only 1 local part
  id_out = proc_config().proc_rank();
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_part_handle( int id, MBEntityHandle& handle_out ) const
{
  // FIXME: assumes only 1 local part
  if ((unsigned)id != proc_config().proc_rank())
    return MB_ENTITY_NOT_FOUND;
  handle_out = partition_sets().front();
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::create_part( MBEntityHandle& set_out )
{
    // mark as invalid so we know that it needs to be updated
  globalPartCount = -1;
  
    // create set representing part
  MBErrorCode rval = mbImpl->create_meshset( MESHSET_SET, set_out );
  if (MB_SUCCESS != rval)
    return rval;
  
    // set tag on set
    // FIXME: need to assign valid global id
  int val = 0;
  rval = mbImpl->tag_set_data( part_tag(), &set_out, 1, &val );
  if (MB_SUCCESS != rval) {
    mbImpl->delete_entities( &set_out, 1 );
    return rval;
  }
  
  if (get_partitioning()) {
    rval = mbImpl->add_entities( get_partitioning(), &set_out, 1 );
    if (MB_SUCCESS != rval) {
      mbImpl->delete_entities( &set_out, 1 );
      return rval;
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::destroy_part( MBEntityHandle part_id )
{
    // mark as invalid so we know that it needs to be updated
  globalPartCount = -1;
  
  MBErrorCode rval;
  if (get_partitioning()) {
    rval = mbImpl->remove_entities( get_partitioning(), &part_id, 1 );
    if (MB_SUCCESS != rval)
      return rval;
  }
  return mbImpl->delete_entities( &part_id, 1 );
}

MBErrorCode MBParallelComm::collective_sync_partition()
{
  int count = partition_sets().size();
  globalPartCount = 0;
  int err = MPI_Allreduce( &count, &globalPartCount, 1, MPI_INT, MPI_SUM, 
                           proc_config().proc_comm() );
  return err ? MB_FAILURE : MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_part_neighbor_ids( MBEntityHandle part,
                                                   int neighbors_out[MAX_SHARING_PROCS],
                                                   int& num_neighbors_out )
{
  MBErrorCode rval;
  MBRange iface;
  rval = get_interface_sets( part, iface );
  if (MB_SUCCESS != rval)
    return rval;
  
  num_neighbors_out = 0;
  int n, j = 0;
  int tmp[MAX_SHARING_PROCS], curr[MAX_SHARING_PROCS];
  int *parts[2] = { neighbors_out, tmp };
  for (MBRange::iterator i = iface.begin(); i != iface.end(); ++i) {
    unsigned char pstat;
    rval = get_sharing_data( *i, curr, NULL, pstat, n);
    if (MB_SUCCESS != rval)
      return rval;
    std::sort( curr, curr+n );
    assert( num_neighbors_out < MAX_SHARING_PROCS );
    int* k = std::set_union( parts[j], parts[j]+num_neighbors_out,
                             curr, curr + n, parts[1-j] );
    j = 1-j;
    num_neighbors_out = k - parts[j];
  }
  if (parts[j] != neighbors_out)
    std::copy( parts[j], parts[j]+num_neighbors_out, neighbors_out );
    
    
    // remove input part from list
  int id;
  rval = get_part_id( part, id );
  if (MB_SUCCESS == rval) 
    num_neighbors_out = std::remove( neighbors_out, neighbors_out+num_neighbors_out, id ) - neighbors_out;
  return rval;
}

MBErrorCode MBParallelComm::get_interface_sets( MBEntityHandle ,
                                                MBRange& iface_sets_out,
                                                int* adj_part_id )
{
    // FIXME : assumes one part per processor.
    // Need to store part iface sets as children to implement
    // this correctly.
  iface_sets_out = interface_sets();

  if (adj_part_id) {
    int part_ids[MAX_SHARING_PROCS], num_parts;
    MBRange::iterator i = iface_sets_out.begin();
    while (i != iface_sets_out.end()) {
      unsigned char pstat;
      MBErrorCode rval = get_sharing_data( *i, part_ids, NULL, pstat, num_parts );
      if (MB_SUCCESS != rval)
        return rval;
      
      if (std::find(part_ids, part_ids+num_parts, *adj_part_id) - part_ids != num_parts)
        ++i;
      else
        i = iface_sets_out.erase( i );
    }
  }
    
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_owning_part( MBEntityHandle handle,
                                             int& owning_part_id,
                                             MBEntityHandle* remote_handle )
{

  // FIXME : assumes one part per proc, and therefore part_id == rank
  
    // If entity is not shared, then we're the owner.
  unsigned char pstat;
  MBErrorCode result = mbImpl->tag_get_data(pstatus_tag(), &handle, 1,
                                            &pstat);
  if (!(pstat & PSTATUS_NOT_OWNED)) {
    owning_part_id = proc_config().proc_rank();
    if (remote_handle)
      *remote_handle = handle;
    return MB_SUCCESS;
  }
  
    // If entity is shared with one other proc, then
    // sharedp_tag will contain a positive value.
  result = mbImpl->tag_get_data( sharedp_tag(), &handle, 1, &owning_part_id );
  if (MB_SUCCESS != result)
    return result;
  if (owning_part_id != -1) {
      // done?
    if (!remote_handle)
      return MB_SUCCESS;
      
      // get handles on remote processors (and this one)
    return mbImpl->tag_get_data( sharedh_tag(), &handle, 1, remote_handle );
  }
  
    // If here, then the entity is shared with at least two other processors.
    // Get the list from the sharedps_tag
  const void* part_id_list = 0;
  result = mbImpl->tag_get_data( sharedps_tag(), &handle, 1, &part_id_list );
  if (MB_SUCCESS != result)
    return result;
  owning_part_id = ((const int*)part_id_list)[0];
 
    // done?
  if (!remote_handle)
    return MB_SUCCESS;
  
    // get remote handles
  const void* handle_list = 0;
  result = mbImpl->tag_get_data( sharedhs_tag(), &handle, 1, &handle_list );
  if (MB_SUCCESS != result)
    return result;
  
  *remote_handle = ((const MBEntityHandle*)handle_list)[0];
  return MB_SUCCESS;
}    

MBErrorCode MBParallelComm::exchange_all_shared_handles(std::vector<std::vector<SharedEntityData> > &result)
{
  MBErrorCode rval;
  int ierr;
  const int tag = 0x4A41534E;
  const MPI_Comm comm = procConfig.proc_comm();
  const int num_proc = buffProcs.size();
  std::vector<MPI_Request> send_req(num_proc), recv_req(num_proc);
  const std::vector<int> procs( buffProcs.begin(), buffProcs.end() );
  
    // get all shared entities
  MBRange all_shared, dum_range;
  rval = mbImpl->get_entities_by_handle(0, all_shared);
  if (MB_SUCCESS != rval)
    return rval;
  rval = get_pstatus_entities(-1, 0, dum_range);
  if (MB_SUCCESS != rval)
    return rval;
  all_shared = all_shared.subtract(dum_range);
  all_shared.erase(all_shared.upper_bound(MBPOLYHEDRON), all_shared.end());
  assert(sharedEnts == all_shared);

    // build up send buffers
  std::vector<std::vector<SharedEntityData> > send_data(buffProcs.size());
  int ent_procs[MAX_SHARING_PROCS];
  MBEntityHandle handles[MAX_SHARING_PROCS];
  int num_sharing;
  SharedEntityData tmp;
  for (MBRange::iterator i = all_shared.begin(); i != all_shared.end(); ++i) {
    tmp.remote = *i; // swap local/remote so they're correct on the remote proc.
    rval = get_owner( *i, tmp.owner );
    if (MB_SUCCESS != rval)
      return rval;

    unsigned char pstat;
    rval = get_sharing_data( *i, ent_procs, handles, pstat, num_sharing );
    for (int j = 0; j < num_sharing; ++j) {
      if (ent_procs[j] == (int)proc_config().proc_rank())
        continue;
      tmp.local = handles[j];
      int ind = get_buffers(ent_procs[j]);
      assert(-1 != ind);
      send_data[ind].push_back( tmp );
    }
  }

    // set up to receive sizes
  std::vector<int> sizes_send(num_proc), sizes_recv(num_proc);
  for (int i = 0; i < num_proc; ++i) {
    ierr = MPI_Irecv( &sizes_recv[i], 1, MPI_INT, procs[i], tag, comm, &recv_req[i] );
    if (ierr) 
      return MB_FILE_WRITE_ERROR;
  }
  
    // send sizes
  for (int i = 0; i < num_proc; ++i) {
    sizes_send[i] = send_data[i].size();
    ierr = MPI_Isend( &sizes_send[i], 1, MPI_INT, buffProcs[i], tag, comm, &send_req[i] );
    if (ierr) 
      return MB_FILE_WRITE_ERROR;
  }
  
    // receive sizes
  std::vector<MPI_Status> stat(num_proc);
  ierr = MPI_Waitall( num_proc, &recv_req[0], &stat[0] );
  if (ierr)
    return MB_FILE_WRITE_ERROR;
  
    // wait until all sizes are sent (clean up pending req's)
  ierr = MPI_Waitall( num_proc, &send_req[0], &stat[0] );
  if (ierr)
    return MB_FILE_WRITE_ERROR;
  
    // set up to receive data
  for (int i = 0; i < num_proc; ++i) {
    result[i].resize( sizes_recv[i] );
    ierr = MPI_Irecv( &result[i][0], 
                      sizeof(SharedEntityData)*sizes_recv[i], 
                      MPI_UNSIGNED_CHAR, 
                      buffProcs[i], tag, comm, &recv_req[i] );
    if (ierr) 
      return MB_FILE_WRITE_ERROR;
  }
  
    // send data
  for (int i = 0; i < num_proc; ++i) {
    ierr = MPI_Isend( &send_data[i][0], 
                      sizeof(SharedEntityData)*sizes_send[i], 
                      MPI_UNSIGNED_CHAR, 
                      buffProcs[i], tag, comm, &send_req[i] );
    if (ierr) 
      return MB_FILE_WRITE_ERROR;
  }
  
    // receive data
  ierr = MPI_Waitall( num_proc, &recv_req[0], &stat[0] );
  if (ierr)
    return MB_FILE_WRITE_ERROR;
  
    // wait until everything is sent to release send buffers
  ierr = MPI_Waitall( num_proc, &send_req[0], &stat[0] );
  if (ierr)
    return MB_FILE_WRITE_ERROR;
  
  return MB_SUCCESS;
}

MBErrorCode MBParallelComm::check_all_shared_handles() 
{
    // get all shared ent data from other procs
  std::vector<std::vector<SharedEntityData> > shents(buffProcs.size());
  MBErrorCode result;
  result = exchange_all_shared_handles(shents);
  if (MB_SUCCESS != result)
    return result;
  else if (shents.empty())
    return MB_SUCCESS;
  
    // now check against what I think data should be
    // get all shared entities
  MBRange all_shared, dum_range;
  result = mbImpl->get_entities_by_handle(0, all_shared);
  if (MB_SUCCESS != result)
    return result;
  result = get_pstatus_entities(-1, 0, dum_range);
  if (MB_SUCCESS != result)
    return result;
  all_shared = all_shared.subtract(dum_range);
  all_shared.erase(all_shared.upper_bound(MBPOLYHEDRON), all_shared.end());

  MBRange bad_ents, local_shared;
  std::vector<SharedEntityData>::iterator vit;
  for (unsigned int i = 0; i < shents.size(); i++) {
    int other_proc = buffProcs[i];
    local_shared = all_shared;
    for (vit = shents[i].begin(); vit != shents[i].end(); vit++) {
      MBEntityHandle localh = vit->local, remoteh = vit->remote, dumh;
      local_shared.erase(localh);
      result = get_remote_handles(true, &localh, &dumh, 1, other_proc, dum_range);
      if (MB_SUCCESS != result || dumh != remoteh) 
        bad_ents.insert(localh);
    }
  }
  
  if (!bad_ents.empty() || !local_shared.empty()) return MB_FAILURE;
  else return MB_SUCCESS;
}

MBErrorCode MBParallelComm::get_shared_entities(int other_proc,
                                                MBRange &shared_ents,
                                                int dim,
                                                const bool iface,
                                                const bool owned_filter) 
{
  shared_ents.clear();
  MBErrorCode result = MB_SUCCESS;
  
    // dimension
  if (-1 != dim) {
    MBDimensionPair dp = MBCN::TypeDimensionMap[dim];
    MBRange dum_range;
    shared_ents.merge(sharedEnts.lower_bound(dp.first), 
                      sharedEnts.upper_bound(dp.second));
  }
  else shared_ents = sharedEnts;

    // filter by iface
  if (iface) {
    result = filter_pstatus(shared_ents, PSTATUS_INTERFACE, PSTATUS_AND);
    RRA("");
  }
  
    // filter by owned
  if (owned_filter) {
    result = filter_pstatus(shared_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT);
    RRA("");
  }

    // filter by proc
  if (-1 != other_proc) {
    result = filter_pstatus(shared_ents, PSTATUS_SHARED, PSTATUS_AND, other_proc);
    RRA("");
  }
  
  return result;
}

#ifdef TEST_PARALLELCOMM

#include <iostream>

#include "MBCore.hpp"
#include "MBParallelComm.hpp"
#include "MBRange.hpp"

#define PM {std::cerr << "Test failed; error message:" << std::endl;\
          std::string errmsg; \
          dynamic_cast<MBCore*>(my_impl)->get_last_error(errmsg); \
          std::cerr << errmsg << std::endl;\
          return 1;}

int main(int argc, char* argv[])
{

    // Check command line arg
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " <mesh_file_name>" << std::endl;
    exit(1);
  }

  const char* file = argv[1];
  MBCore *my_impl = new MBCore(0, 2);
  MBInterface* mbImpl = my_impl;

    // create a communicator class, which will start mpi too
  MBParallelComm pcomm(mbImpl, my_impl->tag_server(), my_impl->sequence_manager());
  MBErrorCode result;

    // load the mesh
  result = mbImpl->load_mesh(file, 0, 0);
  if (MB_SUCCESS != result) return result;

    // get the mesh
  MBRange all_mesh, whole_range;
  result = mbImpl->get_entities_by_dimension(0, 3, all_mesh);
  if (MB_SUCCESS != result) return result;
    
  int buff_size;
  result = pcomm.pack_buffer(all_mesh, false, true, true, false, whole_range, buff_size);
  PM;


    // allocate space in the buffer
  pcomm.buffer_size(buff_size);

    // pack the actual buffer
  int actual_buff_size;
  result = pcomm.pack_buffer(whole_range, false, true, false, false, all_mesh, 
                             actual_buff_size);
  PM;

    // list the entities that got packed
  std::cout << "ENTITIES PACKED:" << std::endl;
  mbImpl->list_entities(all_mesh);

    // get the buffer
  std::vector<unsigned char> tmp_buffer;
  pcomm.take_buffer(tmp_buffer);
    
    // stop and restart MOAB
  delete mbImpl;
  my_impl = new MBCore(1, 2);
  mbImpl = my_impl;
    
    // create a new communicator class, using our old buffer
  MBParallelComm pcomm2(mbImpl, my_impl->tag_server(), my_impl->sequence_manager(),
                        tmp_buffer);

    // unpack the results
  all_mesh.clear();
  result = pcomm2.unpack_buffer(all_mesh, store_remote_handles, from_proc);
  PM;
  
  std::cout << "ENTITIES UNPACKED:" << std::endl;
  mbImpl->list_entities(all_mesh);
  
  std::cout << "Success, processor " << mbImpl->proc_rank() << "." << std::endl;
  
  return 1;
}
#endif
