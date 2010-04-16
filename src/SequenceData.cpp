#include "SequenceData.hpp"
#include "TagServer.hpp"
#include "SysUtil.hpp"
#include "VarLenTag.hpp"
#include <assert.h>

namespace moab {

SequenceData::~SequenceData()
{
  for (int i = -numSequenceData; i <= (int)numTagData; ++i)
    free( arraySet[i] );
  free( arraySet - numSequenceData );
}

void* SequenceData::create_data( int index, int bytes_per_ent, const void* initial_value )
{  
  char* array;
  MALLOC(array, bytes_per_ent * size(), char*);
  if (initial_value)
    SysUtil::setmem( array, initial_value, bytes_per_ent, size() );
  else 
    memset( array, 0, bytes_per_ent * size() );
  
  arraySet[index] = array;
  return array;
}

void* SequenceData::create_sequence_data( int array_num,
                                          int bytes_per_ent,
                                          const void* initial_value )
{
  const int index = -1 - array_num;
  assert( array_num < numSequenceData );
  assert( !arraySet[index] );
  return create_data( index, bytes_per_ent, initial_value );
}


void* SequenceData::create_custom_data( int array_num, size_t total_bytes )
{
  const int index = -1 - array_num;
  assert( array_num < numSequenceData );
  assert( !arraySet[index] );

  void* array;
  MALLOC(array, total_bytes, void*);
  memset( array, 0, total_bytes );
  arraySet[index] = array;
  return array;
}

SequenceData::AdjacencyDataType* SequenceData::allocate_adjacency_data()
{
  assert( !arraySet[0] );
  const size_t s = sizeof(AdjacencyDataType*) * size();
  MALLOC(arraySet[0], s, void* );
  memset( arraySet[0], 0, s );
  return reinterpret_cast<AdjacencyDataType*>(arraySet[0]);
}

void SequenceData::increase_tag_count( unsigned amount )
{
  void** list = arraySet - numSequenceData;
  const size_t size = sizeof(void*) * (numSequenceData + numTagData + amount + 1);
  REALLOC(list, list, (size - sizeof(void*)*amount), size, void**);
//  list = (void**)realloc( list, size );
  arraySet = list + numSequenceData;
  memset( arraySet + numTagData + 1, 0, sizeof(void*) * amount );
  numTagData += amount;
}

void* SequenceData::create_tag_data( TagId tag_num,
                                     int bytes_per_ent,
                                     const void* initial_val )
{
  if (tag_num >= numTagData)
    increase_tag_count( tag_num - numTagData + 1 );
  
  assert( !arraySet[tag_num + 1] );
  return create_data( tag_num + 1, bytes_per_ent, initial_val );
}

SequenceData* SequenceData::subset( EntityHandle start,
                                    EntityHandle end,
                                    const int* sequence_data_sizes ) const
{
  return new SequenceData( this, start, end, sequence_data_sizes );
}

SequenceData::SequenceData( const SequenceData* from,
                            EntityHandle start, 
                            EntityHandle end,
                            const int* sequence_data_sizes )
  : numSequenceData( from->numSequenceData ),
    numTagData( from->numTagData ),
    startHandle( start ),
    endHandle( end )
{
  assert( start <= end );
  assert( from != 0 );
  assert( from->start_handle() <= start );
  assert( from->end_handle() >= end );

  void **array;
  MALLOC(array, sizeof(void*) * (numSequenceData + numTagData + 1), void**);
  arraySet = array + numSequenceData;
  const size_t offset = start - from->start_handle();
  const size_t count = end - start + 1;
  
  for (int i = 0; i < numSequenceData; ++i)
    copy_data_subset( -1 - i, sequence_data_sizes[i], from->get_sequence_data(i), offset, count );
  copy_data_subset( 0, sizeof(AdjacencyDataType*), from->get_adjacency_data(), offset, count );
  for (unsigned i = 1; i <= numTagData; ++i)
    arraySet[i] = 0;
}

void SequenceData::copy_data_subset( int index, 
                                     int size_per_ent, 
                                     const void* source, 
                                     size_t offset, 
                                     size_t count )
{
  if (!source)
    arraySet[index] = 0;
  else {
    MALLOC(arraySet[index], count * size_per_ent, void*);
    memcpy( arraySet[index], 
            (const char*)source + offset * size_per_ent, 
            count * size_per_ent );
  }
}

void SequenceData::move_tag_data( SequenceData* destination, TagServer* tag_server )
{
  assert( destination->start_handle() >= start_handle() );
  assert( destination->end_handle() <= end_handle() );
  const size_t offset = destination->start_handle() - start_handle();
  const size_t count = destination->size();
  if (destination->numTagData < numTagData)
    destination->increase_tag_count( numTagData - destination->numTagData );
  
  for (unsigned i = 1; i <= numTagData; ++i) {
    if (!arraySet[i])
      continue;
    
    const TagInfo* info = tag_server->get_tag_info( TAG_HANDLE_FROM_ID( i-1, MB_TAG_DENSE ) );
    if (!info)
      continue;
    
    const int tag_size = info->get_size();
    if (!destination->arraySet[i])
      MALLOC(destination->arraySet[i], count * tag_size, void*);
    memcpy( destination->arraySet[i], 
            reinterpret_cast<char*>(arraySet[i]) + offset * tag_size,
            count * tag_size );
  }
}

void SequenceData::release_tag_data( const int* tag_sizes, int num_tag_sizes )
{
  assert( num_tag_sizes >= (int)numTagData );
  for (unsigned i = 0; i < numTagData; ++i)
    release_tag_data( i, tag_sizes[i] );
}

void SequenceData::release_tag_data( TagId tag_num, int tag_size )
{
  if (tag_num < numTagData) {
    if (tag_size == MB_VARIABLE_LENGTH && arraySet[tag_num+1]) {
      VarLenTag* iter = reinterpret_cast<VarLenTag*>(arraySet[tag_num+1]);
      VarLenTag* const end = iter + size();
      for (; iter != end; ++iter)
        iter->clear();
    }
    free( arraySet[tag_num+1] );
    arraySet[tag_num+1] = 0;
  }
}

} // namespace moab
