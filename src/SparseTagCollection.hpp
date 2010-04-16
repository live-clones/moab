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


/**********************************************
 * Filename   :     SparseTagCollection.hpp
 *
 * Purpose    :     To store any size data with
 *                  any entity handle
 *
 * Creator    :     Clinton Stimpson
 *
 * Date       :     3 April 2002
 *
 * ********************************************/



#ifndef SPARSE_TAG_COLLECTION_HPP
#define SPARSE_TAG_COLLECTION_HPP

#ifndef IS_BUILDING_MB
#error "SparseTagCollection.hpp isn't supposed to be included into an application"
#endif

#ifdef WIN32
#pragma warning(disable : 4786)
#endif


#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#ifdef HAVE_UNORDERED_MAP
# include STRINGIFY(HAVE_UNORDERED_MAP)
#else
# include <map>
#endif
#include <vector>

#include "moab/Types.hpp"
#include "Internals.hpp"
#include "moab/Range.hpp"
#include "TagInfo.hpp"

namespace moab {

//! allocator for tag data
class SparseTagDataAllocator
{
public:
  //! constructor
  SparseTagDataAllocator(){}
  //! destructor
  ~SparseTagDataAllocator(){}
  //! allocates memory of size and returns pointer
  void* allocate(size_t data_size) { void *array; MALLOC(array, data_size, void*); return array;}
  //! frees the memory
  void destroy(void* p){ free(p); }
};


//! collection of tag data associated with entity ids
class SparseTagCollection
{
public:
  
  //! constructor takes tag data size
  SparseTagCollection(int data_size);
  
  //! destructor
  ~SparseTagCollection();
  
  //! set the tag data for an entity id
  //!\NOTE Will fail with MB_VARIABLE_DATA_LENGTH if called for 
  //!      variable-length tag.
  ErrorCode set_data(const EntityHandle entity_handle, const void* data);

  //! get the tag data for an entity id
  //!\NOTE Will fail with MB_VARIABLE_DATA_LENGTH if called for 
  //!      variable-length tag.
  ErrorCode get_data(const EntityHandle entity_handle, void* data);
  
  //! set variable-length tag data for an entity id
  //!
  //!\NOTE If called for fixed-length tag, size must be either zero or the tag size.
  //!
  //!\NOTE If called with zero size for a variable-length tag, is equivalent
  //!      to remove_data().
  ErrorCode set_data(const EntityHandle entity_handle, const void* data, int length);

  //! get the variable-length data for an entity id
  ErrorCode get_data(const EntityHandle entity_handle, const void*& data, int& length);

  //! removes the data
  ErrorCode remove_data(const EntityHandle entity_handle);

  //! get number of entities of type
  ErrorCode get_number_entities(EntityType type, int& num_entities);
  
  //! get number of entities
  unsigned long get_number_entities()
    { return mData.size(); }

  //! gets all entity handles that match a type and tag
  ErrorCode get_entities(EntityType type, Range &entities);

  //! gets all entity handles that match a tag
  ErrorCode get_entities(Range &entities) const;

  //! gets all entity handles that match a type, tag, tag_value
  ErrorCode get_entities_with_tag_value( const TagInfo& info,
                                           EntityType type, 
                                           Range &entities, 
                                           const void* tag_value,
                                           int value_size);

  //! if this collection contains this entity, return true, otherwise false
  bool contains(const EntityHandle entity) const;
  
  int tag_size() const { return mDataSize; }

protected:
  
  //! hidden constructor
  SparseTagCollection(){}
  
  //! size of the data
  int mDataSize;

  //! allocator for this collection
  SparseTagDataAllocator mAllocator;

  //! map of entity id and tag data
#ifdef HAVE_UNORDERED_MAP
  typedef UNORDERED_MAP_NS::unordered_map<EntityHandle,void*> myMapType;
#else
  typedef std::map<EntityHandle /*entity_handle*/ , void* /*data*/ > myMapType;
#endif

  myMapType mData;
};

inline bool SparseTagCollection::contains(const EntityHandle entity) const
{
  return (mData.find(entity) == mData.end() ? false : true);
}

inline ErrorCode SparseTagCollection::get_entities(Range &entities) const 
{
  for (myMapType::const_iterator mit = mData.begin(); mit != mData.end(); mit++) 
    entities.insert((*mit).first);

  return MB_SUCCESS;
}

} // namespace moab

#endif //SPARSE_TAG_COLLECTION_HPP




