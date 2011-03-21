#ifndef MBIMESH_HPP
#define MBIMESH_HPP

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include <vector>
#include <algorithm>

using namespace moab;

class MBiMesh
{
private:
  bool haveDeletedEntities;
  bool iCreatedInterface;
  std::vector<Tag> setHandleTags, entHandleTags;
public:
  MBiMesh(moab::Interface *mbImpl = NULL);

  virtual ~MBiMesh();
  bool have_deleted_ents( bool reset ) {
    bool result = haveDeletedEntities;
    if (reset)
      haveDeletedEntities = false;
    return result;
  }

  virtual ErrorCode delete_mesh();
  virtual ErrorCode delete_entities( const EntityHandle*, const int );
  virtual ErrorCode delete_entities( const Range& );
  int AdjTable[16];
  moab::Interface *mbImpl;
  
  inline void note_set_handle_tag( Tag );
  inline void note_ent_handle_tag( Tag );
  inline void note_tag_destroyed( Tag );
  inline bool is_set_handle_tag( Tag ) const;
  inline bool is_ent_handle_tag( Tag ) const;
};

static inline MBiMesh *mbimeshi_instance(iMesh_Instance instance) {return reinterpret_cast<MBiMesh*>(instance);}
#define MBIMESHI mbimeshi_instance(instance)
#define MOABI MBIMESHI->mbImpl

inline MBiMesh::MBiMesh(Interface *impl)
        : haveDeletedEntities(false), iCreatedInterface(false), mbImpl(impl)
{
  int tmp_table[] = {
      1, 1, 1, 1,
      1, 0, 2, 2,
      1, 2, 0, 2,
      1, 2, 2, 1
  };
  memcpy(AdjTable, tmp_table, 16*sizeof(int));

  if (!mbImpl) {
    mbImpl = new Core();
    iCreatedInterface = true;
  }
}

inline MBiMesh::~MBiMesh() 
{
  if (iCreatedInterface) delete mbImpl;
}

inline ErrorCode MBiMesh::delete_mesh() {
  haveDeletedEntities = true;
  return mbImpl->delete_mesh();
}

inline ErrorCode MBiMesh::delete_entities( const EntityHandle* a, const int n )
{
  if (n > 0)
    haveDeletedEntities = true;
  return mbImpl->delete_entities( a, n );
}

inline ErrorCode MBiMesh::delete_entities( const Range& r )
{
  if (!r.empty())
    haveDeletedEntities = true;
  return mbImpl->delete_entities( r );
}


void MBiMesh::note_set_handle_tag( Tag t )
{
  std::vector<Tag>::iterator i;
  i = std::lower_bound( entHandleTags.begin(), entHandleTags.end(), t );
  if (i != entHandleTags.end() && *i == t)
    entHandleTags.erase(i);
  i = std::lower_bound( setHandleTags.begin(), setHandleTags.end(), t );
  if (i == setHandleTags.end() || *i != t)
    setHandleTags.insert( i, t );
}

void MBiMesh::note_ent_handle_tag( Tag t )
{
  std::vector<Tag>::iterator i;
  i = std::lower_bound( setHandleTags.begin(), setHandleTags.end(), t );
  if (i != setHandleTags.end() && *i == t)
    setHandleTags.erase(i);
  i = std::lower_bound( entHandleTags.begin(), entHandleTags.end(), t );
  if (i == entHandleTags.end() || *i != t)
    entHandleTags.insert( i, t );
}

void MBiMesh::note_tag_destroyed( Tag t )
{
  std::vector<Tag>::iterator i;
  i = std::lower_bound( setHandleTags.begin(), setHandleTags.end(), t );
  if (i != setHandleTags.end() && *i == t)
    setHandleTags.erase(i);
  i = std::lower_bound( entHandleTags.begin(), entHandleTags.end(), t );
  if (i != entHandleTags.end() && *i == t)
    entHandleTags.erase(i);
}

bool MBiMesh::is_set_handle_tag( Tag t ) const
{
  return std::binary_search( setHandleTags.begin(), setHandleTags.end(), t );
}

bool MBiMesh::is_ent_handle_tag( Tag t ) const
{
  return std::binary_search( entHandleTags.begin(), entHandleTags.end(), t );
}

#endif
