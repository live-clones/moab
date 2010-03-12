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

#include "GeomTopoUtil.hpp"
#include "moab/Range.hpp"
#include "moab/MBTagConventions.hpp"
#include "moab/Interface.hpp"
#include "moab/MBCN.hpp"
#include "Internals.hpp"
#include <assert.h>
#include <iostream>

// Tag name used for saving sense of faces in volumes.
// We assume that the surface occurs in at most two volumes.
// Code will error out if more than two volumes per surface.
// The tag data is a pair of tag handles, representing the
// forward and reverse volumes, respectively.  If a surface
// is non-manifold in a single volume, the same volume will
// be listed for both the forward and reverse slots.
const char GEOM_SENSE_TAG_NAME[] = "GEOM_SENSE_2";

ErrorCode GeomTopoUtil::set_sense( EntityHandle surface,
                                     EntityHandle volume,
                                     bool forward )
{
  ErrorCode rval;
  if (!sense2Tag) {
    rval = mdbImpl->tag_create( GEOM_SENSE_TAG_NAME, 2*sizeof(EntityHandle), 
                           MB_TAG_SPARSE, MB_TYPE_HANDLE, 
                           sense2Tag, 0, true );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  EntityHandle sense_data[2] = {0,0};
  rval = mdbImpl->tag_get_data( sense2Tag, &surface, 1, sense_data );
  if (MB_TAG_NOT_FOUND != rval && MB_SUCCESS != rval)
    return MB_FAILURE;
  
  if (sense_data[!forward] == volume)
    return MB_SUCCESS;
  else if (sense_data[!forward])
    return MB_MULTIPLE_ENTITIES_FOUND;
  
  sense_data[!forward] = volume;
  return mdbImpl->tag_set_data( sense2Tag, &surface, 1, sense_data );
}

ErrorCode GeomTopoUtil::get_sense( EntityHandle surface,
                                     EntityHandle volume,
                                     bool& forward )
{
  ErrorCode rval;
  if (!sense2Tag) {
    rval = mdbImpl->tag_get_handle( GEOM_SENSE_TAG_NAME, sense2Tag );
    if (MB_SUCCESS != rval) {
      sense2Tag = 0;
      return MB_FAILURE;
    }
  }
  
  EntityHandle sense_data[2] = {0,0};
  rval = mdbImpl->tag_get_data( sense2Tag, &surface, 1, sense_data );
  if (MB_SUCCESS != rval)
    return rval;
  
  if (sense_data[0] == volume)
    forward = true;
  else if (sense_data[1] == volume)
    forward = false;
  else
    return MB_ENTITY_NOT_FOUND;
  
  return MB_SUCCESS;
}

ErrorCode GeomTopoUtil::find_geomsets(Range *ranges) 
{
  Tag geom_tag;
  ErrorCode result = mdbImpl->tag_create(GEOM_DIMENSION_TAG_NAME, 4, 
                                           MB_TAG_SPARSE, geom_tag, NULL);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
    return result;
  
    // get all sets with this tag
  Range geom_sets;
  result = mdbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &geom_tag, NULL, 1,
                                                 geom_sets);
  if (MB_SUCCESS != result || geom_sets.empty()) 
    return result;

  result = separate_by_dimension(geom_sets, geomRanges, geom_tag);
  if (MB_SUCCESS != result)
    return result;

  if (ranges) {
    for (int i = 0; i < 4; i++) ranges[i] = geomRanges[i];
  }
  
  return MB_SUCCESS;
}
  
    //! Restore parent/child links between GEOM_TOPO mesh sets
ErrorCode GeomTopoUtil::restore_topology() 
{
  
    // look for geometric topology sets and restore parent/child links between them
    // algorithm:
    // - for each entity of dimension d=D-1..0:
    //   . get d-dimensional entity in entity
    //   . get all (d+1)-dim adjs to that entity
    //   . for each geom entity if dim d+1, if it contains any of the ents,
    //     add it to list of parents
    //   . make parent/child links with parents

  Range entities[4];
  result = find_geomsets(entities);
  if (MB_SUCCESS != result)
    return result;

  std::vector<EntityHandle> parents;
  Range tmp_parents;
  
    // loop over dimensions
  for (int dim = 2; dim >= 0; dim--) {
      // mark entities of next higher dimension with their owners; regenerate tag
      // each dimension so prev dim's tag data goes away
    Tag owner_tag;
    EntityHandle dum_val = 0;
    result = mdbImpl->tag_create("__owner_tag", sizeof(EntityHandle), MB_TAG_DENSE,
                                 MB_TYPE_HANDLE, owner_tag, &dum_val);
    if (MB_SUCCESS != result) continue;
    Range dp1ents;
    std::vector<EntityHandle> owners;
    for (Range::iterator rit = entities[dim+1].begin(); rit != entities[dim+1].end(); rit++) {
      dp1ents.clear();
      result = mdbImpl->get_entities_by_dimension(*rit, dim+1, dp1ents);
      if (MB_SUCCESS != result) continue;
      owners.resize(dp1ents.size());
      std::fill(owners.begin(), owners.end(), *rit);
      result = mdbImpl->tag_set_data(owner_tag, dp1ents, &owners[0]);
      if (MB_SUCCESS != result) continue;
    }
    
    for (Range::iterator d_it = entities[dim].begin(); 
         d_it != entities[dim].end(); d_it++) {
      Range dents;
      result = mdbImpl->get_entities_by_dimension(*d_it, dim, dents);
      if (MB_SUCCESS != result) continue;
      if (dents.empty()) continue;
      
        // get (d+1)-dimensional adjs
      dp1ents.clear();
      result = mdbImpl->get_adjacencies(&(*dents.begin()), 1, dim+1, 
                                        false, dp1ents);
      if (MB_SUCCESS != result || dp1ents.empty()) continue;

        // get owner tags
      parents.resize(dp1ents.size());
      result = mdbImpl->tag_get_data(owner_tag, dp1ents, &parents[0]);
      assert(MB_TAG_NOT_FOUND != result);
      if (MB_SUCCESS != result) continue;
      
        // compress to a range to remove duplicates
      tmp_parents.clear();
      std::copy(parents.begin(), parents.end(), range_inserter(tmp_parents));
      for (Range::iterator pit = tmp_parents.begin(); pit != tmp_parents.end(); pit++) {
        result = mdbImpl->add_parent_child(*pit, *d_it);
        if (MB_SUCCESS != result) return result;
      }
      
        // store surface senses
      if (dim != 2) 
        continue;
      const EntityHandle *conn3, *conn2;
      int len3, len2, err, num, sense, offset;
      for (size_t i = 0; i < parents.size(); ++i) {
        result = mdbImpl->get_connectivity( dp1ents[i], conn3, len3, true );
        if (MB_SUCCESS != result) return result;
        result = mdbImpl->get_connectivity( dents.front(), conn2, len2, true );
        if (MB_SUCCESS != result) return result;
        assert(len2 <= 4);
        err = MBCN::SideNumber( TYPE_FROM_HANDLE(dp1ents[i]), conn3,
                                conn2, len2, dim, num, sense, offset );
        if (err)
          return MB_FAILURE;
        
        result = set_sense( *d_it, parents[i], sense == 1 );
        if (MB_MULTIPLE_ENTITIES_FOUND == result) {
          std::cerr << "Warning: Multiple volumes use surface with same sense." << std::endl
                    << "         Some geometric sense data lost." << std::endl;
        }
        else if (MB_SUCCESS != result) {
          return result;
        }
      }
    }
    
      // now delete owner tag on this dimension, automatically removes tag data
    result = mdbImpl->tag_delete(owner_tag);
    if (MB_SUCCESS != result) return result;
    
  } // dim

  return result;
}

ErrorCode GeomTopoUtil::separate_by_dimension(const Range &geom_sets,
                                                 Range *entities, Tag geom_tag) 
{
  ErrorCode result;
  
  if (0 == geom_tag) {
    
    result = mdbImpl->tag_get_handle(GEOM_DIMENSION_TAG_NAME, geom_tag);
    if (MB_SUCCESS != result)
      return result;
  }

    // get the data for those tags
  std::vector<int> tag_vals(geom_sets.size());
  result = mdbImpl->tag_get_data(geom_tag, geom_sets, &tag_vals[0]);
  if (MB_SUCCESS != result)
    return result;

  Range::const_iterator git;
  std::vector<int>::iterator iit;
  
  for (git = geom_sets.begin(), iit = tag_vals.begin(); git != geom_sets.end(); 
       git++, iit++) {
    if (0 <= *iit && 3 >= *iit) 
      entities[*iit].insert(*git);
    else {
        // assert(false);
        // do nothing for now
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode GeomTopoUtil::construct_vertex_ranges(const Range &geom_sets,
                                                   const Tag verts_tag) 
{
    // construct the vertex range for each entity and put on that tag
  Range *temp_verts, temp_elems;
  ErrorCode result = MB_SUCCESS;
  for (Range::const_iterator it = geom_sets.begin(); it != geom_sets.end(); it++) {
      // make the new range
    temp_verts = new Range();
    assert(NULL != temp_verts);
    temp_elems.clear();
    
      // get all the elements in the set, recursively
    result = mdbImpl->get_entities_by_handle(*it, temp_elems, true);
    if (MB_SUCCESS != result) 
      return result;
    
      // get all the verts of those elements; use get_adjacencies 'cuz it handles ranges better
    result = mdbImpl->get_adjacencies(temp_elems, 0, false, *temp_verts,
                                      Interface::UNION);
    if (MB_SUCCESS != result) 
      return result;

      // store this range as a tag on the entity
    result = mdbImpl->tag_set_data(verts_tag, &(*it), 1, &temp_verts);
    if (MB_SUCCESS != result) 
      return result;
    
  }
  
  return result;
}

  
  
