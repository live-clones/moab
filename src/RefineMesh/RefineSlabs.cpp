
#include "moab/RefineSlabs.hpp"
#include "moab/HalfFacetRep.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/Util.hpp"
#include <iostream>

// turn on debugging asserts
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>

#include <vector>
#include <utility>
#include <limits>
#include <cmath>

namespace moab{

  /**********************
    *  AHF related routines   *
    *********************/
  ErrorCode RefineSlabs::get_adjacencies(const EntityHandle source_entity,
                                          const unsigned int target_dimension,
                                          std::vector<EntityHandle> &target_entities)
  {
    ErrorCode error;
    error = ahf->get_adjacencies(source_entity, target_dimension, target_entities); MB_CHK_ERR(error);

    return MB_SUCCESS;
  }

  ErrorCode RefineSlabs::estimate_allocate_ahf_maps( size_t num_hexes_memory_estimate )
  {
    ErrorCode error;
    error = ahf->resize_hf_maps(3, num_hexes_memory_estimate, num_hexes_memory_estimate);MB_CHK_ERR(error);

    return MB_SUCCESS;
  }

  ErrorCode RefineSlabs::replace_mesh( Entities /*&coarse_hexes*/, Entities &coarse_bndverts, Entities &coarse_opphexes, Entities /*&fine_hexes*/)
  {
    ErrorCode error;
    Entities all_fhexes;

    //First, get the corresponding fine duplicate vertices and then collect all the fine hexes incident on it.
    for (size_t i = 0; i < coarse_bndverts.size(); i++)
      {
        EntityHandle cnode, fnode;
        cnode = coarse_bndverts[i];
        bool is_fine = get_copy(cnode, fnode);

        if (!is_fine)
          MB_SET_ERR(MB_FAILURE, "No fine copy of the coarse boundary vertice");

        Entities fhexes;
       // get_all_hexes(fnode, fhexes, true);
         get_all_hexes(fnode, fhexes, false);
        all_fhexes.insert(all_fhexes.end(), fhexes.begin(), fhexes.end());

        for (size_t j = 0; j < fhexes.size(); j++)
          {
            EntityHandle hex_nodes[8];
            get_hex_nodes(fhexes[j], hex_nodes);
            for (size_t k = 0; k < 8; k++)
              {
                if (hex_nodes[k] == fnode)
                  {
                    hex_nodes[k] = cnode;
                    break;
                  }
              }
            error = mbImpl->set_connectivity(fhexes[j], &hex_nodes[0], 8);MB_CHK_ERR(error);
          }
      }

    std::sort(all_fhexes.begin(), all_fhexes.end());
    all_fhexes.erase(std::unique(all_fhexes.begin(), all_fhexes.end()), all_fhexes.end());

    //step 5: join across the shared coarse-fine boundary
    error = ahf->update_hf_maps_multimesh(3, coarse_opphexes, all_fhexes);MB_CHK_ERR(error);

    return MB_SUCCESS;
  }
  ErrorCode RefineSlabs::replace_node( EntityHandle chex, int node_lid, EntityHandle new_node )
  {
    EntityHandle nodes[8];
    ErrorCode error = get_hex_nodes(chex, nodes);MB_CHK_ERR(error);

    nodes[node_lid] = new_node;
    error = mbImpl->set_connectivity(chex, &nodes[0], 8);MB_CHK_ERR(error);

    return MB_SUCCESS;
  }

  ErrorCode RefineSlabs::update_ahf_maps(Entities /*&shrink_set*/, Entities &new_hexes)
  {
    ErrorCode error = ahf->update_hf_maps_patch(3, new_hexes);MB_CHK_ERR(error);
     return MB_SUCCESS;
  }

  ErrorCode RefineSlabs::update_ahf_maps(Entities &coarse_hexes)
  {
    ErrorCode error = ahf->update_hf_maps_patch(3, coarse_hexes);MB_CHK_ERR(error);
    error = ahf->print_tags(3);MB_CHK_ERR(error);
    return MB_SUCCESS;
  }
  
  ErrorCode RefineSlabs::get_all_hexes( EntityHandle node, Entities &hexes, bool iscoarse )
  {
    // debug
    node_ok(iscoarse, node);
  
    hexes.clear();
    ErrorCode error = ahf->get_up_adjacencies_vert_3d(node, hexes);MB_CHK_ERR(error);

    return MB_SUCCESS;
  }

  ErrorCode RefineSlabs::get_all_quads(EntityHandle node, Entities &quads, bool iscoarse )
  {
    // debug
    node_ok(iscoarse, node);

    quads.clear();
    ErrorCode error = ahf->get_up_adjacencies_vert_2d( node, quads );MB_CHK_ERR(error);

    return MB_SUCCESS;
  }

  ErrorCode RefineSlabs::get_coarse_opp_hexes(Entities &coarse_hexes, Entities &coarse_opphexes)
  {
    ErrorCode error;
    int index = ahf->get_index_in_lmap(coarse_hexes[0]);
    int nfpc = ahf->lConnMap3D[index].num_faces_in_cell;

    for (size_t i = 0; i < coarse_hexes.size(); ++i)
      {
        EntityHandle hex_nodes[8];
        error = get_hex_nodes(coarse_hexes[i], hex_nodes);MB_CHK_ERR(error);

        for (int j = 0; j <nfpc; j++)
          {
            int count = 0;
            for (size_t k = 0; k < 4; k++)
              {
                int id = ahf->lConnMap3D[index].hf2v[j][k];
                if (get_membership(hex_nodes[id]) == SlabData::BOUNDARY)
                  {
                    count += 1;
                  }
              }
            if (count == 4)
              {
               EntityHandle opphex = 0;
                error = ahf->get_neighbor_adjacencies_3d(coarse_hexes[i], j, &opphex);MB_CHK_ERR(error);
                coarse_opphexes.push_back(opphex);
              }
          }
      }

    return MB_SUCCESS;
  }


    // given a SlabEdge (edge and vertex) of a hex, find the other two edges sharing that vertex
  void RefineSlabs::get_adj( const SlabEdge &edge, SlabEdge &adj1, SlabEdge &adj2 )
  {
    // todo AHF, move to AHF
    // this should go into AHF and not be hard-coded here

    adj1.hex = adj2.hex = edge.hex;
    adj1.head_node = adj2.head_node = edge.head_node;    

    // debug. we just call this with coarse entities
    assert(is_registered_hex(edge.hex));
    assert(is_registered_vertex(edge.head_node));
    assert(is_registered_vertex(edge.tail_node));
    
    const int head_index = get_hex_node_index( edge.hex, edge.head_node );
    const int tail_index = get_hex_node_index( edge.hex, edge.tail_node );
    const int edge_lid = edge.edge_lid;
    assert( head_index >= 0 && head_index < 8);
    assert( tail_index >= 0 && tail_index < 8);
    assert( head_index != tail_index );
    assert( edge_lid >= 0 && edge_lid < 12);


    EntityHandle hex_nodes[8];
    get_hex_nodes( edge.hex, hex_nodes);

    // todo AHF, move to AHF
    // 24 different edges and orientations, check for each one
    switch ( head_index )
    {
      case 0:
      switch( tail_index)
      {        
        case 1:
        assert( edge_lid == 0 );
        adj1.tail_node = hex_nodes[3];  adj1.edge_lid = 3;
        adj2.tail_node = hex_nodes[4];  adj2.edge_lid = 4;
        break;
        case 3:
        assert( edge_lid == 3 );
        adj1.tail_node = hex_nodes[1];  adj1.edge_lid = 0;
        adj2.tail_node = hex_nodes[4];  adj2.edge_lid = 4;
        break;
        case 4:
        assert( edge_lid == 4 );
        adj1.tail_node = hex_nodes[3];  adj1.edge_lid = 3;
        adj2.tail_node = hex_nodes[1];  adj2.edge_lid = 0;
        break;
        default:
          assert(0);
      }
      break;

      case 1:
      switch( tail_index)
      {        
        case 0:
        assert( edge_lid == 0 );
        adj1.tail_node = hex_nodes[5];  adj1.edge_lid = 5;
        adj2.tail_node = hex_nodes[2];  adj2.edge_lid = 1;
        break;
        case 5:
        assert( edge_lid == 5 );
        adj1.tail_node = hex_nodes[0];  adj1.edge_lid = 0;
        adj2.tail_node = hex_nodes[2];  adj2.edge_lid = 1;
        break;
        case 2:
        assert( edge_lid == 1 );
        adj1.tail_node = hex_nodes[5];  adj1.edge_lid = 5;
        adj2.tail_node = hex_nodes[0];  adj2.edge_lid = 0;
        break;
        default:
          assert(0);
      }
      break;

      case 2:
      switch( tail_index)
      {        
        case 1:
        assert( edge_lid == 1 );
        adj1.tail_node = hex_nodes[6];  adj1.edge_lid = 6;
        adj2.tail_node = hex_nodes[3];  adj2.edge_lid = 2;
        break;
        case 3:
        assert( edge_lid == 2 );
        adj1.tail_node = hex_nodes[1];  adj1.edge_lid = 1;
        adj2.tail_node = hex_nodes[6];  adj2.edge_lid = 6;
        break;
        case 6:
        assert( edge_lid == 6 );
        adj1.tail_node = hex_nodes[1];  adj1.edge_lid = 1;
        adj2.tail_node = hex_nodes[3];  adj2.edge_lid = 2;
        break;
        default:
          assert(0);
      }
      break;


      case 3:
      switch( tail_index)
      {        
        case 2:
        assert( edge_lid == 2 );
        adj1.tail_node = hex_nodes[0];  adj1.edge_lid = 3;
        adj2.tail_node = hex_nodes[7];  adj2.edge_lid = 7;
        break;
        case 0:
        assert( edge_lid == 3 );
        adj1.tail_node = hex_nodes[2];  adj1.edge_lid = 2;
        adj2.tail_node = hex_nodes[7];  adj2.edge_lid = 7;
        break;
        case 7:
        assert( edge_lid == 7 );
        adj1.tail_node = hex_nodes[2];  adj1.edge_lid = 2;
        adj2.tail_node = hex_nodes[0];  adj2.edge_lid = 3;
        break;
        default:
          assert(0);
      }
      break;


      case 4:
      switch( tail_index)
      {        
        case 0:
        assert( edge_lid == 4 );
        adj1.tail_node = hex_nodes[5];  adj1.edge_lid = 8;
        adj2.tail_node = hex_nodes[7];  adj2.edge_lid = 11;
        break;
        case 5:
        assert( edge_lid == 8 );
        adj1.tail_node = hex_nodes[0];  adj1.edge_lid = 4;
        adj2.tail_node = hex_nodes[7];  adj2.edge_lid = 11;
        break;
        case 7:
        assert( edge_lid == 11);
        adj1.tail_node = hex_nodes[0];  adj1.edge_lid = 4;
        adj2.tail_node = hex_nodes[5];  adj2.edge_lid = 8;
        break;
        default:
          assert(0);
      }
      break;

      case 5:
      switch( tail_index)
      {        
        case 1:
        assert( edge_lid == 5 );
        adj1.tail_node = hex_nodes[4];  adj1.edge_lid = 8;
        adj2.tail_node = hex_nodes[6];  adj2.edge_lid = 9;
        break;
        case 4:
        assert( edge_lid == 8 );
        adj1.tail_node = hex_nodes[1];  adj1.edge_lid = 5;
        adj2.tail_node = hex_nodes[6];  adj2.edge_lid = 9;
        break;
        case 6:
        assert( edge_lid == 9 );
        adj1.tail_node = hex_nodes[1];  adj1.edge_lid = 5;
        adj2.tail_node = hex_nodes[4];  adj2.edge_lid = 8;
        break;
        default:
          assert(0);
      }
      break;



      case 6:
      switch( tail_index)
      {        
        case 2:
        assert( edge_lid == 6 );
        adj1.tail_node = hex_nodes[5];  adj1.edge_lid = 9;
        adj2.tail_node = hex_nodes[7];  adj2.edge_lid = 10;
        break;
        case 5:
        assert( edge_lid == 9 );
        adj1.tail_node = hex_nodes[2];  adj1.edge_lid = 6;
        adj2.tail_node = hex_nodes[7];  adj2.edge_lid = 10;
        break;
        case 7:
        assert( edge_lid == 10);
        adj1.tail_node = hex_nodes[2];  adj1.edge_lid = 6;
        adj2.tail_node = hex_nodes[5];  adj2.edge_lid = 9;
        break;
        default:
          assert(0);
      }
      break;



      case 7:
      switch( tail_index)
      {        
        case 3:
        assert( edge_lid == 7 );
        adj1.tail_node = hex_nodes[4];  adj1.edge_lid = 11;
        adj2.tail_node = hex_nodes[6];  adj2.edge_lid = 10;
        break;
        case 4:
        assert( edge_lid == 11);
        adj1.tail_node = hex_nodes[3];  adj1.edge_lid = 7;
        adj2.tail_node = hex_nodes[6];  adj2.edge_lid = 10;
        break;
        case 6:
        assert( edge_lid == 10);
        adj1.tail_node = hex_nodes[3];  adj1.edge_lid = 7;
        adj2.tail_node = hex_nodes[4];  adj2.edge_lid = 11;
        break;
        default:
          assert(0);
      }
      break;
      default:
        assert(0);
    }
    assert( adj1.edge_lid != adj2.edge_lid && adj1.edge_lid != edge_lid && adj2.edge_lid != edge_lid );
    assert( adj1.head_node == adj2.head_node && adj1.head_node == edge.head_node );
    assert( adj1.tail_node != edge.tail_node && adj1.tail_node != adj2.tail_node && adj2.tail_node != edge.tail_node );
    assert( adj1.head_node != adj1.tail_node );
    assert( adj2.head_node != adj2.tail_node );
  }

  // given a SlabEdge (edge and vertex) of a hex, find the other two edges that are opposite the edge in a quad of the hex
  void RefineSlabs::get_opp( const SlabEdge &edge, SlabEdge &opp1, SlabEdge &opp2 ) 
  {
    // todo AHF, move to AHF
    // this should go into AHF and not be hard-coded here
    opp1.hex = edge.hex;
    opp2.hex = edge.hex;

    const int head_index = get_hex_node_index( edge.hex, edge.head_node );
    const int tail_index = get_hex_node_index( edge.hex, edge.tail_node );
    const int edge_lid = edge.edge_lid;
    assert( head_index >= 0 && head_index < 8);
    assert( tail_index >= 0 && tail_index < 8);
    assert( edge_lid >= 0 && edge_lid < 12);


    EntityHandle hex_nodes[8];
    get_hex_nodes( edge.hex, hex_nodes);

    // todo AHF, move to AHF
    // 24 different edges and orientations, check for each one
    switch ( head_index )
    {
      case 0:
      switch( tail_index)
      {        
        case 1:
        assert( edge_lid == 0 );
        opp1.head_node = hex_nodes[3];  opp1.tail_node = hex_nodes[2];  opp1.edge_lid = 2;
        opp2.head_node = hex_nodes[4];  opp2.tail_node = hex_nodes[5];  opp2.edge_lid = 8;
        break;
        case 3:
        assert( edge_lid == 3 );
        opp1.head_node = hex_nodes[1];  opp1.tail_node = hex_nodes[2];  opp1.edge_lid = 1;
        opp2.head_node = hex_nodes[4];  opp2.tail_node = hex_nodes[7];  opp2.edge_lid = 11;
        break;
        case 4:
        assert( edge_lid == 4 );
        opp1.head_node = hex_nodes[1];  opp1.tail_node = hex_nodes[5];  opp1.edge_lid = 5;
        opp2.head_node = hex_nodes[3];  opp2.tail_node = hex_nodes[7];  opp2.edge_lid = 7;
        break;
        default:
          assert(0);
      }
      break;

      case 1:
      switch( tail_index)
      {        
        case 0:
        assert( edge_lid == 0);
        opp1.head_node = hex_nodes[2];  opp1.tail_node = hex_nodes[3];  opp1.edge_lid = 2;
        opp2.head_node = hex_nodes[5];  opp2.tail_node = hex_nodes[4];  opp2.edge_lid = 8;
        break;
        case 5:
        assert( edge_lid == 5 );
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[4];  opp1.edge_lid = 4;
        opp2.head_node = hex_nodes[2];  opp2.tail_node = hex_nodes[6];  opp2.edge_lid = 6;
        break;
        case 2:
        assert( edge_lid == 1 );
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[3];  opp1.edge_lid = 3;
        opp2.head_node = hex_nodes[5];  opp2.tail_node = hex_nodes[6];  opp2.edge_lid = 9;
        break;
        default:
          assert(0);
      }
      break;

      case 2:
      switch( tail_index)
      {        
        case 1:
        assert( edge_lid == 1 );
        opp1.head_node = hex_nodes[3];  opp1.tail_node = hex_nodes[0];  opp1.edge_lid = 3;
        opp2.head_node = hex_nodes[6];  opp2.tail_node = hex_nodes[5];  opp2.edge_lid = 9;
        break;
        case 3:
        assert( edge_lid == 2 );
        opp1.head_node = hex_nodes[1];  opp1.tail_node = hex_nodes[0];  opp1.edge_lid = 0;
        opp2.head_node = hex_nodes[6];  opp2.tail_node = hex_nodes[7];  opp2.edge_lid = 10;
        break;
        case 6:
        assert( edge_lid == 6 );
        opp1.head_node = hex_nodes[1];  opp1.tail_node = hex_nodes[5];  opp1.edge_lid = 5;
        opp2.head_node = hex_nodes[3];  opp2.tail_node = hex_nodes[7];  opp2.edge_lid = 7;
        break;
        default:
          assert(0);
      }
      break;


      case 3:
      switch( tail_index)
      {        
        case 2:
        assert( edge_lid == 2 );
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[1];  opp1.edge_lid = 0;
        opp2.head_node = hex_nodes[7];  opp2.tail_node = hex_nodes[6];  opp2.edge_lid = 10;
        break;
        case 0:
        assert( edge_lid == 3 );
        opp1.head_node = hex_nodes[2];  opp1.tail_node = hex_nodes[1];  opp1.edge_lid = 1;
        opp2.head_node = hex_nodes[7];  opp2.tail_node = hex_nodes[4];  opp2.edge_lid = 11;
        break;
        case 7:
        assert( edge_lid == 7 );
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[4];  opp1.edge_lid = 4;
        opp2.head_node = hex_nodes[2];  opp2.tail_node = hex_nodes[6];  opp2.edge_lid = 6;
        break;
        default:
          assert(0);
      }
      break;


      case 4:
      switch( tail_index)
      {        
        case 0:
        assert( edge_lid == 4 );
        opp1.head_node = hex_nodes[5];  opp1.tail_node = hex_nodes[1];  opp1.edge_lid = 5;
        opp2.head_node = hex_nodes[7];  opp2.tail_node = hex_nodes[3];  opp2.edge_lid = 7;
        break;
        case 5:
        assert( edge_lid == 8 );
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[1];  opp1.edge_lid = 0;
        opp2.head_node = hex_nodes[7];  opp2.tail_node = hex_nodes[6];  opp2.edge_lid = 10;
        break;
        case 7:
        assert( edge_lid == 11);
        opp1.head_node = hex_nodes[0];  opp1.tail_node = hex_nodes[3];  opp1.edge_lid = 3;
        opp2.head_node = hex_nodes[5];  opp2.tail_node = hex_nodes[6];  opp2.edge_lid = 9;
        break;
        default:
          assert(0);
      }
      break;

      case 5:
      switch( tail_index)
      {        
        case 1:
        assert(edge_lid == 5);
        opp1.head_node = hex_nodes[4];  opp1.tail_node = hex_nodes[0];  opp1.edge_lid = 4;
        opp2.head_node = hex_nodes[6];  opp2.tail_node = hex_nodes[2];  opp2.edge_lid = 6;
        break;
        case 4:
        assert( edge_lid == 8 );
        opp1.head_node = hex_nodes[1];  opp1.tail_node = hex_nodes[0];  opp1.edge_lid = 0;
        opp2.head_node = hex_nodes[6];  opp2.tail_node = hex_nodes[7];  opp2.edge_lid = 10;
        break;
        case 6:
        assert( edge_lid == 9 );
        opp1.head_node = hex_nodes[1];  opp1.tail_node = hex_nodes[2];  opp1.edge_lid = 1;
        opp2.head_node = hex_nodes[4];  opp2.tail_node = hex_nodes[7];  opp2.edge_lid = 11;
        break;
        default:
          assert(0);
      }
      break;



      case 6:
      switch( tail_index)
      {        
        case 2:
        assert( edge_lid == 6 );
        opp1.head_node = hex_nodes[5];  opp1.tail_node = hex_nodes[1];  opp1.edge_lid = 5;
        opp2.head_node = hex_nodes[7];  opp2.tail_node = hex_nodes[3];  opp2.edge_lid = 7; 
        break;
        case 5:
        assert( edge_lid == 9 );
        opp1.head_node = hex_nodes[2];  opp1.tail_node = hex_nodes[1];  opp1.edge_lid = 1;
        opp2.head_node = hex_nodes[7];  opp2.tail_node = hex_nodes[4];  opp2.edge_lid = 11;
        break;
        case 7:
        assert( edge_lid == 10);
        opp1.head_node = hex_nodes[2];  opp1.tail_node = hex_nodes[3];  opp1.edge_lid = 2;
        opp2.head_node = hex_nodes[5];  opp2.tail_node = hex_nodes[4];  opp2.edge_lid = 8;
        break;
        default:
          assert(0);
      }
      break;



      case 7:
      switch( tail_index)
      {        
        case 3:
        assert( edge_lid == 7 );
        opp1.head_node = hex_nodes[4];  opp1.tail_node = hex_nodes[0];  opp1.edge_lid = 4;
        opp2.head_node = hex_nodes[6];  opp2.tail_node = hex_nodes[2];  opp2.edge_lid = 6;
        break;
        case 4:
        assert( edge_lid == 11);
        opp1.head_node = hex_nodes[3];  opp1.tail_node = hex_nodes[0];  opp1.edge_lid = 3;
        opp2.head_node = hex_nodes[6];  opp2.tail_node = hex_nodes[5];  opp2.edge_lid = 9;
        break;
        case 6:
        assert( edge_lid == 10);
        opp1.head_node = hex_nodes[3];  opp1.tail_node = hex_nodes[2];  opp1.edge_lid = 2;
        opp2.head_node = hex_nodes[4];  opp2.tail_node = hex_nodes[5];  opp2.edge_lid = 8;
        break;
        default:
          assert(0);
      }
      break;

      default:
        assert(0);
    }

#ifdef NDEBUG
    {
    // verification
    // head of opposites must be edge-adjacent to the head, and not be the tail
    // tail of opposites must be edge-adjacent to the tail, and not be the head
      SlabEdge hadj1, hadj2;
      get_adj( edge, hadj1, hadj2 );
      assert( opp1.head_node == hadj1.tail_node || opp1.head_node == hadj2.tail_node );
      assert( opp2.head_node == hadj1.tail_node || opp2.head_node == hadj2.tail_node );
      assert( opp1.head_node != opp2.head_node );

      SlabEdge reverse_edge( edge );
      reverse_edge.flip();
      SlabEdge tadj1, tadj2;
      get_adj( edge, tadj1, tadj2 );
      assert( opp1.tail_node == tadj1.tail_node || opp1.tail_node == tadj2.tail_node );
      assert( opp2.tail_node == hadj1.tail_node || opp2.tail_node == tadj2.tail_node );
      assert( opp1.tail_node != opp2.tail_node );
    }
#endif

  }

  ErrorCode RefineSlabs::create_hex( EntityHandle fine_nodes[8], EntityHandle & new_hex )
  {
    ErrorCode error = mbImpl->create_element(MBHEX, fine_nodes, 8, new_hex);
    assert( MB_SUCCESS == error );
    // refinement_AHF
    // tell AHF about the new hex // AHF todo
    // alternatively, as in the pseudocode, AHF can do this once after all the hexes have been created
    //NR: Don't change anything in AHF here.

    //debug
    created_fine_hexes.insert(new_hex);

    return error;
  }

  ErrorCode RefineSlabs::create_node( EntityHandle node, EntityHandle &new_node )
  {
    // get coordinates from node
    const double *xp, *yp, *zp;
    ErrorCode error = mbImpl->get_coords( node, xp, yp, zp );MB_CHK_ERR(error);
   // assert(  MB_SUCCESS == error );


  //  double coord[3];
 //   coord[0] = *xp;
//    coord[1] = *yp;
 //   coord[2] = *zp;
 //   error = mbImpl->create_vertex(coord, new_node);MB_CHK_ERR(error);
   // assert( MB_SUCCESS == error );

    //new strategy for point positioning
    const double eps = 1.0/3.0;
    Entities adjents;
    error = ahf->get_up_adjacencies_vert_3d(node, adjents);
    if (adjents.size() == 1)
      {
        EntityHandle hnodes[8];
        get_hex_nodes(adjents[0], hnodes);
        double centroid[3] = {0,0,0};
        const double *xh, *yh, *zh;
        for (size_t n=0; n<8; ++n)
          {
            error = mbImpl->get_coords(hnodes[n], xh, yh, zh);MB_CHK_ERR(error);
            centroid[0] += *xh; centroid[1] += *yh; centroid[2] += *zh;
          }
        centroid[0] = centroid[0]/8; centroid[1] = centroid[1]/8; centroid[2] = centroid[2]/8;
        double dist = sqrt(((*xp-centroid[0])^2)+((*yp-centroid[1])^2)+((*zp-centroid[2])^2));
        double coords[3] = {0,0,0};
        coords[0] = *xp + (eps/dist) * (centroid[0]-*xp);
        coords[1] = *yp + (eps/dist) * (centroid[1]-*yp);
        coords[2] = *zp + (eps/dist) * (centroid[2]-*zp);
        error = mbImpl->create_vertex(coords, new_node);MB_CHK_ERR(error);
      }
    else
      {
        Entities tnodes[8], hnodes[8], cnodes;
        std::vector<EntityHandle>::iterator it ;
        for (int j=0; j<(int)adjents.size(); j++)
          {
            get_hex_nodes(adjents[j], &hnodes[0]);
            std::sort(hnodes.begin(), hnodes.end());
            if (j==0)
              tnodes.insert(tnodes.end(), hnodes.begin(), hnodes.end());
            else
              {
                it = std::set_intersection(hnodes.begin(), hnodes.end(), tnodes.begin(), tnodes.end(), cnodes.begin());
                cnodes.resize(it-cnodes.begin());
              }
          }
        double centroid[3] = {0,0,0};
        const double *xh, *yh, *zh;
        int numc = cnodes.size();
        if (numc == 1)
          {
            MB_SET_ERR(MB_FAILURE,"Check configuration of the shrink set ");
          }
        else
          {
            for (size_t n=0; n<numc; ++n)
              {
                error = mbImpl->get_coords(cnodes[n], xh, yh, zh);MB_CHK_ERR(error);
                centroid[0] += *xh; centroid[1] += *yh; centroid[2] += *zh;
              }
            centroid[0] = centroid[0]/numc; centroid[1] = centroid[1]/num; centroid[2] = centroid[2]/num;
            double dist = sqrt(((*xp-centroid[0])^2)+((*yp-centroid[1])^2)+((*zp-centroid[2])^2));
            double coords[3] = {0,0,0};
            coords[0] = *xp + (eps/dist) * (centroid[0]-*xp);
            coords[1] = *yp + (eps/dist) * (centroid[1]-*yp);
            coords[2] = *zp + (eps/dist) * (centroid[2]-*zp);
            error = mbImpl->create_vertex(coords, new_node);MB_CHK_ERR(error);
          }
      }

    // refinement_AHF
    // tell AHF about the new node // AHF todo
    // alternatively, as in the pseudocode, AHF can do this once after all the hexes have been created

    //debug
    created_fine_nodes.insert(new_node);

    return error;
  }

  // get the dimension of the geometric object that this mesh entity lies on
  // E.g. 3 if inside the volume, 2 if on its surface, 1 if on an edge of the surface, ...
  int RefineSlabs::get_geometry_dimension( EntityHandle entity_handle )
  {
    // geometry_map[ entity_handle ] == dimension
    // auto it =
    std::map< EntityHandle, int >::iterator
    it = geometry_map.find( entity_handle );
    if ( it == geometry_map.end() )
      return 3;
    return it->second;
  }


  // ========================================================
  // inlined functions placed here for debugging
  #ifndef NDEBUG
  #define MOAB_RefineSlabs_inline inline
  #else
  #define MOAB_RefineSlabs_inline
  #endif

  MOAB_RefineSlabs_inline
  int RefineSlabs::get_hex_edge_index( int a_index, int b_index )
  {
    // call function from AHF, see HalfFacetRef.hpp
    return HalfFacetRep::lConnMap3D[MBHEX].lookup_leids[a_index][b_index];
  }

  //==== SlabEntity methods
  MOAB_RefineSlabs_inline
  bool RefineSlabs::SlabHex::operator==(const SlabHex& rhs) const
  {
    return (entity_handle == rhs.entity_handle);
  }

  MOAB_RefineSlabs_inline
  bool RefineSlabs::SlabNode::operator==(const SlabNode& rhs) const
  {
    return (entity_handle == rhs.entity_handle);
  }


  MOAB_RefineSlabs_inline
  bool RefineSlabs::SlabEdge::operator==(const SlabEdge& rhs) const
  {
    return (head_node == rhs.head_node) &&
    (tail_node == rhs.tail_node);
              // hex and local id need not match
               // && (hex == rhs.hex);
  }

  MOAB_RefineSlabs_inline
  bool RefineSlabs::SlabEdge::nodes_match( const SlabEdge &rhs ) const
  {
    return ( head_node == rhs.head_node && tail_node == rhs.tail_node) 
    || ( head_node == rhs.tail_node && tail_node == rhs.head_node);
  }

  MOAB_RefineSlabs_inline
  bool  RefineSlabs::SlabEdge::directions_match( const SlabEdge &rhs ) const
  {
    return ( head_node == rhs.head_node && tail_node == rhs.tail_node);
  }

  MOAB_RefineSlabs_inline
  void RefineSlabs::SlabEdge::flip()
  {
    EntityHandle tmp = head_node;
    head_node = tail_node;
    tail_node = tmp;
  }

  MOAB_RefineSlabs_inline
  void RefineSlabs::SlabEdge::assignment( const SlabEdge &copy_me )
  { 
    edge_lid  = copy_me.edge_lid;
    head_node = copy_me.head_node;
    tail_node = copy_me.tail_node;
    hex       = copy_me.hex;
  }

  MOAB_RefineSlabs_inline
  void RefineSlabs::SlabData::copy_data( SlabData *copy_me )
  {
    membership = copy_me->membership;
    shrink_membership = copy_me->shrink_membership;
  }



  MOAB_RefineSlabs_inline
  RefineSlabs::HexRefinement *RefineSlabs::get_hex_refinement( EntityHandle coarse_hex )
  {
    SlabData *slab_data = get_slab_data( coarse_hex );
    if ( !slab_data )
      return 0;
    return slab_data->refinement;
  }
  MOAB_RefineSlabs_inline
  RefineSlabs::HexRefinement *RefineSlabs::force_hex_refinement( EntityHandle coarse_hex )
  {
    SlabData *slab_data = force_slab_data( coarse_hex );
    assert( slab_data );
    if ( !slab_data->refinement )
      slab_data->refinement = new HexRefinement;
    return slab_data->refinement;
  }

  MOAB_RefineSlabs_inline
  RefineSlabs::HexCoarsening *RefineSlabs::get_hex_coarsening( EntityHandle fine_entity )
  {
    SlabData *slab_data = get_slab_data( fine_entity );
    if (!slab_data)
      return 0;
    return slab_data->coarsening;
  }

  MOAB_RefineSlabs_inline
  RefineSlabs::HexCoarsening *RefineSlabs::force_hex_coarsening( EntityHandle fine_entity, EntityHandle coarse_hex )
  {
    SlabData *slab_data = force_slab_data( fine_entity );
    assert( slab_data );
    if ( !slab_data->coarsening )
    {
      HexRefinement *hex_refinement = get_hex_refinement( coarse_hex );
      if ( hex_refinement && !hex_refinement->fine_hexes.empty() )
      {
        slab_data->coarsening = get_hex_coarsening( hex_refinement->fine_hexes.front() );
      }
      else
        slab_data->coarsening = new HexCoarsening(coarse_hex);
    }
    assert( slab_data->coarsening );
    return slab_data->coarsening;
  }

  MOAB_RefineSlabs_inline
  void RefineSlabs::add_refined_hex( EntityHandle coarse_hex, EntityHandle fine_hex )
  {
    force_hex_coarsening( fine_hex, coarse_hex );
    force_hex_refinement( coarse_hex )->fine_hexes.push_back( fine_hex );
  }

  MOAB_RefineSlabs_inline
  int RefineSlabs::get_edge_refinement_level(EntityHandle hex, int edge_lid) 
  {
    return force_hex_refinement(hex)->edge_refinement_level[edge_lid]; 
  }
  MOAB_RefineSlabs_inline
  int RefineSlabs::get_edge_refinement_level(const SlabEdge &slab_edge) 
  {
    return get_edge_refinement_level( slab_edge.hex, slab_edge.edge_lid );
  }

  MOAB_RefineSlabs_inline
  int RefineSlabs::increment_edge_refinement_level(SlabEdges &slab_edges)
  {
    int new_val, new_val_check(-999);
    for ( size_t i = 0; i < slab_edges.size(); ++i )
    {
      new_val = ++(force_hex_refinement(slab_edges[i].hex)->edge_refinement_level[slab_edges[i].edge_lid]); 
      assert( new_val_check == -999 || new_val == new_val_check );
      new_val_check = new_val;
    }
    //debug
    if (!slab_edges.empty())
    {
      std::cout << "++edge " << slab_edges[0].head_node << ":" << slab_edges[0].tail_node << " refinement " << new_val << ". ";
      std::cout.flush();
    }
    return new_val;
  }
  MOAB_RefineSlabs_inline
  int RefineSlabs::decrement_edge_refinement_level(SlabEdges &slab_edges)
  {
    int new_val, new_val_check(999);
    for ( size_t i = 0; i < slab_edges.size(); ++i )
    {
      new_val = --(force_hex_refinement(slab_edges[i].hex)->edge_refinement_level[slab_edges[i].edge_lid]); 
      assert( new_val_check == 999 || new_val == new_val_check );
      new_val_check = new_val;
    }
    //debug
    if (!slab_edges.empty())
    {
      std::cout << "--edge " << slab_edges[0].head_node << ":" << slab_edges[0].tail_node << " refinement " << new_val << ". ";
      std::cout.flush();
    }
    return new_val;
  }

  //==== SlabData get/set methods
  MOAB_RefineSlabs_inline
  RefineSlabs::SlabData *RefineSlabs::get_slab_data( EntityHandle entity_handle) 
  { 
    SlabDataIterator it = slab_data_map.find(entity_handle);
    if ( it == slab_data_map.end() )
      return 0;
    return it->second;
  }
  MOAB_RefineSlabs_inline
  RefineSlabs::SlabData *RefineSlabs::force_slab_data( EntityHandle entity_handle ) 
  { 
    SlabData *slab_data = slab_data_map[entity_handle];
    if ( !slab_data )
      slab_data = slab_data_map[entity_handle] = new SlabData;
    return slab_data;
  }

  MOAB_RefineSlabs_inline
  RefineSlabs::SlabData* RefineSlabs::set_coarse_entity( EntityHandle entity )
  {
    SlabData *slab_data = force_slab_data(entity);
    slab_data->is_coarse = true;
    return slab_data;
  }

  MOAB_RefineSlabs_inline
  void RefineSlabs::set_copy( EntityHandle entity, EntityHandle copy )
  {      
    SlabData *slab_data = force_slab_data( entity );
    slab_data->my_copy_good = true;
    slab_data->my_copy = copy;

    SlabData *slab_data_copy = force_slab_data( copy );
    slab_data_copy->my_copy_good = true;
    slab_data_copy->my_copy = entity;
  }
  MOAB_RefineSlabs_inline
  bool RefineSlabs::get_copy( EntityHandle entity, EntityHandle &copy )
  {
    SlabData * slab_data = get_slab_data( entity );
    if (slab_data && slab_data->my_copy_good)
    {
      copy = slab_data->my_copy;
      return true;
    }
    copy = bad_handle;
    return false;
  }


  MOAB_RefineSlabs_inline
  void RefineSlabs::set_shrunk_node( EntityHandle entity, EntityHandle fine )
  {
    SlabData *slab_data = force_slab_data(entity);
    slab_data->mini_me_good = true;
    slab_data->mini_me = fine;

    SlabData *fine_slab = force_slab_data(fine);
    fine_slab->mini_me_good = true;
    fine_slab->mini_me = fine;
  }

  MOAB_RefineSlabs_inline
  bool RefineSlabs::get_shrunk_node( EntityHandle node, EntityHandle &fine_node )
  {
    SlabData *slab_data = get_slab_data(node);
    if (slab_data && slab_data->mini_me_good)
    {
      fine_node = slab_data->mini_me;
      return true;
    }
    fine_node = bad_handle;
    return false;
  }


  MOAB_RefineSlabs_inline
  bool RefineSlabs::is_coarse( EntityHandle entity_handle )
  {
    SlabData *slab_data = get_slab_data(entity_handle);
    if (!slab_data)
      return false;
    return slab_data->is_coarse;
  }
  
  MOAB_RefineSlabs_inline
  RefineSlabs::SlabData::Membership RefineSlabs::get_membership( EntityHandle entity_handle )
  { 
    SlabData *slab_data = get_slab_data(entity_handle);
    if (!slab_data)
      return SlabData::EXTERNAL;
    return slab_data->membership;
  }
  MOAB_RefineSlabs_inline
  void RefineSlabs::set_membership( EntityHandle entity_handle, SlabData::Membership membership )
  { 
    SlabData *slab_data = force_slab_data(entity_handle);
    slab_data->membership = membership;
    if (membership == SlabData::EXTERNAL)
      slab_data->is_coarse = false;
  }


  MOAB_RefineSlabs_inline
  RefineSlabs::SlabData::Membership RefineSlabs::get_shrink_membership( EntityHandle entity_handle )
  { 
    SlabData *slab_data = get_slab_data(entity_handle);
    if (!slab_data)
      return SlabData::EXTERNAL;
    return slab_data->shrink_membership;
  }
  MOAB_RefineSlabs_inline
  void RefineSlabs::set_shrink_membership( EntityHandle entity_handle, SlabData::Membership membership )
  { 
    SlabData *slab_data = force_slab_data(entity_handle);
    assert( slab_data );
    slab_data->shrink_membership = membership;
    if (membership == SlabData::EXTERNAL)
    {
      slab_data->mini_me_good = false;
      slab_data->mini_me = bad_handle;
    }
  }


  // ========================================================
  // Main (non-trivial) RefineSlabs functions

  ErrorCode RefineSlabs::refine_mesh(Entities &coarse_hexes, Entities &coarse_quads, Entities &fine_hexes, Entities &fine_quads)
  {
    bool is_sweep_edge = false;
    NodePair sweep_edge;
    return refine_mesh(coarse_hexes, coarse_quads, fine_hexes, fine_quads, is_sweep_edge, sweep_edge );
  }

  RefineSlabs::RefineSlabs(Core *impl) : file_sequence(0)
  {
    assert(NULL != impl);
    mbImpl = impl;

    ErrorCode error;
    error = initialize();
    if (error != MB_SUCCESS)
    {
      std::cout<<"Error initializing RefineSlabs\n"<<std::endl;
      exit(1);
    }
  }

  RefineSlabs::~RefineSlabs()
  {
    delete ahf;
  }

  ErrorCode RefineSlabs::initialize()
  {
#ifdef USE_AHF
  std::cout << "macro USE_AHF is defined" << std::endl;
#else
  std::cout << "macro USE_AHF is undefined" << std::endl;
#endif

    ErrorCode error;
    ahf = new HalfFacetRep(mbImpl);
    if (!ahf)
      return MB_MEMORY_ALLOCATION_FAILED;

    //Check for mixed entity type
    bool chk_mixed = ahf->check_mixed_entity_type();
    if (chk_mixed)
      MB_SET_ERR(  MB_NOT_IMPLEMENTED, "Encountered a mesh with mixed entity types");

    error = ahf->initialize(); MB_CHK_ERR(error);

    //DBG
    error = ahf->print_tags(3);MB_CHK_ERR(error);

    Range _inverts, _inedges, _infaces, _incells;
    error = ahf->get_entity_ranges(_inverts, _inedges, _infaces, _incells);  MB_CHK_ERR(error);

    // Check for mixed dimensional mesh
    if ((!_inedges.empty() && !_infaces.empty()) ||(!_inedges.empty() &&  !_incells.empty()) || (!_infaces.empty() && !_incells.empty()))
      MB_SET_ERR(MB_NOT_IMPLEMENTED, "Encountered a mixed-dimensional mesh");

    // Only hexes are supported
    //Check for supported entity type
      if (!_incells.empty())
      {
        EntityType type = mbImpl->type_from_handle(_incells[0]);
        if(type != MBHEX)
          MB_SET_ERR(MB_FAILURE, "Encountered a 3D mesh that wasn't hexes, but only hexes are supported");
      }
    else if ((!_inedges.empty())||(!_infaces.empty()))
        MB_SET_ERR(MB_FAILURE, "Encountered a 1D or 2D mesh, but only hexes are supported");

    return MB_SUCCESS;
  }

  /************************************************************
   *     Interface Functions                                  *
   ************************************************************/

  ErrorCode RefineSlabs::refine_mesh(Entities &coarse_hexes, Entities &coarse_quads, Entities &fine_hexes, Entities /*&fine_quads*/, bool is_sweep_edge, NodePair &sweep_edge )
  {
    std::cout << "RefineSlabs::refine_mesh, set of " << coarse_hexes.size() << " hexes" << std::endl;
    write_file_slab( coarse_hexes, "refinement_set" );
    
    // find boundary
    // define which vertices are on the boundary of the coarse_hexes refinement set?
    ErrorCode err = MB_SUCCESS;
    Entities coarse_bndverts, coarse_opphexes;
    err = mark_hex_nodes(coarse_hexes, coarse_bndverts, coarse_opphexes);
    if (err == MB_SUCCESS)
      err = mark_surface_nodes(coarse_quads);

    // find slabs    
    // ideally two parallel sheets of hexes
    // avoid very small slabs, say with one star node, in the first pass
    // do a second pass, allowing one-star nodes
    Slabs slabs;
    if (err == MB_SUCCESS)
      err = find_slabs( coarse_hexes, coarse_quads, is_sweep_edge, sweep_edge, slabs );

    // initialize refinement
    if (err == MB_SUCCESS)
      err = initialize_refinement( coarse_hexes, coarse_quads );

    // pillow slabs
    // refine the hexes of each slab, referencing each hex's current refinement
    if (err == MB_SUCCESS)
      err = pillow_slabs( slabs, fine_hexes );

    // replace mesh
    if (err == MB_SUCCESS)
      err = replace_mesh( coarse_hexes,  coarse_bndverts, coarse_opphexes, fine_hexes);

    delete_slabs( slabs );

    return err;
  }

  void RefineSlabs::delete_slabs( Slabs &slabs )
  {
    for (size_t i = 0; i < slabs.size(); ++i )
    {
      delete slabs[i];
    }
    slabs.clear();
  }

  ErrorCode RefineSlabs::find_slabs(Entities &coarse_hexes, Entities &coarse_quads, bool is_sweep_pair, NodePair &sweep_pair, Slabs &slabs )
  {
    // Slab order strategy:
    //
    // start perpendicular to sweep direction, if specified
    //
    // search for creases
    // search for caps
    //   march outward from existing seed edges
    //
    // search for large spots
    // search for remaining small spots
    
    // rules:
    // march by 2 edges in seed edge direction, to avoid jumping around and leaving even-sized gaps of nodes that won't be pillowed


    // debug
    std::cout << "find_slabs_loc round 1, sweep directions" << std::endl;
    // todo, geometric reasoning to find these, based on caps or somesuch, rather than passed in
    if (is_sweep_pair)
    {
      SlabEdge sweep_edge;
      make_slab_edge( sweep_pair, sweep_edge );
      find_slabs_sweep( sweep_edge, slabs );
    }
    else
      std::cout << "skipping sweep" << std::endl;
    
    // debug 
    std::cout << "find_slabs_loc round 2, all edges, no tiny" << std::endl;

    ErrorCode  err = find_slabs_loc( coarse_hexes, coarse_quads,  slabs, false );

    if (err == MB_SUCCESS)
    {
      // debug 
      std::cout << "find_slabs_loc round 3, all edges, tiny ok" << std::endl;

      err = find_slabs_loc( coarse_hexes, coarse_quads,  slabs, true );
    }
    return err;
  }
  
  ErrorCode RefineSlabs::find_slabs_sweep( SlabEdge &sweep_edge, Slabs &slabs )
  {
    // debug
    size_t start_size = slabs.size();
    
    std::deque< NodePair > candidates;
    NodePair node_pair;
    make_node_pair( sweep_edge, node_pair );
    candidates.push_back( node_pair );
    
    while ( !candidates.empty() )
    {
      
      // pop
      SlabEdge slab_edge;
      make_slab_edge( candidates.front(), slab_edge );
      candidates.pop_front();

      bool is_new_slab = build_slab( slab_edge, slabs, false );
      
      // enqueue the beyond-upper edges,
      // including marching from the one used in a failed attempt to build the slab
      if ( is_new_slab )
      {
        assert( !slabs.empty() );
        Slab *new_slab = slabs.back();
        enqueue_beyond_upper(candidates, new_slab->edge_set );
      }
      else
      {
        EdgeSet edge_set;
        edge_set.insert( NodePair( slab_edge.head_node, slab_edge.tail_node ) );
        enqueue_beyond_upper(candidates, edge_set );
      }
    }
    
    // debug
    size_t end_size = slabs.size();
    std::cout << (end_size - start_size) << " sweep slabs created: ";
    if ( end_size > start_size )
      std::cout << start_size << " to " << end_size - 1;
    std::cout << std::endl;
    
    return MB_SUCCESS;
  }

  // this only makes sense in the structured sweep direction, structured map direction. It may be at odds with finding creases.
  void RefineSlabs::enqueue_beyond_upper( std::deque< NodePair > &candidates, const EdgeSet &edge_set )
  {
    for (EdgeSet::const_iterator i = edge_set.begin(); i != edge_set.end(); ++i )
    {
      SlabEdge slab_edge;
      make_slab_edge( *i, slab_edge );

      // find_beyond_upper
      SlabEdges ortho, upper, parallel, equivalent, equivalent_upper;
      get_adjacent_slab_edges( slab_edge, ortho, upper, parallel, equivalent, equivalent_upper );

      assert( upper.size() < 2 ); // remove assert if we want to experiment with unstructured sets
      for ( size_t u = 0; u < upper.size(); ++u )
      {
        SlabEdge up = upper[u];
        up.flip();
        SlabEdges uortho, uupper, uparallel, uequivalent, uequivalent_upper;
        get_adjacent_slab_edges( up, uortho, uupper, uparallel, uequivalent, uequivalent_upper );
        assert( uupper.size() < 2 ); // remove assert if we want to experiment with unstructured sets
        for ( size_t uu = 0; uu < uupper.size(); ++uu )
        {
          SlabEdge upup = uupper[uu];
          upup.flip();
          NodePair upup_pair;
          make_node_pair( upup, upup_pair );
          candidates.push_back( upup_pair );
        }
      }
    }
  }
  
  ErrorCode RefineSlabs::initialize_refinement( Entities &coarse_hexes, Entities &/*coarse_quads*/ )
  {

    // estimate memory size for AHF
    // this could be made tighter by considering the quads, and the actual slabs found.
    size_t memory_estimate = 8 * coarse_hexes.size() + 8 * coarse_hexes.size();
    //                        refine each hex into 8      transition == outer pillow

    // make a new AHF to store the refinement
    ErrorCode error = estimate_allocate_ahf_maps( memory_estimate );
    assert( error == MB_SUCCESS );

    // copy nodes
    // each node and its copy point to each other
    copy_nodes( coarse_hexes );

    // initial refined hexes = original hexes
    // copy hexes
    // a hex points to its refinement (other way not needed?)
    Entities dup_hexes;
    error = copy_hexes( coarse_hexes, dup_hexes);MB_CHK_ERR(error);

    // ? anything needed for the quads at this point?
    // should we copy them?

    //update ahf maps for the copied entities
    error = update_ahf_maps(dup_hexes);MB_CHK_ERR(error);
    return error;

  }

  ErrorCode RefineSlabs::copy_nodes( Entities &coarse_hexes )
  {
    for (size_t i = 0; i < coarse_hexes.size(); ++i )
    {
      EntityHandle chex = coarse_hexes[i];
      EntityHandle coarse_nodes [8];
      // get 8 coarse nodes of hex
      get_hex_nodes( chex, coarse_nodes);

      for (size_t j = 0; j < 8; ++j)
      {
        // create a new node with the same coordinates as the coarse one
        // associate them with each other
        copy_node( coarse_nodes[j] );          
      }
    }
    return MB_SUCCESS;
  }

  void RefineSlabs::get_copy_nodes( EntityHandle *coarse_nodes, EntityHandle *fine_nodes, int num_nodes)
  {
    for (int i = 0; i < num_nodes; ++i)
    {
      bool success = get_copy( coarse_nodes[i], fine_nodes[i] );      
      assert( success );
      if (!success)
        return;
    }
  }


  ErrorCode RefineSlabs::copy_hexes( Entities &coarse_hexes, Entities &fine_hexes )
  {
    fine_hexes.clear();
    for (size_t i = 0; i < coarse_hexes.size(); ++i )
    {
      EntityHandle chex = coarse_hexes[i];
      EntityHandle coarse_nodes [8], fine_nodes[8]; 
      // todo next level of sophistication
      // generalize to 2nd order hexes with more nodes, etc
      get_hex_nodes( chex, coarse_nodes);
      get_copy_nodes( coarse_nodes, fine_nodes, 8);
      EntityHandle fhex;
      create_hex( fine_nodes, fhex );
      add_refined_hex( chex, fhex );
      fine_hexes.push_back(fhex);
    }
    return MB_SUCCESS;
  }
  
  bool RefineSlabs::build_slab( SlabEdge &slab_edge, Slabs &slabs,  bool tiny_slabs_OK )
  {
    // recursively extend the slab
    // the passed in slab_edge starts the queue
    Slab* slab = new Slab();
    slab->reset();
    
    extend_slab( slab_edge, *slab);
    
    finalize_slab( *slab, tiny_slabs_OK );
    
    // debug
    std::cout << "slab file # " << slabs.size() + 2 << std::endl;
    print_slab( *slab );
    Entities hexes;
    get_pillow_hexes( *slab, hexes );
    write_file_slab( hexes, "slab", hexes.size() );
    
    
    if ( slab->empty() )
    {
      delete slab;
      return false;
    }
    //else
    slabs.push_back( slab );
    return true;
  }


  ErrorCode RefineSlabs::find_slabs_loc( Entities &coarse_hexes, Entities &/*coarse_quads*/, Slabs &slabs, bool tiny_slabs_OK )
  {
    size_t start_size = slabs.size();
    
    for ( size_t c = 0; c < coarse_hexes.size(); ++c )
    {
      EntityHandle hex = coarse_hexes[c];
      SlabEdge slab_edge;
      for (int edge_lid = 0; edge_lid < 12; ++edge_lid)
      {
        // find a seed edge, then recursively extend
        //   find_seed_edge may increment edge_lid until a good one is found, or 12 is reached
        if ( find_seed_edge( hex, edge_lid, slab_edge ) )
        {
          build_slab( slab_edge, slabs, tiny_slabs_OK );
        }
      }
    }

    std::cout << slabs.size();
    std::cout << " slabs total." << std::endl;
    std::cout << (slabs.size() - start_size) << " new slabs, ";
    if (tiny_slabs_OK)
      std::cout << " some could be tiny, a single star node. " << std::endl;
    else
      std::cout << " all big > 1 star node." << std::endl;

    return MB_SUCCESS;
  }

  void RefineSlabs::make_node_pair( const SlabEdge &slab_edge, NodePair &node_pair )
  {
    node_pair.first = slab_edge.head_node;
    node_pair.second = slab_edge.tail_node;
  }

  void RefineSlabs::make_slab_edge( const NodePair &node_pair, SlabEdge &slab_edge ) 
  {
    slab_edge.head_node = node_pair.first;
    slab_edge.tail_node = node_pair.second;
    Entities star_hexes;
    get_all_hexes( slab_edge.head_node, star_hexes, true );
    // find one that contains the tail_node
    for ( size_t h = 0; h < star_hexes.size(); ++h )
    {
      slab_edge.hex = star_hexes[h];
      // we don't need to check if the hex being in the set, since the head was a potential star node, all of its hexes must be in the set
      const int tail_index = get_hex_node_index( slab_edge.hex, slab_edge.tail_node );
      if ( tail_index < 0 )
        continue;
      const int head_index = get_hex_node_index( slab_edge.hex, slab_edge.head_node );
      assert( head_index >= 0 && head_index < 8 );

      slab_edge.edge_lid = get_hex_edge_index ( head_index, tail_index );
      return;
    }
    assert(0);
  }

  void RefineSlabs::finalize_slab( Slab &slab, bool tiny_slabs_OK )
  {
    // if any non-internal ortho edge is already refined, take its hexes out of the set
    // beware creating strange, non-convex sets
    SlabEdges ortho, upper, parallel, equivalent, equivalent_upper;
    bool recheck = true;
    while ( recheck )
    {
      recheck=false;
      for ( EdgeSetIterator e = slab.edge_set.begin(); e != slab.edge_set.end(); /*no increment*/ )
      {
        EntityHandle head_node = e->first;
        assert( slab.get_is_star( head_node ) );
        SlabEdge slab_edge;
        make_slab_edge( *e, slab_edge );
        get_adjacent_slab_edges( slab_edge, ortho, upper, parallel, equivalent, equivalent_upper );
        bool remove_me = false;
        {
          for ( SlabEdges::iterator o = ortho.begin(); o != ortho.end(); ++o )
          {
            EntityHandle other_node = o->tail_node;
            assert( other_node != head_node && o->head_node == head_node );
            // debug
            // std::cout << " is star " << slab.get_is_star( other_node ) << " refinement level " << get_edge_refinement_level( *o ) << std::endl;
            
            if ( !slab.get_is_star( other_node ) && get_edge_refinement_level( *o ) > 1) //zzyk should this be 0 or 1?
            {
              remove_me = true;
              recheck = true;
              break;
            }
          }
        }
        if (remove_me)
        {
          slab.edge_set.erase( e++ );
          slab.set_is_star( head_node, false );
          // decrement edge refinement levels
          decrement_edge_refinement_level( equivalent );
          decrement_edge_refinement_level( equivalent_upper );

          // I haven't incremented the ortho edge refinement levels, so I don't need to increment them
        }
        else
          e++;
      }
    }

    // additional checks
    bool delete_me = false; 

    // avoid very small slabs, say with one star node
    if ( !tiny_slabs_OK && slab.star_nodes.size() < 2 )
    {
      delete_me = true;
    }

    // discard the slab
    if ( delete_me )
    {
      // delete the whole thing
      for ( EdgeSetIterator e = slab.edge_set.begin(); e != slab.edge_set.end(); /*no increment*/ )
      {
        EntityHandle head_node = e->first;
        assert( slab.get_is_star( head_node ) );
        SlabEdge slab_edge;
        make_slab_edge( *e, slab_edge );
        get_adjacent_slab_edges( slab_edge, ortho, upper, parallel, equivalent, equivalent_upper );
        slab.edge_set.erase( e++ );
        slab.set_is_star( head_node, false );
        // decrement edge refinement levels
        decrement_edge_refinement_level( equivalent );
        decrement_edge_refinement_level( equivalent_upper );
      }
      return;
    }

          // I haven't incremented the ortho edge refinement levels, so I don't need to increment them

    // increment the refinement level of the non-internal ortho edges
    {
      for ( EdgeSetIterator e = slab.edge_set.begin(); e != slab.edge_set.end(); ++e )
      {
        EntityHandle head_node = e->first;
        assert( slab.get_is_star( head_node ) );
        SlabEdge slab_edge;
        make_slab_edge( *e, slab_edge );
        get_adjacent_slab_edges( slab_edge, ortho, upper, parallel, equivalent, equivalent_upper );
        for ( SlabEdges::iterator o = ortho.begin(); o != ortho.end(); ++o )
        {
          EntityHandle other_node = o->tail_node;
          assert( other_node != head_node && o->head_node == head_node );
          if ( !slab.get_is_star( other_node ) )
          {
            SlabEdges equivalent_ortho;
            find_equivalent_edges( *o, equivalent_ortho );
            increment_edge_refinement_level( equivalent_ortho );
          }
        }
      }
    }

  }


  void RefineSlabs::extend_slab( SlabEdge &slab_edge, Slab &slab )
  {
    // enqueue the passed in slab_edge
    // while the queue is not empty
    //   pop edge
    //   if the edge is good, then
    //      mark it as being in the set 
    //      enqueue the parallel edges 

    // enqueue passed-in edge
    SlabEdges edge_queue; //todo replace with deque ?
    edge_queue.push_back(slab_edge);

    SlabEdges ortho, upper, parallel, equivalent_edges, equivalent_upper;

    // while the queue is not empty
    // note: edge_queue.size() grows inside the loop
    for( size_t i = 0; i < edge_queue.size(); ++i )
    {
      SlabEdge se ( edge_queue[i] ); // copy constructor

      if ( is_good_slab_edge(se, slab, ortho, upper, parallel, equivalent_edges, equivalent_upper) )
      {
        add_edge( se, slab, equivalent_edges, equivalent_upper );

        // enqueue the parallel edges
        // grow in bfs order, not dfs, in order to get better geometrically shaped slabs
        edge_queue.insert( edge_queue.end(), parallel.begin(), parallel.end() );
      }
    }
  }

  void RefineSlabs::add_edge( SlabEdge &slab_edge, Slab &slab, SlabEdges &equivalent_edges, SlabEdges &equivalent_upper )
  {
    // for the head vertex
    // get all the (coarse) hexes containing the vertex
    // add them to the slab
    slab.add( slab_edge );

    // increment the refinement level of this and all the equivalent edges
    // DO NOT increment_edge_refinement_level(slab_edge); this is already included in the vector below
    increment_edge_refinement_level( equivalent_edges );

    // increment the refinement level of the equivalent upper edges
    increment_edge_refinement_level( equivalent_upper );
    
    // At the end, when the slab is built, we will increment adjacent edges that don't have *both* nodes in the set.
    // that can lead to such edges being refined twice
    // we need to back off if that causes an edge to be refined three times

  }

  void RefineSlabs::Slab::add( SlabEdge &slab_edge )
  {
    set_in_slab( slab_edge, true);
    EntityHandle star_node = slab_edge.head_node;
    set_is_star( star_node, true );
  }

  bool RefineSlabs::Slab::empty() const
  {
    assert( star_nodes.size() == edge_set.size() );
    return star_nodes.empty();
  }


  bool RefineSlabs::find_seed_edge( EntityHandle hex, int &edge_lid, SlabEdge &slab_edge )
  {
    // caller sets initial value of edge_lid
    for ( /*caller initializes edge_lid*/; edge_lid < 12; ++edge_lid)
    {
      if ( hex == EntityHandle(10376293541461622811ul) && edge_lid == 2 )
        std::cout << " special edge " << std::endl;
      
      bool good_edge = 
        is_good_seed_edge( hex, edge_lid, 0, slab_edge ) ||
        is_good_seed_edge( hex, edge_lid, 1, slab_edge );
      if( good_edge )
        return true;
    }
    return false;
  }

  bool RefineSlabs::is_good_slab_edge( const SlabEdge &slab_edge, Slab &slab )
  {
    SlabEdges ortho, upper, parallel, equivalent, equivalent_upper;
    return is_good_slab_edge( slab_edge, slab, ortho, upper, parallel, equivalent, equivalent_upper );
  }
  
  bool RefineSlabs::is_good_slab_edge( const SlabEdge &slab_edge, Slab &slab, std::vector<SlabEdge> &ortho, std::vector<SlabEdge> &upper, 
                                       std::vector<SlabEdge> &parallel, std::vector<SlabEdge> &equivalent, 
                                       std::vector<SlabEdge> &equivalent_upper )
  {
    ortho.clear(); 
    upper.clear();  
    parallel.clear(); 
    equivalent.clear();
    equivalent_upper.clear();

    if 
    ( 
      // skip if head is already in the current set
      slab.get_is_star( slab_edge.head_node) || 

      // skip if we've a twisted sheet and the tail is a star node
      slab.get_is_star( slab_edge.tail_node ) ||

      // skip if head is on the boundary of the refinement set (or outside, but it shouldn't be outside)
      get_membership( slab_edge.head_node ) != SlabData::INTERNAL ||

      // skip if edge has already been refined,
      get_edge_refinement_level(slab_edge) > 0 || 

      // skip if vertex is on the boundary of the geometric volume
      get_geometry_dimension( slab_edge.head_node ) < 3
    )
      return false;

    get_adjacent_slab_edges( slab_edge, ortho, upper, parallel, equivalent, equivalent_upper );
    // ensure no non-ortho edge is already refined
    // ensure upper and ortho edges are not already in the current set!
    if (
         !none_refined( upper ) || 
         slab.any_in_slab(ortho) || slab.any_in_slab(upper)
       )
    {
      ortho.clear(); 
      upper.clear();  
      parallel.clear(); 
      equivalent.clear();
      equivalent_upper.clear();
      return false;
    }

    return true;
  }
  void RefineSlabs::find_equivalent_edges( SlabEdge &slab_edge, std::vector<SlabEdge> & equivalent_edges )
  {
    SlabEdges ortho, upper, parallel, equivalent_upper;
    get_adjacent_slab_edges( slab_edge, ortho, upper, parallel, equivalent_edges, equivalent_upper );
  }


  bool RefineSlabs::is_good_seed_edge( EntityHandle hex, int edge_lid, int node_01, SlabEdge &slab_edge )
  {
    // make SlabEdge from one of the two vertices
    get_edge( hex, edge_lid, node_01, slab_edge );
    Slab slab;

    return is_good_slab_edge( slab_edge, slab );
  }

  bool RefineSlabs::none_refined( std::vector<SlabEdge> slab_edges )
  {
    for (size_t i = 0; i < slab_edges.size(); ++i )
    {
      const SlabEdge &slab_edge = slab_edges[i];
      if (get_edge_refinement_level(slab_edge) > 0)
        return false;
    }
    return true;
  }

  void RefineSlabs::get_adjacent_slab_edges( const SlabEdge &slab_edge, SlabEdges &ortho, 
        SlabEdges &upper, SlabEdges &parallel,
        SlabEdges &equivalent,  SlabEdges &equivalent_upper )
  {
    ortho.clear(); 
    upper.clear();  
    parallel.clear(); 
    equivalent.clear();
    equivalent_upper.clear();

    Entities hexes;
    get_all_hexes( slab_edge.head_node, hexes, true );
    Entities non_sheet_hexes;
    for (size_t h = 0; h < hexes.size(); ++h )
    {
      EntityHandle hex = hexes[h];      
      SlabEdge match;
      if( get_matching_edge( hex, slab_edge, match ) )
      {
        // make sure the returned slab edge is well defined
        assert( match.edge_lid >= 0 );
        assert( match.edge_lid < 12 );
        
        SlabEdge adj1, adj2, opp1, opp2;
        get_adj( match, adj1, adj2 );
        assert( adj1.edge_lid > -1 );
        assert( adj2.edge_lid > -1 );
        add_unique( ortho, adj1 );
        add_unique( ortho, adj2 );
        get_opp( match, opp1, opp2 );
        add_unique( parallel, opp1 );
        add_unique( parallel, opp2 );

        equivalent.push_back( match );
      }
      else
        non_sheet_hexes.push_back( hex );
    }
    for (size_t h = 0; h < non_sheet_hexes.size(); ++h )
    {
      // the ortho slab edges are the edges of the upper hexes containing the head vertex, that are not orthogonal to it
      // e.g. in a structured mesh, there should be only one such edge
      EntityHandle hex = non_sheet_hexes[h];
      //get_matching_node( slab_edge.head_node, hex, node_lid );
      SlabEdge star[3];
      get_star_edges( hex, slab_edge.head_node, star );
      for ( size_t e = 0; e < 3; ++e )
        if ( unique( ortho, star[e] ) ) //i.e. if its not an ortho edge
        {
          add_unique( upper, star[e] );
          equivalent_upper.push_back( star[e] );
        }
    }
  }

  int RefineSlabs::get_hex_node_index( EntityHandle hex, EntityHandle node )
  {
    EntityHandle hex_nodes[8];
    get_hex_nodes( hex, hex_nodes);
    for ( int n = 0; n < 12; ++n )
      if ( hex_nodes[n] == node )
        return n;
    return -1;
  }

  void RefineSlabs::get_star_edges( EntityHandle hex, EntityHandle node, SlabEdge star[3] )
  {
    int i = 0;
    for (size_t e = 0; e < 12; ++e )
    {
      // array range check
      assert( i < 3 );
      assert( i >= 0 );

      get_edge( hex, e, 0, star[i] );
      if ( star[i].head_node == node )
      {
        // match
        assert( i < 3 );
        ++i;
      }
      else if ( star[i].tail_node == node )
      {
        // match flipped
        assert( i < 3 );
        star[i].flip();
        ++i;
      }
      // all done if we've found all three
      if ( i == 3 )
        return;
    }
    // error if we didn't find all three
    assert( 0 );
  }


  bool RefineSlabs::get_matching_edge( EntityHandle hex, const SlabEdge &slab_edge, SlabEdge &match )
  {
    for (int e = 0; e < 12; ++e)
    {
      get_edge( hex, e, 0, match );
      if ( slab_edge.nodes_match(match) )
      { 
        if( !slab_edge.directions_match(match) )
        {
          match.flip();
        }
        return true;
      }
    }
    match.edge_lid = -1;
    match.hex = bad_handle;
    match.head_node = bad_handle;
    match.tail_node = bad_handle;
    return false;
  }
  bool RefineSlabs::unique( const SlabEdges &edges, const SlabEdge &edge )
  {
    for (size_t i = 0; i < edges.size(); ++i )
    {
      if ( edges[i].nodes_match( edge ))
      {
        // should we also require directions to match?
        if ( !edges[i].directions_match( edge ) )
        {
          // this shouldn't happen, because the caller is careful to orient them the same way
          std::cout << "SlabEdge nodes match, but not their directions" << std::endl; 
          assert(0);
        }
        return false;
      }
    }
   return true;
  }
  bool RefineSlabs::add_unique( SlabEdges &edges, const SlabEdge &edge )
  {
    if (unique( edges, edge ))
    {
      edges.push_back( edge );
      return true;
    }
    return false;
  }
  
  ErrorCode RefineSlabs::pillow_slabs( Slabs &slabs, Entities &new_hexes)
  {
    for (size_t i = 0; i < slabs.size(); ++i)
    {
      ErrorCode err_i = pillow_slab( *slabs[i], new_hexes );
      if (err_i != MB_SUCCESS)
        return err_i;
    }
    return MB_SUCCESS;
  }

  void RefineSlabs::get_pillow_hexes( const Slab& slab, Entities &pillowhexes )
  {
    Entities star_hexes;
    for ( EntitySet::const_iterator n = slab.star_nodes.begin(); n != slab.star_nodes.end(); ++n )
    {
      EntityHandle star_node = *n;
      get_all_hexes( star_node, star_hexes, true );
      pillowhexes.insert( pillowhexes.end(), star_hexes.begin(), star_hexes.end() );
    }

    // uniquify the slab, removing redundant hexes
    std::sort( pillowhexes.begin(), pillowhexes.end() );
    pillowhexes.erase( std::unique( pillowhexes.begin(), pillowhexes.end() ), pillowhexes.end() );

  }

  void RefineSlabs::shrink_mark_slab( Entities &slab, bool iscoarse )
  { 
    for (size_t i = 0; i < slab.size(); ++i)
    {
      EntityHandle hex = slab[i];
      set_shrink_membership( hex, SlabData::INTERNAL );
    }

    // set internal and boundary nodes of the shrink set
    for (size_t i = 0; i < slab.size(); ++i)
    {
      EntityHandle chex = slab[i];
      EntityHandle hex_nodes[8];
      get_hex_nodes( chex, hex_nodes );
      for ( size_t n = 0; n < 8; ++n )
      {
        EntityHandle hex_node = hex_nodes[n];

        // only check/set it once
        if ( get_shrink_membership(hex_node) == SlabData::EXTERNAL )
        {
          bool is_internal = true;

          // check if its internal to the refinement set
          if (get_membership(hex_node) != SlabData::INTERNAL ) 
            is_internal = false;

          // check if all its hexes are in the set
          if (is_internal)
          {
            Entities node_hexes;
            get_all_hexes( hex_node, node_hexes, iscoarse );
            for ( size_t k = 0; k < node_hexes.size(); ++k )
            {
              EntityHandle node_hex = node_hexes[k];
              if ( get_shrink_membership( node_hex ) == SlabData::EXTERNAL )
              {
                is_internal = false;
                break;
              }
            }
          }
          set_shrink_membership( hex_node, is_internal ? SlabData::INTERNAL : SlabData::BOUNDARY );
        }
      }
    }
  }

  void RefineSlabs::remove_shrink_mark_slab( Entities &slab, bool /*is_coarse*/ )
  {
    for (size_t i = 0; i < slab.size(); ++i)
    {
      EntityHandle chex = slab[i];
      set_shrink_membership( chex, SlabData::EXTERNAL );
      EntityHandle hex_nodes[8];
      get_hex_nodes( chex, hex_nodes );
      for ( size_t n = 0; n < 8; ++n )
      {
        set_shrink_membership( hex_nodes[n], SlabData::EXTERNAL );
      }
    }
  }

  void RefineSlabs::shrink_mark_coarse_slab( Entities &slab )
  { 
    shrink_mark_slab( slab, true );
  }

  void RefineSlabs::shrink_mark_fine_slab( Entities &slab, Entities &shrink_set)
  {
    // todo next level of sophistication
    // the tricky part is tracking the edges and quads of the fine mesh that are on the coarse mesh quad and hex
    // for now, just put *all* the fine hexes into the set

    // choose the set of hexes

    // for every coarse hex
    for (size_t i = 0; i < slab.size(); ++i)
    {
      EntityHandle chex = slab[i];
      SlabData *coarse_data = get_slab_data( chex );
      assert( coarse_data );      
      HexRefinement *crefine = coarse_data->refinement;
      assert( crefine );
      Entities &fine_hexes = crefine->fine_hexes; // alias

      // for all its fine hexes
      for (size_t j = 0; j < fine_hexes.size(); ++j)
      {
        EntityHandle fine_hex = fine_hexes[j];

        // todo next level of sophistication
        // additional containment tests here
        // check the coarse_owner of every fine hex, see if that entity is internal or boundary to the set

        shrink_set.push_back( fine_hex );
      }
    }

    // mark the nodes for boundary or internal, as before
    shrink_mark_slab( shrink_set, false );
  }

  void RefineSlabs::pillow_hexes( Entities &shrink_set, Entities &new_hexes )
  {
    Entities finehexes;

    // shrink set and new hexes are both fine-mesh entities
    // split vertices into edges
    // gather boundary entities
    for (size_t i = 0; i < shrink_set.size(); ++i)
    {
      EntityHandle chex = shrink_set[i];
      EntityHandle hex_nodes[8];
      // these are fine hexes and nodes
      get_hex_nodes( chex, hex_nodes ); 
      for ( size_t n = 0; n < 8; ++n )
      {
        EntityHandle node = hex_nodes[n];        
        EntityHandle shrunk_node; 
        if ( get_shrink_membership( node ) == SlabData::BOUNDARY )
        {
          if ( !get_shrunk_node( node, shrunk_node ) ) 
          {
            shrunk_node = shrink_node( node );
          // todo next level of sophistication
          // pick some new geometry for the copied node, say partway towards the hex center
            move_into_shrinkset( node );
          }
          assert( shrunk_node );
          replace_node( chex, n, shrunk_node);
        }
      } 

      // split edges into quads
      // not needed explicitly

      // split quads into hexes
      for (size_t q = 0; q < 6; ++q)
      {
        // get quad nodes, with outward facing normal
        EntityHandle quad_nodes[8];
        get_quad_nodes( chex, hex_nodes, q, quad_nodes );
        bool make_hex = true;
        for ( size_t n = 0; n < 4; n++ )
        {
          EntityHandle shrunk_node;
          if ( get_shrunk_node( quad_nodes[n], shrunk_node ) )
            quad_nodes[n+4] = shrunk_node;
          else 
          {
            make_hex = false;
            break;
          }
        }
        // if all nodes had a copy, make a hex out of the copy and the original
        if ( make_hex )
        {
          EntityHandle fhex;
          create_hex( quad_nodes, fhex );
          add_refined_hex(chex, fhex);
          finehexes.push_back(fhex);
          new_hexes.push_back(fhex);
        }
      }
    }

    // update the ahf maps for the new hexes
    //DBG
    std::cout<<"Shrink set and new hexes"<<std::endl;
    for (int i=0; i<shrink_set.size(); i++)
      {
        EntityHandle hex_nodes[8];
        get_hex_nodes( shrink_set[i], hex_nodes );
        std::cout<<"shrinkset["<<i<<"] = "<<shrink_set[i]<<", conn = ["<<hex_nodes[0]<<", "<<hex_nodes[1]<<", "<<hex_nodes[2]<<", "<<hex_nodes[3]<<", "<<hex_nodes[4]<<", "<<hex_nodes[5]<<", "<<hex_nodes[6]<<", "<<hex_nodes[7]<<" ]"<<std::endl;
      }

    for (int i=0; i<finehexes.size(); i++)
      {
        EntityHandle hex_nodes[8];
        get_hex_nodes(finehexes[i], hex_nodes );
        std::cout<<"finehexes["<<i<<"] = "<<finehexes[i]<<", conn = ["<<hex_nodes[0]<<", "<<hex_nodes[1]<<", "<<hex_nodes[2]<<", "<<hex_nodes[3]<<", "<<hex_nodes[4]<<", "<<hex_nodes[5]<<", "<<hex_nodes[6]<<", "<<hex_nodes[7]<<" ]"<<std::endl;
      }

    ErrorCode err = update_ahf_maps(shrink_set, finehexes);
    assert(err == MB_SUCCESS);
  }
  void RefineSlabs::get_quad_nodes( EntityHandle /*hex*/, const EntityHandle hex_nodes[8], int face_lid, EntityHandle* quad_nodes )
  {
    int face_node_ids[4];
    int num_face_nodes = 4;
    // see CN.hpp and MBEntityType.h
    EntityType quad_type;
    CN::SubEntityNodeIndices( moab::MBHEX, 8, 2, face_lid, quad_type, num_face_nodes, face_node_ids);
    assert( quad_type == moab::MBQUAD );
    assert( num_face_nodes >= 4 );
    for (size_t i = 0; i < 4; ++i)
      quad_nodes[i] = hex_nodes[ face_node_ids[i] ];
  }


  ErrorCode RefineSlabs::get_quad_nodes( EntityHandle quad, EntityHandle quad_nodes[4] )
  {
    int num_nodes = 0;
    const bool corners_only = true;
    const EntityHandle *const_quad_nodes;
    ErrorCode error = mbImpl->get_connectivity(quad, const_quad_nodes, num_nodes, corners_only ); MB_CHK_ERR(error);
    assert(num_nodes == 4);

    // transfer and cast away const
    for (size_t i = 0; i < 4; ++i)
    {
      quad_nodes[i] = *const_cast<EntityHandle*>(&const_quad_nodes[i]);
    }

    return MB_SUCCESS;
  }


  ErrorCode RefineSlabs::get_hex_nodes( EntityHandle hex, EntityHandle hex_nodes[8] )
  {
    // debug
    hex_ok( hex );
    
    int num_nodes = 0;
    const bool corners_only = true;
    const EntityHandle *const_hex_nodes;
    ErrorCode error = mbImpl->get_connectivity(hex, const_hex_nodes, num_nodes, corners_only ); MB_CHK_ERR(error);
    assert( num_nodes == 8 );

    // transfer and cast away const
    for (size_t i = 0; i < 8; ++i)
    {
      hex_nodes[i] = *const_cast<EntityHandle*>(&const_hex_nodes[i]);
    }

    return MB_SUCCESS;
  }

  void RefineSlabs::get_edge( EntityHandle hex, int edge_lid, int node_01, SlabEdge &slab_edge )
  {
    // get handles of nodes in the hex, correctly ordered
    EntityHandle hex_nodes [8];
    get_hex_nodes( hex, hex_nodes );

    // CN.hpp
    // Use the cannonical connectivity of a hex to get the indices of the vertices in the edge_lid
    EntityType edge_type;
    int edge_nodes[2];
    int num_edge_nodes;
    CN::SubEntityNodeIndices( moab::MBHEX, 8, 1, edge_lid, edge_type, num_edge_nodes, edge_nodes );
    assert( edge_type == moab::MBEDGE );
    assert( num_edge_nodes == 2 );
    slab_edge.hex = hex;
    slab_edge.edge_lid = edge_lid;
    if (node_01 == 0)
    {
      slab_edge.head_node = hex_nodes[ edge_nodes[0] ];
      slab_edge.tail_node = hex_nodes[ edge_nodes[1] ];
    }
    else
    {
      slab_edge.head_node = hex_nodes[ edge_nodes[1] ];
      slab_edge.tail_node = hex_nodes[ edge_nodes[0] ];
    }
  }

  bool RefineSlabs::is_equal( const SlabEntity *m1, const SlabEntity *m2 ) const
  {
    // relies on these type overloading operator == 

    if ( dynamic_cast<const SlabHex*>(m1) && dynamic_cast<const SlabHex*>(m2) )
    {
      return *dynamic_cast<const SlabHex*>(m1) == *dynamic_cast<const SlabHex*>(m2);
    }
    // unused
    // if ( dynamic_cast<const SlabQuad*>(m1) && dynamic_cast<const SlabQuad*>(m2) )
    // {
    //   return *dynamic_cast<const SlabQuad*>(m1) == *dynamic_cast<const SlabQuad*>(m2);
    // }
    if ( dynamic_cast<const SlabEdge*>(m1) && dynamic_cast<const SlabEdge*>(m2) )
    {
      return *dynamic_cast<const SlabEdge*>(m1) == *dynamic_cast<const SlabEdge*>(m2);
    }
    if ( dynamic_cast<const SlabNode*>(m1) && dynamic_cast<const SlabNode*>(m2) )
    {
      return *dynamic_cast<const SlabNode*>(m1) == *dynamic_cast<const SlabNode*>(m2);
    }
    // BadEntity or different classes
    return false;
  }


  ErrorCode RefineSlabs::pillow_slab( Slab &slab, Entities &new_hexes )
  {

    Entities coarse_hexes;
    get_pillow_hexes( slab, coarse_hexes );

    // mark each hex as being in the set    
    shrink_mark_coarse_slab( coarse_hexes );

    // convert coarse set to fine set of hexes to pillow
    // find internal and boundary fine hexes
    Entities shrink_set;
    shrink_mark_fine_slab( coarse_hexes, shrink_set);

    // pillow, splitting vertices, edges, then quads into hexes
   // Entities new_hexes;
    pillow_hexes( shrink_set, new_hexes );

    // cleanup 
    remove_shrink_mark_slab( coarse_hexes, true );
    remove_shrink_mark_slab( shrink_set, false );
    remove_shrink_mark_slab( new_hexes, false );

    return MB_SUCCESS;
  }

  ErrorCode RefineSlabs::mark_hex_nodes(Entities &coarse_hexes, Entities &coarse_bndverts, Entities &coarse_opphexes)
  {
    ErrorCode error;

    // mark all hexes as being in the set
    for (size_t i = 0; i < coarse_hexes.size(); ++i)
    {
      EntityHandle hex = coarse_hexes[i];
      set_coarse_hex( hex );
    }

    // mark selected nodes as being on the boundary of the refinement hex set
    // a node is internal iff all its hexes are in the set, and the node is internal to the geometric volume
    // geometric-surface nodes are marked as on the boundary in mark_surface_nodes
    for (size_t i = 0; i < coarse_hexes.size(); ++i)
    {
      EntityHandle hex = coarse_hexes[i];
      EntityHandle hex_nodes [8];
      error = get_hex_nodes( hex, hex_nodes );MB_CHK_ERR(error);
      for (size_t j = 0; j < 8; ++j)
      {
        EntityHandle node = hex_nodes[j];
        // process each node only once
        if ( !is_coarse(node) )
        {
          bool is_internal(true);

          // check if internal to the geometry
          if ( get_geometry_dimension(node) < 3)
            is_internal = false;

          // check if internal to the set of hexes
          if (is_internal)
          {
            Entities hexes;
            get_all_hexes( node, hexes, true );
            for ( size_t k = 0; k < hexes.size(); ++k )
            {
              EntityHandle h = hexes[k];
              if ( !is_coarse(h) )
              {
                is_internal = false;
                break;
              }
            }
          }

          set_coarse_node( node );
          set_membership( node, is_internal ? SlabData::INTERNAL : SlabData::BOUNDARY );
          if (!is_internal)
            coarse_bndverts.push_back(node);
        }
      }
    }

    //Return the coarse set boundary vertices and the opposite hexes incident on the coarse set boundary hexes.
    error = get_coarse_opp_hexes(coarse_hexes, coarse_opphexes);MB_CHK_ERR(error);

    return MB_SUCCESS;
  }

  ErrorCode RefineSlabs::mark_surface_nodes(Entities &coarse_quads)
  {
    // mark all the quads as being on the surface == a coarse entity
    for (size_t i = 0; i < coarse_quads.size(); ++i)
    {
      EntityHandle quad = coarse_quads[i];
      set_coarse_quad( quad );
    }

    // mark selected nodes as being internal to the set
    // a node is internal iff all of quads that are owned by a geometric surface 
    // are in the coarse_quads set
    for (size_t i = 0; i < coarse_quads.size(); ++i)
    {
      EntityHandle quad = coarse_quads[i];
      EntityHandle quad_nodes [4];
      get_quad_nodes( quad, quad_nodes );
      for (size_t j = 0; j < 4; ++j)
      {
        EntityHandle node = quad_nodes[j];
        SlabData::Membership membership = get_membership( node );
        // its boundary if it was on the geometric surface
        if ( ( membership == SlabData::BOUNDARY ) && ( get_geometry_dimension(node) < 3 ) )
        {
          Entities quads;
          get_all_quads( node, quads );
          bool is_internal(true);
          for ( size_t k = 0; k < quads.size(); ++k )
          {
            EntityHandle q = quads[k];
            if ( (get_geometry_dimension(q)==2) && !is_coarse(q))
            {
              // some of its faces are not in the set, so consider it boundary
              is_internal = false;
              break;
            }
          }
          set_membership( node, is_internal ? SlabData::INTERNAL : SlabData::BOUNDARY );
        }
      }
    }
    return MB_SUCCESS;
  }

  EntityHandle RefineSlabs::copy_node( EntityHandle coarse_node )
  {
    SlabData *slab_data = force_slab_data(coarse_node);
    if (!slab_data->my_copy_good)
    {
      EntityHandle fine_copy;
      create_node( coarse_node, fine_copy );
      set_copy( coarse_node, fine_copy ); // sets slab_data->my_copy
      SlabData *copy_slab = force_slab_data(fine_copy);
      copy_slab->copy_data( slab_data );
      copy_slab->is_coarse = false;
    }
    assert( slab_data->is_coarse );
    return slab_data->my_copy;
  }

  EntityHandle RefineSlabs::shrink_node( EntityHandle fine_node )
  {
    SlabData *slab_data = force_slab_data(fine_node);
    assert( slab_data->is_coarse == false );
    EntityHandle fine_copy;
    create_node( fine_node, fine_copy );
    set_shrunk_node( fine_node, fine_copy );
    SlabData *copy_slab = force_slab_data(fine_copy);
    copy_slab->copy_data( slab_data );
    copy_slab->is_coarse = false;

    copy_slab->coarsening = slab_data->coarsening;

    // the coarse mesh should still point to the original fine_node, not its copy
    return fine_copy;
  }

  // coord = 0
  void RefineSlabs::zero_coord( double *coord )
  {
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 0;
  }
  // coord += add
  void RefineSlabs::add_coord( double *coord, double *add )
  {
    coord[0] += add[0];
    coord[1] += add[1];
    coord[2] += add[2];
  }
  // c = a*x + b*y
  void RefineSlabs::axpby_coord( double *coord, double a, double *x, double b, double *y )
  {
    coord[0] += a * x[0] + b * y[0];
    coord[1] += a * x[1] + b * y[1];
    coord[2] += a * x[2] + b * y[2];
  }

  void RefineSlabs::move_into_shrinkset( EntityHandle node )
  {
    // get the hexes around the node
    // get an average interior and exterior direction, towards the centers of hexes
    // Use the interior or the negative of the exterior, depending on which subset had fewer hexes.
    Entities hexes;
    get_all_hexes( node, hexes, false );
    size_t num_int(0), num_ext(0);
    double ave_int[3], ave_ext[3];
    zero_coord( ave_int );
    zero_coord( ave_ext );

    Coord center_coord; // from Util.hpp
    double center[3];
    for ( size_t h = 0; h < hexes.size(); ++h )
    {
      // get center coordinate of hex
      EntityHandle hex = hexes[h];
      Util::centroid( mbImpl, hex, center_coord );
      center[0] = center_coord.x;
      center[1] = center_coord.y;
      center[2] = center_coord.z;

      if ( get_shrink_membership(hex) == SlabData::INTERNAL )
      {
        ++num_int;
        add_coord( ave_int, center );
      }
      else
      {
        ++num_ext;
        add_coord( ave_ext, center );
      }
    }
    double new_pos[3], node_pos[3];
    zero_coord( new_pos );
    mbImpl->get_coords( &node, 1, node_pos );
    const double frac = 0.2; // fraction towards average hex center, probably 1/2 an edge length
    if ( (num_int <= num_ext) && num_int )
    {
      axpby_coord( new_pos,   1. - frac, node_pos,   frac / num_int, ave_int );
    }
    else if ( num_ext )
    {
      axpby_coord( new_pos,   1. - frac, node_pos,  -frac / num_ext, ave_ext );
    }
    else
      // no hexes, this shouldn't happen
      return;
    mbImpl->set_coords( &node, 1, new_pos );
  }



  void RefineSlabs::Slab::reset()
  {
    edge_set.clear();
    star_nodes.clear();
  }
  bool RefineSlabs::Slab::any_in_slab( const SlabEdges &slab_edges )
  {      
    for (size_t i=0; i<slab_edges.size(); ++i )
    {
      if ( get_in_slab(slab_edges[i]) )
      {
        return true;
      }
    }
    return false;
  }
  bool RefineSlabs::Slab::get_in_slab( const SlabEdge &slab_edge )
  {
    NodePair node_pair(slab_edge.head_node, slab_edge.tail_node);
    EdgeSetIterator i = edge_set.find(node_pair);
    return (i != edge_set.end());
  }
  void RefineSlabs::Slab::set_in_slab( const SlabEdge &slab_edge, bool new_value )
  {
    NodePair node_pair(slab_edge.head_node, slab_edge.tail_node);
    // add?
    if (new_value)
    {
      edge_set.insert( node_pair );
    }
    // erase?
    else
    {
      edge_set.erase( node_pair );
    }
  }

  bool RefineSlabs::Slab::get_is_star( const EntityHandle node )
  {
    EntitySet::iterator i = star_nodes.find(node);
    return (i != star_nodes.end());
  }
  void RefineSlabs::Slab::set_is_star( EntityHandle node, bool new_value )
  {
    // add?
    if (new_value)
    {
      star_nodes.insert( node );
    }
    // erase?
    else
    {
      star_nodes.erase( node );
    }
  }

  void RefineSlabs::register_boundary_nodes( Entities &surface_nodes, Entities &curve_nodes, Entities &vertex_nodes )
  {
    for ( size_t i = 0; i < surface_nodes.size(); ++i )
      geometry_map[ surface_nodes[i] ] = 2;
    for ( size_t i = 0; i < curve_nodes.size(); ++i )
      geometry_map[ curve_nodes[i] ] = 1;
    for ( size_t i = 0; i < vertex_nodes.size(); ++i )
      geometry_map[ vertex_nodes[i] ] = 0;
  }


  void RefineSlabs::register_entity_handles(EntityHandle *hexes, int num_hexes, EntityHandle *vertices, int num_vertices)
  {
    registered_hexes = hexes;
    registered_vertices = vertices;
    registered_num_hexes = num_hexes;
    registered_num_vertices = num_vertices;
  }
  bool RefineSlabs::is_registered_hex( EntityHandle hex )
  {
    // return std::find( registered_hexes.begin(), registered_hexes.end(), hex );
    for (int i = 0; i < registered_num_hexes; ++i)
      if ( registered_hexes[i] == hex )
        return true;
    return false;
  }
  bool RefineSlabs::is_registered_vertex( EntityHandle vertex )
  {
    // return std::find( registered_vertices.begin(), registered_vertices.end(), vertex );
    for (int i = 0; i < registered_num_vertices; ++i)
      if ( registered_vertices[i] == vertex )
        return true;
    return false;
  }

  bool RefineSlabs::node_ok( bool iscoarse, EntityHandle node )
  {
    if ( iscoarse && !is_registered_vertex( node ) )
    {
      std::cout << "bad coarse node" << std::endl;
      assert(0);
      return false;
    }
    else if ( !iscoarse && created_fine_nodes.find(node) == created_fine_nodes.end() )
    {
      std::cout << "bad fine node" << std::endl;
      assert(0);
      return false;
    }
    return true;
  }
  bool RefineSlabs::hex_ok( EntityHandle hex )
  {
    if ( created_fine_hexes.find(hex) != created_fine_hexes.end() || is_registered_hex( hex ) )
    {
      return true;
    }
    
    std::cout << "bad hex" << std::endl;
    assert(0);
    return false;
  }

  void RefineSlabs::print_slab( const Slab & slab )
  {
    Entities hexes;
    get_pillow_hexes( slab, hexes );
    print_slab( hexes );
    print_slab( slab.star_nodes );
    std::cout << std::endl;
  }
  void RefineSlabs::print_slab( const EntitySet &star_nodes )
  {
    using namespace std;
    cout << "Slab with " << star_nodes.size() << " star_nodes: ";
    for ( EntitySet::const_iterator i = star_nodes.begin(); i != star_nodes.end(); ++i )
    {
      cout << mbImpl->id_from_handle( *i ) << " ";
    }
    cout << endl;
  }
  void RefineSlabs::print_slab( const Entities & slab_hexes )
  {
    using namespace std;
    cout << "Slab with " << slab_hexes.size() << " hexes: ";
    for ( size_t i = 0; i < slab_hexes.size(); ++i )
    {
      cout << mbImpl->id_from_handle(slab_hexes[i]) << " ";
      // cout << slab_hexes[i] << " ";
    }
    cout << endl;
  }

  ErrorCode RefineSlabs::write_file_slab( Entities &slab_hexes, std::string fname, int fname_version )
  {
    // see sset.cpp
    // 
    // see Interface.hpp
    // virtual ErrorCode write_file( const char* file_name,
    //                               const char* file_type = 0,
    //                               const char* options = 0,
    //                               const EntityHandle* output_sets = 0,
    //                               int num_output_sets = 0,
    //                               const Tag* tag_list = 0,
    //                               int num_tags = 0 ) = 0;

    if ( slab_hexes.empty() )
      return MB_SUCCESS;

    using namespace std;

    const bool write_vtk(false);
    const bool write_exo(true);

    // Get MOAB instance
    Interface* mb = mbImpl;

    // file name
    stringstream ss;
    if ( file_sequence < 10 )
      ss << 0;
    ss << file_sequence++ << "_";
    stringstream ss2;
      ss2 << "_" << fname_version;
    fname = ss.str() + fname + ss2.str();

    // entity handle to a set of mesh elements?
    EntityHandle newSet;
    ErrorCode rval = mb->create_meshset(MESHSET_SET, newSet);MB_CHK_ERR(rval);

    Range someElems;
    std::copy(slab_hexes.rbegin(), slab_hexes.rend(), range_inserter( someElems ) );
    rval = mb->add_entities(newSet, someElems);MB_CHK_ERR(rval);

    // write a vtk file
    if (write_vtk)
    {
      string fnamevtk(fname);
      fnamevtk += ".vtk";
      mb->write_file(fnamevtk.c_str(), 0, 0, &newSet, 1);
    }

    // write an exodus file
    if (write_exo)
    {
      string fnameexo(fname);
      fnameexo += ".exo";

      Tag mtag;
      rval=mb->tag_get_handle("MATERIAL_SET", mtag);MB_CHK_ERR(rval);

      int val=100;
      rval = mb->tag_set_data(mtag, &newSet, 1, &val);MB_CHK_ERR(rval);      

      rval = mb->write_file(fnameexo.c_str(), 0, 0, &newSet, 1);MB_CHK_ERR(rval);
    }
    return MB_SUCCESS;
  }

}//namesapce moab

