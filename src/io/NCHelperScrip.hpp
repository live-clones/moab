/*
 * NCHelperScrip.hpp
 *  Purpose       : Climate NC file helper for Scrip grid
 */

#ifndef SRC_IO_NCHELPERSCRIP_HPP_
#define SRC_IO_NCHELPERSCRIP_HPP_

#include "NCHelper.hpp"

namespace moab
{

class NCHelperScrip : public ScdNCHelper
{
  public:
    NCHelperScrip( ReadNC* readNC, int fileId, const FileOptions& opts, EntityHandle fileSet )
        : ScdNCHelper( readNC, fileId, opts, fileSet ), grid_corners( 0 ), grid_size( 0 ), nLocalCells( 0 ),
          degrees( true )
    {
    }
    static bool can_read_file( ReadNC* readNC, int fileId );

    ErrorCode create_mesh( Range& faces );

#ifdef MOAB_HAVE_MPI
    //! Redistribute local cells after trivial partition (e.g. Zoltan partition, if applicable)
    ErrorCode redistribute_local_cells( int start_cell_index );
#endif

  private:
    virtual ErrorCode init_mesh_vals();
    virtual std::string get_mesh_type_name()
    {
        return "SCRIP";
    }

    int grid_corners;  // number of vertices per cell
    int grid_size;
    int nLocalCells;      // in parallel, number of local cells, initially, and after repartition
    Range localGidCells;  // will store the ids after repartitioning;
    bool degrees;         // if false, it means it is radians
};

} /* namespace moab */

#endif /* SRC_IO_NCHELPERSCRIP_HPP_ */
