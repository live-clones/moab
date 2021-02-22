/*
 * =====================================================================================
 *
 *       Filename:  Remapper.hpp
 *
 *    Description:  Interface to the a general remapping capability on arbitrary topology
 *                  that performs both mesh intersection between a source and target grid,
 *                  with arbitrary decompositions. The intersections can then be used to
 *                  either evaluate interpolation weights or to perform high-order
 *                  conservative remapping of solutions defined on the source grid.
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *
 * =====================================================================================
 */

#ifndef MB_REMAPPER_HPP
#define MB_REMAPPER_HPP

#include <string>

#include "moab/Interface.hpp"
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif

// Tempest includes
#ifdef MOAB_HAVE_TEMPESTREMAP
#include "netcdfcpp.h"
#include "TempestRemapAPI.h"
#else
#error "This tool depends on TempestRemap library. Reconfigure using --with-tempestremap"
#endif

namespace moab
{

class Remapper
{
  public:
#ifdef MOAB_HAVE_MPI
    Remapper( moab::Interface* mbInt, moab::ParallelComm* pcomm = NULL ) : m_interface( mbInt ), m_pcomm( pcomm )
#else
    Remapper( moab::Interface* mbInt ) : m_interface( mbInt )
#endif
    {
    }

    virtual ~Remapper()
    {
#ifdef MOAB_HAVE_MPI
        m_pcomm = NULL;
#endif
        m_interface = NULL;
    }

    enum IntersectionContext
    {
        DEFAULT      = -1,
        SourceMesh   = 0,
        TargetMesh   = 1,
        OverlapMesh  = 2,
        CoveringMesh = 3
    };

    moab::Interface* get_interface()
    {
        return m_interface;
    }

#ifdef MOAB_HAVE_MPI
    moab::ParallelComm* get_parallel_communicator()
    {
        return m_pcomm;
    }
#endif

    ErrorCode LoadNativeMesh( std::string filename, moab::EntityHandle& meshset, const char* readopts = 0 )
    {
#ifdef MOAB_HAVE_MPI
        size_t lastindex      = filename.find_last_of( "." );
        std::string extension = filename.substr( lastindex + 1, filename.size() );
        std::string opts      = "";
        if( m_pcomm->size() > 1 )
        {
            if( extension != "h5m" )
                opts = std::string( "PARALLEL=BCAST_DELETE;PARTITION=TRIVIAL;PARALLEL_RESOLVE_SHARED_ENTS" );
            else
                opts = std::string( "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_"
                                    "RESOLVE_SHARED_ENTS" );
        }

        if( readopts )
        {
            if( opts.size() )
                opts = opts + ";" + std::string( readopts );
            else
                opts = std::string( readopts );
        }

        if( !m_pcomm->rank() ) std::cout << "Reading file (" << filename << ") with options = [" << opts << "]\n";
#else
        const std::string opts = std::string( ( readopts ? readopts : "" ) );
        std::cout << "Reading file (" << filename << ") with options = [" << opts << "]\n";
#endif
        return m_interface->load_file( filename.c_str(), &meshset, opts.c_str() );
    }

  protected:
    // member data
    Interface* m_interface;

#ifdef MOAB_HAVE_MPI
    ParallelComm* m_pcomm;
#endif
};

}  // namespace moab

#endif /* MB_REMAPPER_HPP */
