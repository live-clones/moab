/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 *
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */

// If Microsoft compiler, then WIN32
#ifndef WIN32
#define WIN32
#endif

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "MBTagConventions.hpp"
#include "moab/ReaderWriterSet.hpp"
#include "moab/ReorderTool.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <set>
#include <list>
#include <cstdlib>
#include <algorithm>
#ifndef WIN32
#include <sys/times.h>
#include <limits.h>
#include <unistd.h>
#endif
#include <ctime>
#ifdef MOAB_HAVE_MPI
#include "moab/ParallelComm.hpp"
#endif

#ifdef MOAB_HAVE_TEMPESTREMAP
#include "moab/Remapping/TempestRemapper.hpp"
#endif

#include <cstdio>

/* Exit values */
#define USAGE_ERROR   1
#define READ_ERROR    2
#define WRITE_ERROR   3
#define OTHER_ERROR   4
#define ENT_NOT_FOUND 5

using namespace moab;

static void print_usage( const char* name, std::ostream& stream )
{
    stream << "Usage: " << name
           << " [-a <sat_file>|-A] [-t] [subset options] [-f format] <input_file> [<input_file2> "
              "...] <output_file>"
           << std::endl
           << "\t-f <format>    - Specify output file format" << std::endl
           << "\t-a <acis_file> - ACIS SAT file dumped by .cub reader (same as \"-o "
              "SAT_FILE=acis_file\""
           << std::endl
           << "\t-A             - .cub file reader should not dump a SAT file (depricated default)" << std::endl
           << "\t-o option      - Specify write option." << std::endl
           << "\t-O option      - Specify read option." << std::endl
           << "\t-t             - Time read and write of files." << std::endl
           << "\t-g             - Enable verbose/debug output." << std::endl
           << "\t-h             - Print this help text and exit." << std::endl
           << "\t-l             - List available file formats and exit." << std::endl
           << "\t-I <dim>       - Generate internal entities of specified dimension." << std::endl
#ifdef MOAB_HAVE_MPI
           << "\t-P             - Append processor ID to output file name" << std::endl
           << "\t-p             - Replace '%' with processor ID in input and output file name" << std::endl
           << "\t-M[0|1|2]      - Read/write in parallel, optionally also doing "
              "resolve_shared_ents (1) and exchange_ghosts (2)"
           << std::endl
           << "\t-z <file>      - Read metis partition information corresponding to an MPAS grid "
              "file and create h5m partition file"
           << std::endl
#endif

#ifdef MOAB_HAVE_TEMPESTREMAP
           << "\t-B             - Use TempestRemap exodus file reader and convert to MOAB format" << std::endl
           << "\t-b             - Convert MOAB mesh to TempestRemap exodus file writer" << std::endl
           << "\t-i             - Name of the global DoF tag to use with mbtempest" << std::endl
           << "\t-r             - Order of field DoF (discretization) data; FV=1,SE=[1,N]" << std::endl
#endif
           << "\t--             - treat all subsequent options as file names" << std::endl
           << "\t                 (allows file names beginning with '-')" << std::endl
           << "  subset options: " << std::endl
           << "\tEach of the following options should be followed by " << std::endl
           << "\ta list of ids.  IDs may be separated with commas.  " << std::endl
           << "\tRanges of IDs may be specified with a '-' between " << std::endl
           << "\ttwo values.  The list may not contain spaces." << std::endl
           << "\t-v  - Volume" << std::endl
           << "\t-s  - Surface" << std::endl
           << "\t-c  - Curve" << std::endl
           << "\t-V  - Vertex" << std::endl
           << "\t-m  - Material set (block)" << std::endl
           << "\t-d  - Dirichlet set (nodeset)" << std::endl
           << "\t-n  - Neumann set (sideset)" << std::endl
           << "\t-D  - Parallel partitioning set (PARALLEL_PARTITION)" << std::endl
           << "\tThe presence of one or more of the following flags limits " << std::endl
           << "\tthe exported mesh to only elements of the corresponding " << std::endl
           << "\tdimension.  Vertices are always exported." << std::endl
           << "\t-1  - Edges " << std::endl
           << "\t-2  - Tri, Quad, Polygon " << std::endl
           << "\t-3  - Tet, Hex, Prism, etc. " << std::endl;
}

static void print_help( const char* name )
{
    std::cout << " This program can be used to convert between mesh file\n"
                 " formats, extract a subset of a mesh file to a separate\n"
                 " file, or both.  The type of file to write is determined\n"
                 " from the file extension (e.g. \".vtk\") portion of the\n"
                 " output file name.\n"
                 " \n"
                 " While MOAB strives to export and import all data from\n"
                 " each supported file format, most file formats do\n"
                 " not support MOAB's entire data model.  Thus MOAB cannot\n"
                 " guarantee lossless conversion for any file formats\n"
                 " other than the native HDF5 representation.\n"
                 "\n";

    print_usage( name, std::cout );
    exit( 0 );
}

static void usage_error( const char* name )
{
    print_usage( name, std::cerr );
#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    exit( USAGE_ERROR );
}

static void list_formats( Interface* );
static bool parse_id_list( const char* string, std::set< int >& );
static void print_id_list( const char*, std::ostream& stream, const std::set< int >& list );
static void reset_times();
static void write_times( std::ostream& stream );
static void remove_entities_from_sets( Interface* gMB, Range& dead_entities, Range& empty_sets );
static void remove_from_vector( std::vector< EntityHandle >& vect, const Range& ents_to_remove );
static bool make_opts_string( std::vector< std::string > options, std::string& result );
static std::string percent_subst( const std::string& s, int val );

static int process_partition_file( Interface* gMB, std::string& metis_partition_file );

int main( int argc, char* argv[] )
{
    int proc_id = 0;
#ifdef MOAB_HAVE_MPI
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
#endif

#ifdef MOAB_HAVE_TEMPESTREMAP
    bool tempestin = false, tempestout = false;
#endif

    Core core;
    Interface* gMB = &core;
    ErrorCode result;
    Range range;

#if( defined( MOAB_HAVE_MPI ) && defined( MOAB_HAVE_TEMPESTREMAP ) )
    moab::ParallelComm* pcomm = new moab::ParallelComm( gMB, MPI_COMM_WORLD, 0 );
#endif

    bool append_rank        = false;
    bool percent_rank_subst = false;
    int i, dim;
    std::list< std::string >::iterator j;
    bool dims[4]       = { false, false, false, false };
    const char* format = NULL;    // output file format
    std::list< std::string > in;  // input file name list
    std::string out;              // output file name
    bool verbose = false;
    std::set< int > geom[4], mesh[4];      // user-specified IDs
    std::vector< EntityHandle > set_list;  // list of user-specified sets to write
    std::vector< std::string > write_opts, read_opts;
    std::string metis_partition_file;
#ifdef MOAB_HAVE_TEMPESTREMAP
    std::string globalid_tag_name;
    int spectral_order = 1;
#endif

    const char* const mesh_tag_names[] = { DIRICHLET_SET_TAG_NAME, NEUMANN_SET_TAG_NAME, MATERIAL_SET_TAG_NAME,
                                           PARALLEL_PARTITION_TAG_NAME };
    const char* const geom_names[]     = { "VERTEX", "CURVE", "SURFACE", "VOLUME" };

    // scan arguments
    bool do_flag     = true;
    bool print_times = false;
    bool generate[]  = { false, false, false };
    bool pval;
    bool parallel = false, resolve_shared = false, exchange_ghosts = false;
    for( i = 1; i < argc; i++ )
    {
        if( !argv[i][0] ) usage_error( argv[0] );

        if( do_flag && argv[i][0] == '-' )
        {
            if( !argv[i][1] || ( argv[i][1] != 'M' && argv[i][2] ) ) usage_error( argv[0] );

            switch( argv[i][1] )
            {
                    // do flag arguments:
                case '-':
                    do_flag = false;
                    break;
                case 'g':
                    verbose = true;
                    break;
                case 't':
                    print_times = true;
                    break;
                case 'A':
                    break;
                case 'h':
                case 'H':
                    print_help( argv[0] );
                    break;
                case 'l':
                    list_formats( gMB );
                    break;
#ifdef MOAB_HAVE_MPI
                case 'P':
                    append_rank = true;
                    break;
                case 'p':
                    percent_rank_subst = true;
                    break;
                case 'M':
                    parallel = true;
                    if( argv[i][2] == '1' || argv[i][2] == '2' ) resolve_shared = true;
                    if( argv[i][2] == '2' ) exchange_ghosts = true;
#endif
#ifdef MOAB_HAVE_TEMPESTREMAP
                case 'B':
                    tempestin = true;
                    break;
                case 'b':
                    tempestout = true;
                    break;
#endif
                case '1':
                case '2':
                case '3':
                    dims[argv[i][1] - '0'] = true;
                    break;
                    // do options that require additional args:
                default:
                    ++i;
                    if( i == argc || argv[i][0] == '-' )
                    {
                        std::cerr << "Expected argument following " << argv[i - 1] << std::endl;
                        usage_error( argv[0] );
                    }
                    if( argv[i - 1][1] == 'I' )
                    {
                        dim = atoi( argv[i] );
                        if( dim < 1 || dim > 2 )
                        {
                            std::cerr << "Invalid dimension value following -I" << std::endl;
                            usage_error( argv[0] );
                        }
                        generate[dim] = true;
                        continue;
                    }
                    pval = false;
                    switch( argv[i - 1][1] )
                    {
                        case 'a':
                            read_opts.push_back( std::string( "SAT_FILE=" ) + argv[i] );
                            pval = true;
                            break;
                        case 'f':
                            format = argv[i];
                            pval   = true;
                            break;
                        case 'o':
                            write_opts.push_back( argv[i] );
                            pval = true;
                            break;
                        case 'O':
                            read_opts.push_back( argv[i] );
                            pval = true;
                            break;
#ifdef MOAB_HAVE_TEMPESTREMAP
                        case 'i':
                            globalid_tag_name = std::string( argv[i] );
                            pval              = true;
                            break;
                        case 'r':
                            spectral_order = atoi( argv[i] );
                            pval           = true;
                            break;
#endif
                        case 'v':
                            pval = parse_id_list( argv[i], geom[3] );
                            break;
                        case 's':
                            pval = parse_id_list( argv[i], geom[2] );
                            break;
                        case 'c':
                            pval = parse_id_list( argv[i], geom[1] );
                            break;
                        case 'V':
                            pval = parse_id_list( argv[i], geom[0] );
                            break;
                        case 'D':
                            pval = parse_id_list( argv[i], mesh[3] );
                            break;
                        case 'm':
                            pval = parse_id_list( argv[i], mesh[2] );
                            break;
                        case 'n':
                            pval = parse_id_list( argv[i], mesh[1] );
                            break;
                        case 'd':
                            pval = parse_id_list( argv[i], mesh[0] );
                            break;
                        case 'z':
                            metis_partition_file = argv[i];
                            pval                 = true;
                            break;
                        default:
                            std::cerr << "Invalid option: " << argv[i] << std::endl;
                    }

                    if( !pval )
                    {
                        std::cerr << "Invalid flag or flag value: " << argv[i - 1] << " " << argv[i] << std::endl;
                        usage_error( argv[0] );
                    }
            }
        }
        // do file names
        else
        {
            in.push_back( argv[i] );
        }
    }
    if( in.size() < 2 )
    {
        std::cerr << "No output file name specified." << std::endl;
        usage_error( argv[0] );
    }
    // output file name is the last one specified
    out = in.back();
    in.pop_back();

    if( append_rank )
    {
        std::ostringstream mod;
        mod << out << "." << proc_id;
        out = mod.str();
    }

    if( percent_rank_subst )
    {
        for( j = in.begin(); j != in.end(); ++j )
            *j = percent_subst( *j, proc_id );
        out = percent_subst( out, proc_id );
    }

    // construct options string from individual options
    std::string read_options, write_options;
    if( parallel )
    {
        read_opts.push_back( "PARALLEL=READ_PART" );
        read_opts.push_back( "PARTITION=PARALLEL_PARTITION" );
        if( !append_rank && !percent_rank_subst ) write_opts.push_back( "PARALLEL=WRITE_PART" );
    }
    if( resolve_shared ) read_opts.push_back( "PARALLEL_RESOLVE_SHARED_ENTS" );
    if( exchange_ghosts ) read_opts.push_back( "PARALLEL_GHOSTS=3.0.1" );

    if( !make_opts_string( read_opts, read_options ) || !make_opts_string( write_opts, write_options ) )
    {
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        return USAGE_ERROR;
    }

    if( !metis_partition_file.empty() )
    {
        if( ( in.size() != 1 ) || ( proc_id != 0 ) )
        {
            std::cerr << " mpas partition allows only one input file, in serial conversion\n";
#ifdef MOAB_HAVE_MPI
            MPI_Finalize();
#endif
            return USAGE_ERROR;
        }
    }

    Tag id_tag = gMB->globalId_tag();

    // Read the input file.
#ifdef MOAB_HAVE_TEMPESTREMAP
    if( tempestin && in.size() > 1 )
    {
        std::cerr << " we can read only one tempest files at a time\n";
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        return USAGE_ERROR;
    }

#ifdef MOAB_HAVE_MPI
    TempestRemapper* remapper = new moab::TempestRemapper( gMB, pcomm );
#else
    TempestRemapper* remapper = new moab::TempestRemapper( gMB );
#endif

    bool use_overlap_context = false;
    Tag srcParentTag, tgtParentTag;

#endif
    for( j = in.begin(); j != in.end(); ++j )
    {
        std::string inFileName = *j;

        reset_times();

#ifdef MOAB_HAVE_TEMPESTREMAP

        remapper->meshValidate = false;
        // remapper->constructEdgeMap = true;
        remapper->initialize();

        if( tempestin )
        {
            // convert
            result = remapper->LoadMesh( moab::Remapper::SourceMesh, inFileName, moab::TempestRemapper::DEFAULT );MB_CHK_ERR( result );

            Mesh* tempestMesh = remapper->GetMesh( moab::Remapper::SourceMesh );
            tempestMesh->RemoveZeroEdges();
            tempestMesh->RemoveCoincidentNodes();

            // Load the meshes and validate
            result = remapper->ConvertTempestMesh( moab::Remapper::SourceMesh );

            // Check if we are converting a RLL grid
            NcFile ncInput( inFileName.c_str(), NcFile::ReadOnly );
            bool isRectilinearGrid = false;

            NcError error_temp( NcError::silent_nonfatal );
            // get the attribute
            NcAtt* attRectilinear = ncInput.get_att( "rectilinear" );

            // If rectilinear attribute present, mark it
            std::vector< int > vecDimSizes( 3, 0 );
            Tag rectilinearTag;
            // Tag data contains: guessed mesh type,     mesh size1,     mesh size 2
            //          Example:  CS(0)/ICO(1)/ICOD(2),  num_elements,   num_nodes
            //                 :       RLL(3),           num_lat,        num_lon
            result = gMB->tag_get_handle( "ClimateMetadata", 3, MB_TYPE_INTEGER, rectilinearTag,
                                          MB_TAG_SPARSE | MB_TAG_CREAT, vecDimSizes.data() );MB_CHK_SET_ERR( result, "can't create rectilinear sizes tag" );

            if( attRectilinear != nullptr )
            {
                isRectilinearGrid = true;

                // Obtain rectilinear attributes (dimension sizes)
                NcAtt* attRectilinearDim0Size = ncInput.get_att( "rectilinear_dim0_size" );
                NcAtt* attRectilinearDim1Size = ncInput.get_att( "rectilinear_dim1_size" );

                if( attRectilinearDim0Size == nullptr )
                {
                    _EXCEPTIONT( "Missing attribute \"rectilinear_dim0_size\"" );
                }
                if( attRectilinearDim1Size == nullptr )
                {
                    _EXCEPTIONT( "Missing attribute \"rectilinear_dim1_size\"" );
                }

                int nDim0Size = attRectilinearDim0Size->as_int( 0 );
                int nDim1Size = attRectilinearDim1Size->as_int( 0 );

                // Obtain rectilinear attributes (dimension names)
                NcAtt* attRectilinearDim0Name = ncInput.get_att( "rectilinear_dim0_name" );
                NcAtt* attRectilinearDim1Name = ncInput.get_att( "rectilinear_dim1_name" );

                if( attRectilinearDim0Name == nullptr )
                {
                    _EXCEPTIONT( "Missing attribute \"rectilinear_dim0_name\"" );
                }
                if( attRectilinearDim1Name == nullptr )
                {
                    _EXCEPTIONT( "Missing attribute \"rectilinear_dim1_name\"" );
                }

                std::string strDim0Name = attRectilinearDim0Name->as_string( 0 );
                std::string strDim1Name = attRectilinearDim1Name->as_string( 0 );

                std::map< std::string, int > vecDimNameSizes;
                // Push rectilinear attributes into array
                vecDimNameSizes[strDim0Name] = nDim0Size;
                vecDimNameSizes[strDim1Name] = nDim1Size;
                vecDimSizes[0]               = static_cast< int >( moab::TempestRemapper::RLL );
                vecDimSizes[1]               = vecDimNameSizes["lat"];
                vecDimSizes[2]               = vecDimNameSizes["lon"];

                printf( "Rectilinear RLL mesh size: (lat) %d X (lon) %d\n", vecDimSizes[1], vecDimSizes[2] );

                moab::EntityHandle mSet = 0;
                // mSet   = remapper->GetMeshSet( moab::Remapper::SourceMesh );
                result = gMB->tag_set_data( rectilinearTag, &mSet, 1, vecDimSizes.data() );MB_CHK_ERR( result );
            }
            else
            {
                const Range& elems = remapper->GetMeshEntities( moab::Remapper::SourceMesh );
                bool isQuads       = elems.all_of_type( moab::MBQUAD );
                bool isTris        = elems.all_of_type( moab::MBTRI );
                // vecDimSizes[0] = static_cast< int >( remapper->GetMeshType( moab::Remapper::SourceMesh );
                vecDimSizes[0] = ( isQuads ? static_cast< int >( moab::TempestRemapper::CS )
                                           : ( isTris ? static_cast< int >( moab::TempestRemapper::ICO )
                                                      : static_cast< int >( moab::TempestRemapper::ICOD ) ) );
                vecDimSizes[1] = elems.size();
                vecDimSizes[2] = remapper->GetMeshVertices( moab::Remapper::SourceMesh ).size();

                switch( vecDimSizes[0] )
                {
                    case 0:
                        printf( "Cubed-Sphere mesh: %d (elems), %d (nodes)\n", vecDimSizes[1], vecDimSizes[2] );
                        break;
                    case 2:
                        printf( "Icosahedral (triangular) mesh: %d (elems), %d (nodes)\n", vecDimSizes[1],
                                vecDimSizes[2] );
                        break;
                    case 3:
                    default:
                        printf( "Polygonal mesh: %d (elems), %d (nodes)\n", vecDimSizes[1], vecDimSizes[2] );
                        break;
                }

                moab::EntityHandle mSet = 0;
                // mSet   = remapper->GetMeshSet( moab::Remapper::SourceMesh );
                result = gMB->tag_set_data( rectilinearTag, &mSet, 1, vecDimSizes.data() );MB_CHK_ERR( result );
            }

            ncInput.close();

            const size_t nOverlapFaces = tempestMesh->faces.size();
            if( tempestMesh->vecSourceFaceIx.size() == nOverlapFaces &&
                tempestMesh->vecSourceFaceIx.size() == nOverlapFaces )
            {
                int defaultInt      = -1;
                use_overlap_context = true;
                // Check if our MOAB mesh has RED and BLUE tags; this would indicate we are
                // converting an overlap grid
                result = gMB->tag_get_handle( "TargetParent", 1, MB_TYPE_INTEGER, tgtParentTag,
                                              MB_TAG_DENSE | MB_TAG_CREAT, &defaultInt );MB_CHK_SET_ERR( result, "can't create target parent tag" );

                result = gMB->tag_get_handle( "SourceParent", 1, MB_TYPE_INTEGER, srcParentTag,
                                              MB_TAG_DENSE | MB_TAG_CREAT, &defaultInt );MB_CHK_SET_ERR( result, "can't create source parent tag" );

                const Range& faces = remapper->GetMeshEntities( moab::Remapper::SourceMesh );

                std::vector< int > gids( faces.size() ), srcpar( faces.size() ), tgtpar( faces.size() );
                result = gMB->tag_get_data( id_tag, faces, &gids[0] );MB_CHK_ERR( result );

                for( unsigned ii = 0; ii < faces.size(); ++ii )
                {
                    srcpar[ii] = tempestMesh->vecSourceFaceIx[gids[ii] - 1];
                    tgtpar[ii] = tempestMesh->vecTargetFaceIx[gids[ii] - 1];
                }

                result = gMB->tag_set_data( srcParentTag, faces, &srcpar[0] );MB_CHK_ERR( result );
                result = gMB->tag_set_data( tgtParentTag, faces, &tgtpar[0] );MB_CHK_ERR( result );

                srcpar.clear();
                tgtpar.clear();
                gids.clear();
            }
        }
        else if( tempestout )
        {
            moab::EntityHandle& srcmesh = remapper->GetMeshSet( moab::Remapper::SourceMesh );
            moab::EntityHandle& ovmesh  = remapper->GetMeshSet( moab::Remapper::OverlapMesh );

            // load the mesh in MOAB format
            std::vector< int > metadata;
            result = remapper->LoadNativeMesh( *j, srcmesh, metadata );MB_CHK_ERR( result );

            // Check if our MOAB mesh has RED and BLUE tags; this would indicate we are converting
            // an overlap grid
            ErrorCode rval1 = gMB->tag_get_handle( "SourceParent", srcParentTag );
            ErrorCode rval2 = gMB->tag_get_handle( "TargetParent", tgtParentTag );
            if( rval1 == MB_SUCCESS && rval2 == MB_SUCCESS )
            {
                use_overlap_context = true;
                ovmesh              = srcmesh;

                Tag countTag;
                result = gMB->tag_get_handle( "Counting", countTag );
                // std::vector<int> count_ids()

                // Load the meshes and validate
                Tag order;
                ReorderTool reorder_tool( &core );
                result = reorder_tool.handle_order_from_int_tag( srcParentTag, -1, order );MB_CHK_ERR( result );
                result = reorder_tool.reorder_entities( order );MB_CHK_ERR( result );
                result = gMB->tag_delete( order );MB_CHK_ERR( result );
                result = remapper->ConvertMeshToTempest( moab::Remapper::OverlapMesh );MB_CHK_ERR( result );
            }
            else
            {
                if( metadata[0] == static_cast< int >( moab::TempestRemapper::RLL ) )
                {
                    assert( metadata.size() );
                    std::cout << "Converting a RLL mesh with rectilinear dimension: " << metadata[0] << " X "
                              << metadata[1] << std::endl;
                }

                // Convert the mesh and validate
                result = remapper->ConvertMeshToTempest( moab::Remapper::SourceMesh );MB_CHK_ERR( result );
            }
        }
        else
            result = gMB->load_file( j->c_str(), 0, read_options.c_str() );
#else
        result = gMB->load_file( j->c_str(), 0, read_options.c_str() );
#endif
        if( MB_SUCCESS != result )
        {
            std::cerr << "Failed to load \"" << *j << "\"." << std::endl;
            std::cerr << "Error code: " << gMB->get_error_string( result ) << " (" << result << ")" << std::endl;
            std::string message;
            if( MB_SUCCESS == gMB->get_last_error( message ) && !message.empty() )
                std::cerr << "Error message: " << message << std::endl;
#ifdef MOAB_HAVE_MPI
            MPI_Finalize();
#endif
            return READ_ERROR;
        }
        if( !proc_id ) std::cerr << "Read \"" << *j << "\"" << std::endl;
        if( print_times && !proc_id ) write_times( std::cout );
    }

    // Determine if the user has specified any geometry sets to write
    bool have_geom = false;
    for( dim = 0; dim <= 3; ++dim )
    {
        if( !geom[dim].empty() ) have_geom = true;
        if( verbose ) print_id_list( geom_names[dim], std::cout, geom[dim] );
    }

    // True if the user has specified any sets to write
    bool have_sets = have_geom;

    // Get geometry tags
    Tag dim_tag;
    if( have_geom )
    {
        if( id_tag == 0 )
        {
            std::cerr << "No ID tag defined." << std::endl;
            have_geom = false;
        }
        result = gMB->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, dim_tag );
        if( MB_SUCCESS != result )
        {
            std::cerr << "No geometry tag defined." << std::endl;
            have_geom = false;
        }
    }

    // Get geometry sets
    if( have_geom )
    {
        int id_val;
        Tag tags[]         = { id_tag, dim_tag };
        const void* vals[] = { &id_val, &dim };
        for( dim = 0; dim <= 3; ++dim )
        {
            int init_count = set_list.size();
            for( std::set< int >::iterator iter = geom[dim].begin(); iter != geom[dim].end(); ++iter )
            {
                id_val = *iter;
                range.clear();
                result = gMB->get_entities_by_type_and_tag( 0, MBENTITYSET, tags, vals, 2, range );
                if( MB_SUCCESS != result || range.empty() )
                {
                    range.clear();
                    std::cerr << geom_names[dim] << " " << id_val << " not found.\n";
                }
                std::copy( range.begin(), range.end(), std::back_inserter( set_list ) );
            }

            if( verbose )
                std::cout << "Found " << ( set_list.size() - init_count ) << ' ' << geom_names[dim] << " sets"
                          << std::endl;
        }
    }

    // Get mesh groupings
    for( i = 0; i < 4; ++i )
    {
        if( verbose ) print_id_list( mesh_tag_names[i], std::cout, mesh[i] );

        if( mesh[i].empty() ) continue;
        have_sets = true;

        // Get tag
        Tag tag;
        result = gMB->tag_get_handle( mesh_tag_names[i], 1, MB_TYPE_INTEGER, tag );
        if( MB_SUCCESS != result )
        {
            std::cerr << "Tag not found: " << mesh_tag_names[i] << std::endl;
            continue;
        }

        // get entity sets
        int init_count = set_list.size();
        for( std::set< int >::iterator iter = mesh[i].begin(); iter != mesh[i].end(); ++iter )
        {
            range.clear();
            const void* vals[] = { &*iter };
            result             = gMB->get_entities_by_type_and_tag( 0, MBENTITYSET, &tag, vals, 1, range );
            if( MB_SUCCESS != result || range.empty() )
            {
                range.clear();
                std::cerr << mesh_tag_names[i] << " " << *iter << " not found.\n";
            }
            std::copy( range.begin(), range.end(), std::back_inserter( set_list ) );
        }

        if( verbose )
            std::cout << "Found " << ( set_list.size() - init_count ) << ' ' << mesh_tag_names[i] << " sets"
                      << std::endl;
    }

    // Check if output is limited to certain dimensions of elements
    bool bydim = false;
    for( dim = 1; dim < 4; ++dim )
        if( dims[dim] ) bydim = true;

    // Check conflicting input
    if( bydim )
    {
        if( generate[1] && !dims[1] )
        {
            std::cerr << "Warning: Request to generate 1D internal entities but not export them." << std::endl;
            generate[1] = false;
        }
        if( generate[2] && !dims[2] )
        {
            std::cerr << "Warning: Request to generate 2D internal entities but not export them." << std::endl;
            generate[2] = false;
        }
    }

    // Generate any internal entities
    if( generate[1] || generate[2] )
    {
        EntityHandle all_mesh    = 0;
        const EntityHandle* sets = &all_mesh;
        int num_sets             = 1;
        if( have_sets )
        {
            num_sets = set_list.size();
            sets     = &set_list[0];
        }
        for( i = 0; i < num_sets; ++i )
        {
            Range dim3, dim2, adj;
            gMB->get_entities_by_dimension( sets[i], 3, dim3, true );
            if( generate[1] )
            {
                gMB->get_entities_by_dimension( sets[i], 2, dim2, true );
                gMB->get_adjacencies( dim3, 1, true, adj, Interface::UNION );
                gMB->get_adjacencies( dim2, 1, true, adj, Interface::UNION );
            }
            if( generate[2] ) { gMB->get_adjacencies( dim3, 2, true, adj, Interface::UNION ); }
            if( sets[i] ) gMB->add_entities( sets[i], adj );
        }
    }

    // Delete any entities not of the dimensions to be exported
    if( bydim )
    {
        // Get list of dead elements
        Range dead_entities, tmp_range;
        for( dim = 1; dim <= 3; ++dim )
        {
            if( dims[dim] ) continue;
            gMB->get_entities_by_dimension( 0, dim, tmp_range );
            dead_entities.merge( tmp_range );
        }
        // Remove dead entities from all sets, and add all
        // empty sets to list of dead entities.
        Range empty_sets;
        remove_entities_from_sets( gMB, dead_entities, empty_sets );
        while( !empty_sets.empty() )
        {
            if( !set_list.empty() ) remove_from_vector( set_list, empty_sets );
            dead_entities.merge( empty_sets );
            tmp_range.clear();
            remove_entities_from_sets( gMB, empty_sets, tmp_range );
            empty_sets = subtract( tmp_range, dead_entities );
        }
        // Destroy dead entities
        gMB->delete_entities( dead_entities );
    }

    // If user specified sets to write, but none were found, exit.
    if( have_sets && set_list.empty() )
    {
        std::cerr << "Nothing to write." << std::endl;
#ifdef MOAB_HAVE_MPI
        MPI_Finalize();
#endif
        return ENT_NOT_FOUND;
    }

    // interpret the mpas partition file created by gpmetis
    if( !metis_partition_file.empty() )
    {
        int err = process_partition_file( gMB, metis_partition_file );
        if( err )
        {
            std::cerr << "Failed to process partition file \"" << metis_partition_file << "\"." << std::endl;
#ifdef MOAB_HAVE_MPI
            MPI_Finalize();
#endif
            return WRITE_ERROR;
        }
    }
    if( verbose )
    {
        if( have_sets )
            std::cout << "Found " << set_list.size() << " specified sets to write (total)." << std::endl;
        else
            std::cout << "No sets specifed.  Writing entire mesh." << std::endl;
    }

    // Write the output file
    reset_times();
#ifdef MOAB_HAVE_TEMPESTREMAP
    Range faces;
    Mesh* tempestMesh =
        remapper->GetMesh( ( use_overlap_context ? moab::Remapper::OverlapMesh : moab::Remapper::SourceMesh ) );
    moab::EntityHandle& srcmesh =
        remapper->GetMeshSet( ( use_overlap_context ? moab::Remapper::OverlapMesh : moab::Remapper::SourceMesh ) );
    result = gMB->get_entities_by_dimension( srcmesh, 2, faces );MB_CHK_ERR( result );
    int ntot_elements = 0, nelements = faces.size();
#ifdef MOAB_HAVE_MPI
    int ierr = MPI_Allreduce( &nelements, &ntot_elements, 1, MPI_INT, MPI_SUM, pcomm->comm() );
    if( ierr != 0 ) MB_CHK_SET_ERR( MB_FAILURE, "MPI_Allreduce failed to get total source elements" );
#else
    ntot_elements             = nelements;
#endif

    Tag gidTag = gMB->globalId_tag();
    std::vector< int > gids( faces.size() );
    result = gMB->tag_get_data( gidTag, faces, &gids[0] );MB_CHK_ERR( result );

    if( faces.size() > 1 && gids[0] == gids[1] && !use_overlap_context )
    {
#ifdef MOAB_HAVE_MPI
        result = pcomm->assign_global_ids( srcmesh, 2, 1, false );MB_CHK_ERR( result );
#else
        result = remapper->assign_vertex_element_IDs( gidTag, srcmesh, 2, 1 );MB_CHK_ERR( result );
        result = remapper->assign_vertex_element_IDs( gidTag, srcmesh, 0, 1 );MB_CHK_ERR( result );
#endif
    }

    // VSM: If user requested explicitly for some metadata, we need to generate the DoF ID tag
    // and set the appropriate numbering based on specified discretization order
    // Useful only for SE meshes with GLL DoFs
    if( spectral_order > 1 && globalid_tag_name.size() > 1 )
    {
        result = remapper->GenerateMeshMetadata( *tempestMesh, ntot_elements, faces, NULL, globalid_tag_name,
                                                 spectral_order );MB_CHK_ERR( result );
    }

    if( tempestout )
    {
        // Check if our MOAB mesh has RED and BLUE tags; this would indicate we are converting an
        // overlap grid
        if( use_overlap_context && false )
        {
            const int nOverlapFaces = faces.size();
            // Overlap mesh: resize the source and target connection arrays
            tempestMesh->vecSourceFaceIx.resize( nOverlapFaces );  // 0-based indices corresponding to source mesh
            tempestMesh->vecTargetFaceIx.resize( nOverlapFaces );  // 0-based indices corresponding to target mesh
            result = gMB->tag_get_data( srcParentTag, faces, &tempestMesh->vecSourceFaceIx[0] );MB_CHK_ERR( result );
            result = gMB->tag_get_data( tgtParentTag, faces, &tempestMesh->vecTargetFaceIx[0] );MB_CHK_ERR( result );
        }
        // Write out the mesh using TempestRemap
        tempestMesh->Write( out, NcFile::Netcdf4 );
    }
    else
    {
#endif

        if( have_sets )
            result = gMB->write_file( out.c_str(), format, write_options.c_str(), &set_list[0], set_list.size() );
        else
            result = gMB->write_file( out.c_str(), format, write_options.c_str() );
        if( MB_SUCCESS != result )
        {
            std::cerr << "Failed to write \"" << out << "\"." << std::endl;
            std::cerr << "Error code: " << gMB->get_error_string( result ) << " (" << result << ")" << std::endl;
            std::string message;
            if( MB_SUCCESS == gMB->get_last_error( message ) && !message.empty() )
                std::cerr << "Error message: " << message << std::endl;
#ifdef MOAB_HAVE_MPI
            MPI_Finalize();
#endif
            return WRITE_ERROR;
        }
#ifdef MOAB_HAVE_TEMPESTREMAP
    }
    delete remapper;
#endif

    if( !proc_id ) std::cerr << "Wrote \"" << out << "\"" << std::endl;
    if( print_times && !proc_id ) write_times( std::cout );

#ifdef MOAB_HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}

bool parse_id_list( const char* string, std::set< int >& results )
{
    bool okay   = true;
    char* mystr = strdup( string );
    for( const char* ptr = strtok( mystr, "," ); ptr; ptr = strtok( 0, "," ) )
    {
        char* endptr;
        long val = strtol( ptr, &endptr, 0 );
        if( endptr == ptr || val <= 0 )
        {
            std::cerr << "Not a valid id: " << ptr << std::endl;
            okay = false;
            break;
        }

        long val2 = val;
        if( *endptr == '-' )
        {
            const char* sptr = endptr + 1;
            val2             = strtol( sptr, &endptr, 0 );
            if( endptr == sptr || val2 <= 0 )
            {
                std::cerr << "Not a valid id: " << sptr << std::endl;
                okay = false;
                break;
            }
            if( val2 < val )
            {
                std::cerr << "Invalid id range: " << ptr << std::endl;
                okay = false;
                break;
            }
        }

        if( *endptr )
        {
            std::cerr << "Unexpected character: " << *endptr << std::endl;
            okay = false;
            break;
        }

        for( ; val <= val2; ++val )
            if( !results.insert( (int)val ).second ) std::cerr << "Warning: duplicate Id: " << val << std::endl;
    }

    free( mystr );
    return okay;
}

void print_id_list( const char* head, std::ostream& stream, const std::set< int >& list )
{
    stream << head << ": ";

    if( list.empty() )
    {
        stream << "(none)" << std::endl;
        return;
    }

    int start, prev;
    std::set< int >::const_iterator iter = list.begin();
    start = prev = *( iter++ );
    for( ;; )
    {
        if( iter == list.end() || *iter != 1 + prev )
        {
            stream << start;
            if( prev != start ) stream << '-' << prev;
            if( iter == list.end() ) break;
            stream << ", ";
            start = *iter;
        }
        prev = *( iter++ );
    }

    stream << std::endl;
}

static void print_time( int clk_per_sec, const char* prefix, clock_t ticks, std::ostream& stream )
{
    ticks *= clk_per_sec / 100;
    clock_t centi   = ticks % 100;
    clock_t seconds = ticks / 100;
    stream << prefix;
    if( seconds < 120 ) { stream << ( ticks / 100 ) << "." << centi << "s" << std::endl; }
    else
    {
        clock_t minutes = ( seconds / 60 ) % 60;
        clock_t hours   = ( seconds / 3600 );
        seconds %= 60;
        if( hours ) stream << hours << "h";
        if( minutes ) stream << minutes << "m";
        if( seconds || centi ) stream << seconds << "." << centi << "s";
        stream << " (" << ( ticks / 100 ) << "." << centi << "s)" << std::endl;
    }
}

clock_t usr_time, sys_time, abs_time;

#ifdef WIN32

void reset_times()
{
    abs_time = clock();
}

void write_times( std::ostream& stream )
{
    clock_t abs_tm = clock();
    print_time( CLOCKS_PER_SEC, "  ", abs_tm - abs_time, stream );
    abs_time = abs_tm;
}

#else

void reset_times()
{
    tms timebuf;
    abs_time = times( &timebuf );
    usr_time = timebuf.tms_utime;
    sys_time = timebuf.tms_stime;
}

void write_times( std::ostream& stream )
{
    clock_t usr_tm, sys_tm, abs_tm;
    tms timebuf;
    abs_tm = times( &timebuf );
    usr_tm = timebuf.tms_utime;
    sys_tm = timebuf.tms_stime;
    print_time( sysconf( _SC_CLK_TCK ), "  real:   ", abs_tm - abs_time, stream );
    print_time( sysconf( _SC_CLK_TCK ), "  user:   ", usr_tm - usr_time, stream );
    print_time( sysconf( _SC_CLK_TCK ), "  system: ", sys_tm - sys_time, stream );
    abs_time = abs_tm;
    usr_time = usr_tm;
    sys_time = sys_tm;
}

#endif

bool make_opts_string( std::vector< std::string > options, std::string& opts )
{
    opts.clear();
    if( options.empty() ) return true;

    // choose a separator character
    std::vector< std::string >::const_iterator i;
    char separator             = '\0';
    const char* alt_separators = ";+,:\t\n";
    for( const char* sep_ptr = alt_separators; *sep_ptr; ++sep_ptr )
    {
        bool seen = false;
        for( i = options.begin(); i != options.end(); ++i )
            if( i->find( *sep_ptr, 0 ) != std::string::npos )
            {
                seen = true;
                break;
            }
        if( !seen )
        {
            separator = *sep_ptr;
            break;
        }
    }
    if( !separator )
    {
        std::cerr << "Error: cannot find separator character for options string" << std::endl;
        return false;
    }
    if( separator != ';' )
    {
        opts = ";";
        opts += separator;
    }

    // concatenate options
    i = options.begin();
    opts += *i;
    for( ++i; i != options.end(); ++i )
    {
        opts += separator;
        opts += *i;
    }

    return true;
}

void list_formats( Interface* gMB )
{
    const char iface_name[] = "ReaderWriterSet";
    ErrorCode err;
    ReaderWriterSet* set = 0;
    ReaderWriterSet::iterator i;
    std::ostream& str = std::cout;

    // get ReaderWriterSet
    err = gMB->query_interface( set );
    if( err != MB_SUCCESS || !set )
    {
        std::cerr << "Internal error:  Interface \"" << iface_name << "\" not available.\n";
        exit( OTHER_ERROR );
    }

    // get field width for format description
    size_t w = 0;
    for( i = set->begin(); i != set->end(); ++i )
        if( i->description().length() > w ) w = i->description().length();

    // write table header
    str << "Format  " << std::setw( w ) << std::left << "Description"
        << "  Read  Write  File Name Suffixes\n"
        << "------  " << std::setw( w ) << std::setfill( '-' ) << "" << std::setfill( ' ' )
        << "  ----  -----  ------------------\n";

    // write table data
    for( i = set->begin(); i != set->end(); ++i )
    {
        std::vector< std::string > ext;
        i->get_extensions( ext );
        str << std::setw( 6 ) << i->name() << "  " << std::setw( w ) << std::left << i->description() << "  "
            << ( i->have_reader() ? " yes" : "  no" ) << "  " << ( i->have_writer() ? "  yes" : "   no" ) << " ";
        for( std::vector< std::string >::iterator j = ext.begin(); j != ext.end(); ++j )
            str << " " << *j;
        str << std::endl;
    }
    str << std::endl;

    gMB->release_interface( set );
    exit( 0 );
}

void remove_entities_from_sets( Interface* gMB, Range& dead_entities, Range& empty_sets )
{
    empty_sets.clear();
    Range sets;
    gMB->get_entities_by_type( 0, MBENTITYSET, sets );
    for( Range::iterator i = sets.begin(); i != sets.end(); ++i )
    {
        Range set_contents;
        gMB->get_entities_by_handle( *i, set_contents, false );
        set_contents = intersect( set_contents, dead_entities );
        gMB->remove_entities( *i, set_contents );
        set_contents.clear();
        gMB->get_entities_by_handle( *i, set_contents, false );
        if( set_contents.empty() ) empty_sets.insert( *i );
    }
}

void remove_from_vector( std::vector< EntityHandle >& vect, const Range& ents_to_remove )
{
    Range::const_iterator i;
    std::vector< EntityHandle >::iterator j;
    for( i = ents_to_remove.begin(); i != ents_to_remove.end(); ++i )
    {
        j = std::find( vect.begin(), vect.end(), *i );
        if( j != vect.end() ) vect.erase( j );
    }
}

std::string percent_subst( const std::string& s, int val )
{
    if( s.empty() ) return s;

    size_t j = s.find( '%' );
    if( j == std::string::npos ) return s;

    std::ostringstream st;
    st << s.substr( 0, j );
    st << val;

    size_t i;
    while( ( i = s.find( '%', j + 1 ) ) != std::string::npos )
    {
        st << s.substr( j, i - j );
        st << val;
        j = i;
    }
    st << s.substr( j + 1 );
    return st.str();
}

int process_partition_file( Interface* mb, std::string& metis_partition_file )
{
    // how many faces in the file ? how do we make sure it is an mpas file?
    // mpas atmosphere files can be downloaded from here
    // https://mpas-dev.github.io/atmosphere/atmosphere_meshes.html
    Range faces;
    ErrorCode rval = mb->get_entities_by_dimension( 0, 2, faces );MB_CHK_ERR( rval );
    std::cout << " MPAS model has " << faces.size() << " polygons\n";

    // read the partition file
    std::ifstream partfile;
    partfile.open( metis_partition_file.c_str() );
    std::string line;
    std::vector< int > parts;
    parts.resize( faces.size(), -1 );
    int i = 0;
    if( partfile.is_open() )
    {
        while( getline( partfile, line ) )
        {
            // cout << line << '\n';
            parts[i++] = atoi( line.c_str() );
            if( i > (int)faces.size() )
            {
                std::cerr << " too many lines in partition file \n. bail out \n";
                return 1;
            }
        }
        partfile.close();
    }
    std::vector< int >::iterator pmax = max_element( parts.begin(), parts.end() );
    std::vector< int >::iterator pmin = min_element( parts.begin(), parts.end() );
    if( *pmin <= -1 )
    {
        std::cerr << " partition file is incomplete, *pmin is -1 !! \n";
        return 1;
    }
    std::cout << " partitions range: " << *pmin << " " << *pmax << "\n";
    Tag part_set_tag;
    int dum_id = -1;
    rval = mb->tag_get_handle( "PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, part_set_tag, MB_TAG_SPARSE | MB_TAG_CREAT,
                               &dum_id );MB_CHK_ERR( rval );

    // get any sets already with this tag, and clear them
    // remove the parallel partition sets if they exist
    Range tagged_sets;
    rval = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &part_set_tag, NULL, 1, tagged_sets, Interface::UNION );MB_CHK_ERR( rval );
    if( !tagged_sets.empty() )
    {
        rval = mb->clear_meshset( tagged_sets );MB_CHK_ERR( rval );
        rval = mb->tag_delete_data( part_set_tag, tagged_sets );MB_CHK_ERR( rval );
    }
    Tag gid;
    rval = mb->tag_get_handle( "GLOBAL_ID", gid );MB_CHK_ERR( rval );
    int num_sets = *pmax + 1;
    if( *pmin != 0 )
    {
        std::cout << " problem reading parts; min is not 0 \n";
        return 1;
    }
    for( i = 0; i < num_sets; i++ )
    {
        EntityHandle new_set;
        rval = mb->create_meshset( MESHSET_SET, new_set );MB_CHK_ERR( rval );
        tagged_sets.insert( new_set );
    }
    int* dum_ids = new int[num_sets];
    for( i = 0; i < num_sets; i++ )
        dum_ids[i] = i;

    rval = mb->tag_set_data( part_set_tag, tagged_sets, dum_ids );MB_CHK_ERR( rval );
    delete[] dum_ids;

    std::vector< int > gids;
    int num_faces = (int)faces.size();
    gids.resize( num_faces );
    rval = mb->tag_get_data( gid, faces, &gids[0] );MB_CHK_ERR( rval );

    for( int j = 0; j < num_faces; j++ )
    {
        int eid         = gids[j];
        EntityHandle eh = faces[j];
        int partition   = parts[eid - 1];
        if( partition < 0 || partition >= num_sets )
        {
            std::cout << " wrong partition number \n";
            return 1;
        }
        rval = mb->add_entities( tagged_sets[partition], &eh, 1 );MB_CHK_ERR( rval );
    }
    return 0;
}
