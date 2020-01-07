/*
 * =====================================================================================
 *
 *       Filename:  TempestOnlineMap.hpp
 *
 *    Description:  Interface to the TempestRemap library to compute the consistent,
 *                  and accurate high-order conservative remapping weights for overlap
 *                  grids on the sphere in climate simulations.
 *
 *         Author:  Vijay S. Mahadevan (vijaysm), mahadevan@anl.gov
 *
 * =====================================================================================
 */

#include "Announce.h"
#include "DataArray3D.h"
#include "FiniteElementTools.h"
#include "SparseMatrix.h"
#include "STLStringHelper.h"

#include "moab/Remapping/TempestOnlineMap.hpp"
#include "DebugOutput.hpp"

#include <fstream>
#include <cmath>
#include <cstdlib>

///////////////////////////////////////////////////////////////////////////////

// #define VERBOSE
// #define VVERBOSE

void LinearRemapFVtoGLL(
    const Mesh & meshInput,
    const Mesh & meshOutput,
    const Mesh & meshOverlap,
    const DataArray3D<int> & dataGLLNodes,
    const DataArray3D<double> & dataGLLJacobian,
    const DataArray1D<double> & dataGLLNodalArea,
    int nOrder,
    OfflineMap & mapRemap,
    int nMonotoneType,
    bool fContinuous,
    bool fNoConservation
);

void LinearRemapFVtoGLL_Volumetric(
    const Mesh & meshInput,
    const Mesh & meshOutput,
    const Mesh & meshOverlap,
    const DataArray3D<int> & dataGLLNodes,
    const DataArray3D<double> & dataGLLJacobian,
    const DataArray1D<double> & dataGLLNodalArea,
    int nOrder,
    OfflineMap & mapRemap,
    int nMonotoneType,
    bool fContinuous,
    bool fNoConservation
);

///////////////////////////////////////////////////////////////////////////////

moab::TempestOnlineMap::TempestOnlineMap ( moab::TempestRemapper* remapper ) : OfflineMap(), m_remapper ( remapper )
{
    // Get the references for the MOAB core objects
    m_interface = m_remapper->get_interface();
#ifdef MOAB_HAVE_MPI
    m_pcomm = m_remapper->get_parallel_communicator();
#endif

    // Update the references to the meshes
    m_meshInput = remapper->GetMesh ( moab::Remapper::SourceMesh );
    m_meshInputCov = remapper->GetCoveringMesh();
    m_meshOutput = remapper->GetMesh ( moab::Remapper::TargetMesh );
    m_meshOverlap = remapper->GetMesh ( moab::Remapper::IntersectedMesh );

    is_parallel = false;
    is_root = true;
    rank = 0;
    size = 1;
#ifdef MOAB_HAVE_MPI
    int flagInit;
    MPI_Initialized( &flagInit );
    if (flagInit) {
        is_parallel = true;
        assert(m_pcomm != NULL);
        rank = m_pcomm->rank();
        size = m_pcomm->size();
        is_root = (rank == 0);
    }
#endif

    // Compute and store the total number of source and target DoFs corresponding
    // to number of rows and columns in the mapping.

    // Initialize dimension information from file
    std::vector<std::string> dimNames(1);
    std::vector<int> dimSizes(1);
    dimNames[0] = "num_elem";

    dimSizes[0] = m_meshInputCov->faces.size();
    this->InitializeSourceDimensions(dimNames, dimSizes);
    dimSizes[0] = m_meshOutput->faces.size();
    this->InitializeTargetDimensions(dimNames, dimSizes);

    // Build a matrix of source and target discretization so that we know how to assign
    // the global DoFs in parallel for the mapping weights
    // For example, FV->FV: rows X cols = faces_source X faces_target
}

///////////////////////////////////////////////////////////////////////////////

moab::TempestOnlineMap::~TempestOnlineMap()
{
    m_interface = NULL;
#ifdef MOAB_HAVE_MPI
    m_pcomm = NULL;
#endif
    m_meshInput = NULL;
    m_meshOutput = NULL;
    m_meshOverlap = NULL;
}

///////////////////////////////////////////////////////////////////////////////

static void ParseVariableList (
    const std::string & strVariables,
    std::vector< std::string > & vecVariableStrings
)
{
    unsigned iVarBegin = 0;
    unsigned iVarCurrent = 0;

    // Parse variable name
    for ( ;; )
    {
        if ( ( iVarCurrent >= strVariables.length() ) ||
                ( strVariables[iVarCurrent] == ',' ) ||
                ( strVariables[iVarCurrent] == ' ' )
           )
        {
            if ( iVarCurrent == iVarBegin )
            {
                if ( iVarCurrent >= strVariables.length() )
                {
                    break;
                }
                continue;
            }

            vecVariableStrings.push_back (
                strVariables.substr ( iVarBegin, iVarCurrent - iVarBegin ) );

            iVarBegin = iVarCurrent + 1;
        }

        iVarCurrent++;
    }
}

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::SetDOFmapTags(const std::string srcDofTagName, const std::string tgtDofTagName)
{
    moab::ErrorCode rval;

    int tagSize = 0;
    tagSize = (m_eInputType == DiscretizationType_FV ? 1 : m_nDofsPEl_Src*m_nDofsPEl_Src);
    rval = m_interface->tag_get_handle ( srcDofTagName.c_str(), tagSize, MB_TYPE_INTEGER,
                             this->m_dofTagSrc, MB_TAG_ANY );

    if (rval == moab::MB_TAG_NOT_FOUND && m_eInputType != DiscretizationType_FV)
    {
        int ntot_elements = 0, nelements = m_remapper->m_source_entities.size();
#ifdef MOAB_HAVE_MPI
        int ierr = MPI_Allreduce(&nelements, &ntot_elements, 1, MPI_INT, MPI_SUM, m_pcomm->comm());
        if (ierr !=0) MB_CHK_SET_ERR(MB_FAILURE, "MPI_Allreduce failed to get total source elements");
#else
        ntot_elements = nelements;
#endif

        rval = m_remapper->GenerateCSMeshMetadata(ntot_elements, 
                                                m_remapper->m_covering_source_entities, 
                                                &m_remapper->m_source_entities, 
                                                srcDofTagName, m_nDofsPEl_Src);MB_CHK_ERR(rval);

        rval = m_interface->tag_get_handle ( srcDofTagName.c_str(), m_nDofsPEl_Src * m_nDofsPEl_Src, 
                                        MB_TYPE_INTEGER,
                                        this->m_dofTagSrc, MB_TAG_ANY);MB_CHK_ERR(rval);
    }
    else MB_CHK_ERR(rval);

    tagSize = (m_eOutputType == DiscretizationType_FV ? 1 : m_nDofsPEl_Dest*m_nDofsPEl_Dest);
    rval = m_interface->tag_get_handle ( tgtDofTagName.c_str(), tagSize, MB_TYPE_INTEGER,
                             this->m_dofTagDest, MB_TAG_ANY );
    if (rval == moab::MB_TAG_NOT_FOUND && m_eOutputType != DiscretizationType_FV)
    {
        int ntot_elements = 0, nelements = m_remapper->m_target_entities.size();
#ifdef MOAB_HAVE_MPI
        int ierr = MPI_Allreduce(&nelements, &ntot_elements, 1, MPI_INT, MPI_SUM, m_pcomm->comm());
        if (ierr !=0) MB_CHK_SET_ERR(MB_FAILURE, "MPI_Allreduce failed to get total source elements");
#else
        ntot_elements = nelements;
#endif

        rval = m_remapper->GenerateCSMeshMetadata(ntot_elements, 
                                                    m_remapper->m_target_entities, NULL, 
                                                    tgtDofTagName, m_nDofsPEl_Dest);MB_CHK_ERR(rval);

        rval = m_interface->tag_get_handle ( tgtDofTagName.c_str(), m_nDofsPEl_Dest * m_nDofsPEl_Dest, 
                                        MB_TYPE_INTEGER,
                                        this->m_dofTagDest, MB_TAG_ANY);MB_CHK_ERR(rval);
    }
    else MB_CHK_ERR(rval);

    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::SetDOFmapAssociation(DiscretizationType srcType, bool isSrcContinuous, 
                                                            DataArray3D<int>* srcdataGLLNodes, DataArray3D<int>* srcdataGLLNodesSrc,
                                                            DiscretizationType destType, bool isTgtContinuous, DataArray3D<int>* tgtdataGLLNodes)
{
    moab::ErrorCode rval;
    std::vector<bool> dgll_cgll_row_ldofmap, dgll_cgll_col_ldofmap, dgll_cgll_covcol_ldofmap;
    std::vector<unsigned> src_soln_gdofs, locsrc_soln_gdofs, tgt_soln_gdofs;

    // We are assuming that these are element based tags that are sized: np * np
    m_srcDiscType = srcType;
    m_destDiscType = destType;

    bool vprint = is_root && false;

#ifdef VVERBOSE
    {
        src_soln_gdofs.resize(m_remapper->m_covering_source_entities.size()*m_nDofsPEl_Src*m_nDofsPEl_Src, -1);
        rval = m_interface->tag_get_data ( m_dofTagSrc, m_remapper->m_covering_source_entities, &src_soln_gdofs[0] );MB_CHK_ERR(rval);
        locsrc_soln_gdofs.resize(m_remapper->m_source_entities.size()*m_nDofsPEl_Src*m_nDofsPEl_Src);
        rval = m_interface->tag_get_data ( m_dofTagSrc, m_remapper->m_source_entities, &locsrc_soln_gdofs[0] );MB_CHK_ERR(rval);
        tgt_soln_gdofs.resize(m_remapper->m_target_entities.size()*m_nDofsPEl_Dest*m_nDofsPEl_Dest);
        rval = m_interface->tag_get_data ( m_dofTagDest, m_remapper->m_target_entities, &tgt_soln_gdofs[0] );MB_CHK_ERR(rval);

        if (is_root)
        {
            {
                std::ofstream output_file ( "sourcecov-gids-0.txt" );
                output_file << "I, GDOF\n";
                for (unsigned i=0; i < src_soln_gdofs.size(); ++i)
                    output_file << i << ", " << src_soln_gdofs[i] << "\n";

                output_file << "ELEMID, IDOF, LDOF, GDOF, NDOF\n";
                m_nTotDofs_SrcCov=0;
                if (isSrcContinuous) dgll_cgll_covcol_ldofmap.resize (m_remapper->m_covering_source_entities.size() * m_nDofsPEl_Src * m_nDofsPEl_Src, false);
                for ( unsigned j = 0; j < m_remapper->m_covering_source_entities.size(); j++ )
                {
                    for ( int p = 0; p < m_nDofsPEl_Src; p++ )
                    {
                        for ( int q = 0; q < m_nDofsPEl_Src; q++)
                        {
                            const int ldof = (*srcdataGLLNodes)[p][q][j] - 1;
                            const int idof = j * m_nDofsPEl_Src * m_nDofsPEl_Src + p * m_nDofsPEl_Src + q;
                            if ( isSrcContinuous && !dgll_cgll_covcol_ldofmap[ldof] ) {
                                m_nTotDofs_SrcCov++;
                                dgll_cgll_covcol_ldofmap[ldof] = true;
                            }
                            output_file << m_remapper->lid_to_gid_covsrc[j] << ", " <<  idof << ", " << ldof << ", " << src_soln_gdofs[idof] << ", " << m_nTotDofs_SrcCov << "\n";
                        }
                    }
                }
                output_file.flush(); // required here
                output_file.close();
                dgll_cgll_covcol_ldofmap.clear();
            }

            {
                std::ofstream output_file ( "source-gids-0.txt" );
                output_file << "I, GDOF\n";
                for (unsigned i=0; i < locsrc_soln_gdofs.size(); ++i)
                    output_file << i << ", " << locsrc_soln_gdofs[i] << "\n";

                output_file << "ELEMID, IDOF, LDOF, GDOF, NDOF\n";
                m_nTotDofs_Src=0;
                if (isSrcContinuous) dgll_cgll_col_ldofmap.resize (m_remapper->m_source_entities.size() * m_nDofsPEl_Src * m_nDofsPEl_Src, false);
                for ( unsigned j = 0; j < m_remapper->m_source_entities.size(); j++ )
                {
                    for ( int p = 0; p < m_nDofsPEl_Src; p++ )
                    {
                        for ( int q = 0; q < m_nDofsPEl_Src; q++)
                        {
                            const int ldof = (*srcdataGLLNodesSrc)[p][q][j] - 1;
                            const int idof = j * m_nDofsPEl_Src * m_nDofsPEl_Src + p * m_nDofsPEl_Src + q;
                            if ( isSrcContinuous && !dgll_cgll_col_ldofmap[ldof] ) {
                                m_nTotDofs_Src++;
                                dgll_cgll_col_ldofmap[ldof] = true;
                            }
                            output_file << m_remapper->lid_to_gid_src[j] << ", " <<  idof << ", " << ldof << ", " << locsrc_soln_gdofs[idof] << ", " << m_nTotDofs_Src << "\n";
                        }
                    }
                }
                output_file.flush(); // required here
                output_file.close();
                dgll_cgll_col_ldofmap.clear();
            }

            {
                std::ofstream output_file ( "target-gids-0.txt" );
                output_file << "I, GDOF\n";
                for (unsigned i=0; i < tgt_soln_gdofs.size(); ++i)
                    output_file << i << ", " << tgt_soln_gdofs[i] << "\n";

                output_file << "ELEMID, IDOF, GDOF, NDOF\n";
                m_nTotDofs_Dest=0;

                for (unsigned i=0; i < tgt_soln_gdofs.size(); ++i) {
                    output_file << m_remapper->lid_to_gid_tgt[i] << ", " <<  i << ", " << tgt_soln_gdofs[i] << ", " << m_nTotDofs_Dest << "\n";
                    m_nTotDofs_Dest++;
                }

                output_file.flush(); // required here
                output_file.close();
            }
        }
        else
        {
            {
                std::ofstream output_file ( "sourcecov-gids-1.txt" );
                output_file << "I, GDOF\n";
                for (unsigned i=0; i < src_soln_gdofs.size(); ++i)
                    output_file << i << ", " << src_soln_gdofs[i] << "\n";

                output_file << "ELEMID, IDOF, LDOF, GDOF, NDOF\n";
                m_nTotDofs_SrcCov=0;
                if (isSrcContinuous) dgll_cgll_covcol_ldofmap.resize (m_remapper->m_covering_source_entities.size() * m_nDofsPEl_Src * m_nDofsPEl_Src, false);
                for ( unsigned j = 0; j < m_remapper->m_covering_source_entities.size(); j++ )
                {
                    for ( int p = 0; p < m_nDofsPEl_Src; p++ )
                    {
                        for ( int q = 0; q < m_nDofsPEl_Src; q++)
                        {
                            const int ldof = (*srcdataGLLNodes)[p][q][j] - 1;
                            const int idof = j * m_nDofsPEl_Src * m_nDofsPEl_Src + p * m_nDofsPEl_Src + q;
                            if ( isSrcContinuous && !dgll_cgll_covcol_ldofmap[ldof] ) {
                                m_nTotDofs_SrcCov++;
                                dgll_cgll_covcol_ldofmap[ldof] = true;
                            }
                            output_file << m_remapper->lid_to_gid_covsrc[j] << ", " <<  idof << ", " << ldof << ", " << src_soln_gdofs[idof] << ", " << m_nTotDofs_SrcCov << "\n";
                        }
                    }
                }
                output_file.flush(); // required here
                output_file.close();
                dgll_cgll_covcol_ldofmap.clear();
            }

            {
                std::ofstream output_file ( "source-gids-1.txt" );
                output_file << "I, GDOF\n";
                for (unsigned i=0; i < locsrc_soln_gdofs.size(); ++i)
                    output_file << i << ", " << locsrc_soln_gdofs[i] << "\n";

                output_file << "ELEMID, IDOF, LDOF, GDOF, NDOF\n";
                m_nTotDofs_Src=0;
                if (isSrcContinuous) dgll_cgll_col_ldofmap.resize (m_remapper->m_source_entities.size() * m_nDofsPEl_Src * m_nDofsPEl_Src, false);
                for ( unsigned j = 0; j < m_remapper->m_source_entities.size(); j++ )
                {
                    for ( int p = 0; p < m_nDofsPEl_Src; p++ )
                    {
                        for ( int q = 0; q < m_nDofsPEl_Src; q++)
                        {
                            const int ldof = (*srcdataGLLNodesSrc)[p][q][j] - 1;
                            const int idof = j * m_nDofsPEl_Src * m_nDofsPEl_Src + p * m_nDofsPEl_Src + q;
                            if ( isSrcContinuous && !dgll_cgll_col_ldofmap[ldof] ) {
                                m_nTotDofs_Src++;
                                dgll_cgll_col_ldofmap[ldof] = true;
                            }
                            output_file << m_remapper->lid_to_gid_src[j] << ", " <<  idof << ", " << ldof << ", " << locsrc_soln_gdofs[idof] << ", " << m_nTotDofs_Src << "\n";
                        }
                    }
                }
                output_file.flush(); // required here
                output_file.close();
                dgll_cgll_col_ldofmap.clear();
            }

            {
                std::ofstream output_file ( "target-gids-1.txt" );
                output_file << "I, GDOF\n";
                for (unsigned i=0; i < tgt_soln_gdofs.size(); ++i)
                    output_file << i << ", " << tgt_soln_gdofs[i] << "\n";

                output_file << "ELEMID, IDOF, GDOF, NDOF\n";
                m_nTotDofs_Dest=0;

                for (unsigned i=0; i < tgt_soln_gdofs.size(); ++i) {
                    output_file << m_remapper->lid_to_gid_tgt[i] << ", " <<  i << ", " << tgt_soln_gdofs[i] << ", " << m_nTotDofs_Dest << "\n";
                    m_nTotDofs_Dest++;
                }

                output_file.flush(); // required here
                output_file.close();
            }
        }
    }
#endif

    // Now compute the mapping and store it for the covering mesh
    int srcTagSize = (m_eInputType == DiscretizationType_FV ? 1 : m_nDofsPEl_Src * m_nDofsPEl_Src);
    if (m_remapper->point_cloud_source)
    {
        assert(m_nDofsPEl_Src == 1);
        col_dofmap.resize (m_remapper->m_covering_source_vertices.size(), UINT_MAX);
        col_ldofmap.resize (m_remapper->m_covering_source_vertices.size(), UINT_MAX);
        col_gdofmap.resize (m_remapper->m_covering_source_vertices.size(), UINT_MAX);
        src_soln_gdofs.resize(m_remapper->m_covering_source_vertices.size(), UINT_MAX);
        rval = m_interface->tag_get_data ( m_dofTagSrc, m_remapper->m_covering_source_vertices, &src_soln_gdofs[0] );MB_CHK_ERR(rval);
        srcTagSize = 1;
    }
    else
    {
        col_dofmap.resize (m_remapper->m_covering_source_entities.size() * srcTagSize, UINT_MAX);
        col_ldofmap.resize (m_remapper->m_covering_source_entities.size() * srcTagSize, UINT_MAX);
        col_gdofmap.resize (m_remapper->m_covering_source_entities.size() * srcTagSize, UINT_MAX);
        src_soln_gdofs.resize(m_remapper->m_covering_source_entities.size() * srcTagSize, UINT_MAX);
        rval = m_interface->tag_get_data ( m_dofTagSrc, m_remapper->m_covering_source_entities, &src_soln_gdofs[0] );MB_CHK_ERR(rval);
    }
    m_nTotDofs_SrcCov = 0;
    if (srcdataGLLNodes == NULL) { /* we only have a mapping for elements as DoFs */
        for (unsigned i=0; i < col_dofmap.size(); ++i) {
            col_dofmap[i] = i;
            col_ldofmap[i] = i;
            col_gdofmap[i] = src_soln_gdofs[i]-1;
            if (vprint) std::cout << "Col: " << i << ", " << col_gdofmap[i] << "\n";
            m_nTotDofs_SrcCov++;
        }
    }
    else {
        if (isSrcContinuous) dgll_cgll_covcol_ldofmap.resize (m_remapper->m_covering_source_entities.size() * srcTagSize, false);
        // Put these remap coefficients into the SparseMatrix map
        for ( unsigned j = 0; j < m_remapper->m_covering_source_entities.size(); j++ )
        {
            for ( int p = 0; p < m_nDofsPEl_Src; p++ )
            {
                for ( int q = 0; q < m_nDofsPEl_Src; q++)
                {
                    const int ldof = (*srcdataGLLNodes)[p][q][j] - 1;
                    const int idof = j * srcTagSize + p * m_nDofsPEl_Src + q;
                    if ( isSrcContinuous && !dgll_cgll_covcol_ldofmap[ldof] ) {
                        m_nTotDofs_SrcCov++;
                        dgll_cgll_covcol_ldofmap[ldof] = true;
                    }
                    if ( !isSrcContinuous ) m_nTotDofs_SrcCov++;
                    col_dofmap[ idof ] = ldof;
                    col_ldofmap[ ldof ] = idof;
                    assert(src_soln_gdofs[idof] > 0);
                    col_gdofmap[ idof ] = src_soln_gdofs[idof] - 1;
                    if (vprint) std::cout << "Col: " << m_remapper->lid_to_gid_covsrc[j] << ", " <<  idof << ", " << ldof << ", " << col_gdofmap[idof] << ", " << m_nTotDofs_SrcCov << "\n";
                }
            }
        }
    }

    if (m_remapper->point_cloud_source)
    {
        assert(m_nDofsPEl_Src == 1);
        srccol_dofmap.resize (m_remapper->m_source_vertices.size(), UINT_MAX);
        srccol_ldofmap.resize (m_remapper->m_source_vertices.size(), UINT_MAX);
        srccol_gdofmap.resize (m_remapper->m_source_vertices.size(), UINT_MAX);
        locsrc_soln_gdofs.resize(m_remapper->m_source_vertices.size(), UINT_MAX);
        rval = m_interface->tag_get_data ( m_dofTagSrc, m_remapper->m_source_vertices, &locsrc_soln_gdofs[0] );MB_CHK_ERR(rval);
    }
    else
    {
        srccol_dofmap.resize (m_remapper->m_source_entities.size() * srcTagSize, UINT_MAX);
        srccol_ldofmap.resize (m_remapper->m_source_entities.size() * srcTagSize, UINT_MAX);
        srccol_gdofmap.resize (m_remapper->m_source_entities.size() * srcTagSize, UINT_MAX);
        locsrc_soln_gdofs.resize(m_remapper->m_source_entities.size() * srcTagSize, UINT_MAX);
        rval = m_interface->tag_get_data ( m_dofTagSrc, m_remapper->m_source_entities, &locsrc_soln_gdofs[0] );MB_CHK_ERR(rval);
    }

    // Now compute the mapping and store it for the original source mesh
    m_nTotDofs_Src = 0;
    if (srcdataGLLNodesSrc == NULL) { /* we only have a mapping for elements as DoFs */
        for (unsigned i=0; i < srccol_dofmap.size(); ++i) {
            srccol_dofmap[i] = i;
            srccol_ldofmap[i] = i;
            srccol_gdofmap[i] = locsrc_soln_gdofs[i]-1;
            m_nTotDofs_Src++;
        }
    }
    else {
        if (isSrcContinuous) dgll_cgll_col_ldofmap.resize(m_remapper->m_source_entities.size() * srcTagSize, false);
        // Put these remap coefficients into the SparseMatrix map
        for ( unsigned j = 0; j < m_remapper->m_source_entities.size(); j++ )
        {
            for ( int p = 0; p < m_nDofsPEl_Src; p++ )
            {
                for ( int q = 0; q < m_nDofsPEl_Src; q++ )
                {
                    const int ldof = (*srcdataGLLNodesSrc)[p][q][j] - 1;
                    const int idof = j * srcTagSize + p * m_nDofsPEl_Src + q;
                    if ( isSrcContinuous && !dgll_cgll_col_ldofmap[ldof] ) {
                        m_nTotDofs_Src++;
                        dgll_cgll_col_ldofmap[ldof] = true;
                    }
                    if ( !isSrcContinuous ) m_nTotDofs_Src++;
                    srccol_dofmap[ idof ] = ldof;
                    srccol_ldofmap[ ldof ] = idof;
                    srccol_gdofmap[ idof ] = locsrc_soln_gdofs[idof] - 1;
                }
            }
        }
    }

    int tgtTagSize = (m_eOutputType == DiscretizationType_FV ? 1 : m_nDofsPEl_Dest * m_nDofsPEl_Dest);
    if (m_remapper->point_cloud_target)
    {
        assert(m_nDofsPEl_Dest == 1);
        row_dofmap.resize (m_remapper->m_target_vertices.size(), UINT_MAX);
        row_ldofmap.resize (m_remapper->m_target_vertices.size(), UINT_MAX);
        row_gdofmap.resize (m_remapper->m_target_vertices.size(), UINT_MAX);
        tgt_soln_gdofs.resize(m_remapper->m_target_vertices.size(), UINT_MAX);
        rval = m_interface->tag_get_data ( m_dofTagDest, m_remapper->m_target_vertices, &tgt_soln_gdofs[0] );MB_CHK_ERR(rval);
        tgtTagSize = 1;
    }
    else
    {
        row_dofmap.resize (m_remapper->m_target_entities.size() * tgtTagSize, UINT_MAX);
        row_ldofmap.resize (m_remapper->m_target_entities.size() * tgtTagSize, UINT_MAX);
        row_gdofmap.resize (m_remapper->m_target_entities.size() * tgtTagSize, UINT_MAX);
        tgt_soln_gdofs.resize(m_remapper->m_target_entities.size() * tgtTagSize, UINT_MAX);
        rval = m_interface->tag_get_data ( m_dofTagDest, m_remapper->m_target_entities, &tgt_soln_gdofs[0] );MB_CHK_ERR(rval);
    }

    // Now compute the mapping and store it for the target mesh
    // To access the GID for each row: row_gdofmap [ row_ldofmap [ 0 : local_ndofs ] ] = GDOF
    m_nTotDofs_Dest = 0;
    if (tgtdataGLLNodes == NULL) { /* we only have a mapping for elements as DoFs */
        for (unsigned i=0; i < row_dofmap.size(); ++i) {
            row_dofmap[i] = i;
            row_ldofmap[i] = i;
            row_gdofmap[i] = tgt_soln_gdofs[i]-1;
            if (vprint) std::cout << "Row: " << i << ", " << row_gdofmap[i] << "\n";
            m_nTotDofs_Dest++;
        }
    }
    else {
        if (isTgtContinuous) dgll_cgll_row_ldofmap.resize (m_remapper->m_target_entities.size() * tgtTagSize, false);
        // Put these remap coefficients into the SparseMatrix map
        for ( unsigned j = 0; j < m_remapper->m_target_entities.size(); j++ )
        {
            for ( int p = 0; p < m_nDofsPEl_Dest; p++ )
            {
                for ( int q = 0; q < m_nDofsPEl_Dest; q++ )
                {
                    const int ldof = (*tgtdataGLLNodes)[p][q][j] - 1;
                    const int idof = j * tgtTagSize + p * m_nDofsPEl_Dest + q;
                    if ( isTgtContinuous && !dgll_cgll_row_ldofmap[ldof] ) {
                        m_nTotDofs_Dest++;
                        dgll_cgll_row_ldofmap[ldof] = true;
                    }
                    if ( !isTgtContinuous ) m_nTotDofs_Dest++;
                    row_dofmap[ idof ] = ldof;
                    row_ldofmap[ ldof ] = idof;
                    row_gdofmap[ idof ] = tgt_soln_gdofs[idof] - 1;
                    if (vprint) std::cout << "Row: " << idof << ", " << ldof << ", " << tgt_soln_gdofs[idof] - 1 << "\n";
                }
            }
        }
    }

    // Let us also allocate the local representation of the sparse matrix
#if defined(MOAB_HAVE_EIGEN) && defined(VERBOSE)
    if (vprint)
    {
        std::cout << "[" << rank << "]" << "DoFs: row = " << m_nTotDofs_Dest << ", " << row_dofmap.size() << ", col = " << m_nTotDofs_Src << ", " << m_nTotDofs_SrcCov << ", " << col_dofmap.size() << "\n";
        // std::cout << "Max col_dofmap: " << maxcol << ", Min col_dofmap" << mincol << "\n";
    }
#endif

    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::GenerateRemappingWeights ( std::string strInputType, std::string strOutputType,
        const int nPin, const int nPout,
        bool fBubble, int fMonotoneTypeID,
        bool fVolumetric, bool fNoConservation, bool fNoCheck,
        const std::string srcDofTagName, const std::string tgtDofTagName,
        const std::string strVariables,
        const std::string strInputData, const std::string strOutputData,
        const std::string strNColName, const bool fOutputDouble,
        const std::string strPreserveVariables, const bool fPreserveAll, const double dFillValueOverride,
        const bool fInputConcave, const bool fOutputConcave )
{
    NcError error ( NcError::silent_nonfatal );

    moab::DebugOutput dbgprint ( std::cout, rank, 0 );
    dbgprint.set_prefix("[TempestOnlineMap]: ");
    moab::ErrorCode rval;

    const bool m_bPointCloudSource = (m_remapper->point_cloud_source);
    const bool m_bPointCloudTarget = (m_remapper->point_cloud_target);
    const bool m_bPointCloud = m_bPointCloudSource || m_bPointCloudTarget;

    try
    {
        // Check command line parameters (data arguments)
        if ( ( strInputData != "" ) && ( strOutputData == "" ) )
        {
            _EXCEPTIONT ( "--in_data specified without --out_data" );
        }
        if ( ( strInputData == "" ) && ( strOutputData != "" ) )
        {
            _EXCEPTIONT ( "--out_data specified without --in_data" );
        }

        // Check command line parameters (data type arguments)
        STLStringHelper::ToLower ( strInputType );
        STLStringHelper::ToLower ( strOutputType );

        DiscretizationType eInputType;
        DiscretizationType eOutputType;

        if ( strInputType == "fv" )
        {
            eInputType = DiscretizationType_FV;
        }
        else if ( strInputType == "cgll" )
        {
            eInputType = DiscretizationType_CGLL;
        }
        else if ( strInputType == "dgll" )
        {
            eInputType = DiscretizationType_DGLL;
        }
        else if ( strInputType == "pcloud" )
        {
            eInputType = DiscretizationType_PCLOUD;
        }
        else
        {
            _EXCEPTION1 ( "Invalid \"in_type\" value (%s), expected [fv|cgll|dgll]",
                          strInputType.c_str() );
        }

        if ( strOutputType == "fv" )
        {
            eOutputType = DiscretizationType_FV;
        }
        else if ( strOutputType == "cgll" )
        {
            eOutputType = DiscretizationType_CGLL;
        }
        else if ( strOutputType == "dgll" )
        {
            eOutputType = DiscretizationType_DGLL;
        }
        else if ( strOutputType == "pcloud" )
        {
            eOutputType = DiscretizationType_PCLOUD;
        }
        else
        {
            _EXCEPTION1 ( "Invalid \"out_type\" value (%s), expected [fv|cgll|dgll]",
                          strOutputType.c_str() );
        }

        // Monotonicity flags
        int nMonotoneType = fMonotoneTypeID;

        // Parse variable list
        std::vector< std::string > vecVariableStrings;
        ParseVariableList ( strVariables, vecVariableStrings );

        // Parse preserve variable list
        std::vector< std::string > vecPreserveVariableStrings;
        ParseVariableList ( strPreserveVariables, vecPreserveVariableStrings );

        if ( fPreserveAll && ( vecPreserveVariableStrings.size() != 0 ) )
        {
            _EXCEPTIONT ( "--preserveall and --preserve cannot both be specified" );
        }

        m_nDofsPEl_Src = nPin;
        m_nDofsPEl_Dest = nPout;

        m_bConserved = !fNoConservation;
        m_iMonotonicity = fMonotoneTypeID;
        m_eInputType = eInputType;
        m_eOutputType = eOutputType;

        rval = SetDOFmapTags(srcDofTagName, tgtDofTagName);MB_CHK_ERR(rval);

        double dTotalAreaInput = 0.0, dTotalAreaOutput = 0.0;
        if (!m_bPointCloudSource)
        {
            // Calculate Face areas
            if ( is_root ) dbgprint.printf ( 0, "Calculating input mesh Face areas\n" );
            double dTotalAreaInput_loc = m_meshInput->CalculateFaceAreas(fInputConcave);
            dTotalAreaInput = dTotalAreaInput_loc;
    #ifdef MOAB_HAVE_MPI
            if (m_pcomm) MPI_Reduce ( &dTotalAreaInput_loc, &dTotalAreaInput, 1, MPI_DOUBLE, MPI_SUM, 0, m_pcomm->comm() );
    #endif
            if ( is_root ) dbgprint.printf ( 0, "Input Mesh Geometric Area: %1.15e\n", dTotalAreaInput );

            // Input mesh areas
            m_meshInputCov->CalculateFaceAreas(fInputConcave);
            if ( eInputType == DiscretizationType_FV )
            {
                this->SetSourceAreas ( m_meshInputCov->vecFaceArea );
                if (m_meshInputCov->vecMask.IsAttached()) {
                    this->SetSourceMask(m_meshInputCov->vecMask);
                }
            }
        }

        if (!m_bPointCloudTarget)
        {
            // Calculate Face areas
            if ( is_root ) dbgprint.printf ( 0, "Calculating output mesh Face areas\n" );
            double dTotalAreaOutput_loc = m_meshOutput->CalculateFaceAreas(fOutputConcave);
            dTotalAreaOutput = dTotalAreaOutput_loc;
    #ifdef MOAB_HAVE_MPI
            if (m_pcomm) MPI_Reduce ( &dTotalAreaOutput_loc, &dTotalAreaOutput, 1, MPI_DOUBLE, MPI_SUM, 0, m_pcomm->comm() );
    #endif
            if ( is_root ) dbgprint.printf ( 0, "Output Mesh Geometric Area: %1.15e\n", dTotalAreaOutput );

            // Output mesh areas
            if ( eOutputType == DiscretizationType_FV )
            {
                this->SetTargetAreas ( m_meshOutput->vecFaceArea );
                if (m_meshOutput->vecMask.IsAttached()) {
                    this->SetTargetMask(m_meshOutput->vecMask);
                }
            }
        }

        if (!m_bPointCloud)
        {
            // Verify that overlap mesh is in the correct order
            int ixSourceFaceMax = ( -1 );
            int ixTargetFaceMax = ( -1 );

            if ( m_meshOverlap->vecSourceFaceIx.size() !=
                    m_meshOverlap->vecTargetFaceIx.size()
               )
            {
                _EXCEPTIONT ( "Invalid overlap mesh:\n"
                              "    Possible mesh file corruption?" );
            }

            for ( unsigned i = 0; i < m_meshOverlap->vecSourceFaceIx.size(); i++ )
            {
                if ( m_meshOverlap->vecSourceFaceIx[i] + 1 > ixSourceFaceMax )
                {
                    ixSourceFaceMax = m_meshOverlap->vecSourceFaceIx[i] + 1;
                }
                if ( m_meshOverlap->vecTargetFaceIx[i] + 1 > ixTargetFaceMax )
                {
                    ixTargetFaceMax = m_meshOverlap->vecTargetFaceIx[i] + 1;
                }
            }

            // Check for forward correspondence in overlap mesh
            if ( m_meshInput->faces.size() - ixSourceFaceMax == 0 )
            {
                if ( is_root ) dbgprint.printf ( 0, "Overlap mesh forward correspondence found\n" );
            }
            else if ( m_meshOutput->faces.size() - ixSourceFaceMax == 0 )
            {   // Check for reverse correspondence in overlap mesh
                if ( is_root ) dbgprint.printf ( 0, "Overlap mesh reverse correspondence found (reversing)\n" );

                // Reorder overlap mesh
                m_meshOverlap->ExchangeFirstAndSecondMesh();
            }
            // else
            // {   // No correspondence found
            //     _EXCEPTION4 ( "Invalid overlap mesh:\n"
            //                   "    No correspondence found with input and output meshes (%i,%i) vs (%i,%i)",
            //                   m_meshInputCov->faces.size(), m_meshOutput->faces.size(), ixSourceFaceMax, ixTargetFaceMax );
            // }


            // Calculate Face areas
            if ( is_root ) dbgprint.printf ( 0, "Calculating overlap mesh Face areas\n" );
            double dTotalAreaOverlap_loc = m_meshOverlap->CalculateFaceAreas(false);
            double dTotalAreaOverlap = dTotalAreaOverlap_loc;
    #ifdef MOAB_HAVE_MPI
            if (m_pcomm) MPI_Reduce ( &dTotalAreaOverlap_loc, &dTotalAreaOverlap, 1, MPI_DOUBLE, MPI_SUM, 0, m_pcomm->comm() );
    #endif
            if ( is_root ) dbgprint.printf ( 0, "Overlap Mesh Area: %1.15e\n", dTotalAreaOverlap );

            // Partial cover
            if ( fabs ( dTotalAreaOverlap - dTotalAreaInput ) > 1.0e-10 )
            {
                if ( !fNoCheck )
                {
                    if ( is_root ) dbgprint.printf ( 0, "WARNING: Significant mismatch between overlap mesh area "
                                                                "and input mesh area.\n  Automatically enabling --nocheck\n" );
                    fNoCheck = true;
                }
            }

            /*
                // Recalculate input mesh area from overlap mesh
                if (fabs(dTotalAreaOverlap - dTotalAreaInput) > 1.0e-10) {
                    dbgprint.printf(0, "Overlap mesh only covers a sub-area of the sphere\n");
                    dbgprint.printf(0, "Recalculating source mesh areas\n");
                    dTotalAreaInput = m_meshInput->CalculateFaceAreasFromOverlap(m_meshOverlap);
                    dbgprint.printf(0, "New Input Mesh Geometric Area: %1.15e\n", dTotalAreaInput);
                }
            */
        }

        // Finite volume input / Finite volume output
        if ( ( eInputType  == DiscretizationType_FV ) &&
                ( eOutputType == DiscretizationType_FV )
           )
        {
            if(m_meshInputCov->faces.size()>0)
            {

              // Generate reverse node array and edge map
              m_meshInputCov->ConstructReverseNodeArray();
              m_meshInputCov->ConstructEdgeMap();

              // Initialize coordinates for map
              this->InitializeSourceCoordinatesFromMeshFV ( *m_meshInputCov );
              this->InitializeTargetCoordinatesFromMeshFV ( *m_meshOutput );

              // Finite volume input / Finite element output
              rval = this->SetDOFmapAssociation(eInputType, false, NULL, NULL, eOutputType, false, NULL);MB_CHK_ERR(rval);

              // Construct remap
              if ( is_root ) dbgprint.printf ( 0, "Calculating remap weights\n" );
              LinearRemapFVtoFV_Tempest_MOAB ( nPin );
            }
        }
        else if ( eInputType == DiscretizationType_FV )
        {
            DataArray3D<double> dataGLLJacobian;

            if ( is_root ) dbgprint.printf ( 0, "Generating output mesh meta data\n" );
            double dNumericalArea_loc =
                GenerateMetaData (
                    *m_meshOutput,
                    nPout,
                    fBubble,
                    dataGLLNodesDest,
                    dataGLLJacobian );

            double dNumericalArea = dNumericalArea_loc;
#ifdef MOAB_HAVE_MPI
            if (m_pcomm) MPI_Allreduce ( &dNumericalArea_loc, &dNumericalArea, 1, MPI_DOUBLE, MPI_SUM, m_pcomm->comm() );
#endif
            if ( is_root ) dbgprint.printf ( 0, "Output Mesh Numerical Area: %1.15e\n", dNumericalArea );

            // Initialize coordinates for map
            this->InitializeSourceCoordinatesFromMeshFV ( *m_meshInputCov );
            this->InitializeTargetCoordinatesFromMeshFE (
                *m_meshOutput, nPout, dataGLLNodesDest );

            // Generate the continuous Jacobian
            bool fContinuous = ( eOutputType == DiscretizationType_CGLL );

            if ( eOutputType == DiscretizationType_CGLL )
            {
                GenerateUniqueJacobian (
                    dataGLLNodesDest,
                    dataGLLJacobian,
                    this->GetTargetAreas() );
            }
            else
            {
                GenerateDiscontinuousJacobian (
                    dataGLLJacobian,
                    this->GetTargetAreas() );
            }

            // Generate reverse node array and edge map
            m_meshInputCov->ConstructReverseNodeArray();
            m_meshInputCov->ConstructEdgeMap();

            // Finite volume input / Finite element output
            rval = this->SetDOFmapAssociation(eInputType, false, NULL, NULL, 
                eOutputType, (eOutputType == DiscretizationType_CGLL), &dataGLLNodesDest);MB_CHK_ERR(rval);

            // Generate remap weights
            if ( is_root ) dbgprint.printf ( 0, "Calculating remap weights\n" );

            if ( fVolumetric )
            {
                LinearRemapFVtoGLL_Volumetric (
                    *m_meshInputCov,
                    *m_meshOutput,
                    *m_meshOverlap,
                    dataGLLNodesDest,
                    dataGLLJacobian,
                    this->GetTargetAreas(),
                    nPin,
                    *this,
                    nMonotoneType,
                    fContinuous,
                    fNoConservation );
            }
            else
            {
                LinearRemapFVtoGLL (
                    *m_meshInputCov,
                    *m_meshOutput,
                    *m_meshOverlap,
                    dataGLLNodesDest,
                    dataGLLJacobian,
                    this->GetTargetAreas(),
                    nPin,
                    *this,
                    nMonotoneType,
                    fContinuous,
                    fNoConservation );
            }
            
        }
        else if ( ( eInputType  == DiscretizationType_PCLOUD ) ||
                ( eOutputType == DiscretizationType_PCLOUD )
           )
        {
            DataArray3D<double> dataGLLJacobian;
            if (!m_bPointCloudSource)
            {
                // Generate reverse node array and edge map
                m_meshInputCov->ConstructReverseNodeArray();
                m_meshInputCov->ConstructEdgeMap();

                // Initialize coordinates for map
                if (eInputType == DiscretizationType_FV) {
                    this->InitializeSourceCoordinatesFromMeshFV ( *m_meshInputCov );
                }
                else
                {
                    if ( is_root ) dbgprint.printf ( 0, "Generating input mesh meta data\n" );
                    DataArray3D<double> dataGLLJacobianSrc;
                    GenerateMetaData (
                        *m_meshInputCov,
                        nPin,
                        fBubble,
                        dataGLLNodesSrcCov,
                        dataGLLJacobian );
                    GenerateMetaData (
                        *m_meshInput,
                        nPin,
                        fBubble,
                        dataGLLNodesSrc,
                        dataGLLJacobianSrc );
                }
            }
            else { /* Source is a point cloud dataset */

            }

            if (!m_bPointCloudTarget)
            {
                // Generate reverse node array and edge map
                m_meshOutput->ConstructReverseNodeArray();
                m_meshOutput->ConstructEdgeMap();

                // Initialize coordinates for map
                if (eOutputType == DiscretizationType_FV) {
                    this->InitializeSourceCoordinatesFromMeshFV ( *m_meshOutput );
                }
                else
                {
                    if ( is_root ) dbgprint.printf ( 0, "Generating output mesh meta data\n" );
                    GenerateMetaData (
                        *m_meshOutput,
                        nPout,
                        fBubble,
                        dataGLLNodesDest,
                        dataGLLJacobian );
                }
            }
            else { /* Target is a point cloud dataset */

            }

            // Finite volume input / Finite element output
            rval = this->SetDOFmapAssociation(eInputType, (eInputType == DiscretizationType_CGLL), 
                                                (m_bPointCloudSource || eInputType == DiscretizationType_FV ? NULL : &dataGLLNodesSrcCov), 
                                                (m_bPointCloudSource || eInputType == DiscretizationType_FV ? NULL : &dataGLLNodesSrc), 
                                                eOutputType, (eOutputType == DiscretizationType_CGLL), 
                                                (m_bPointCloudTarget ? NULL : &dataGLLNodesDest));MB_CHK_ERR(rval);

            // Construct remap
            if ( is_root ) dbgprint.printf ( 0, "Calculating remap weights with Nearest-Neighbor method\n" );
            rval = LinearRemapNN_MOAB(true /*use_GID_matching*/, false /*strict_check*/ );MB_CHK_ERR(rval);
        }
        else if (
            ( eInputType != DiscretizationType_FV ) &&
            ( eOutputType == DiscretizationType_FV )
        )
        {
            DataArray3D<double> dataGLLJacobianSrc, dataGLLJacobian;

            if ( is_root ) dbgprint.printf ( 0, "Generating input mesh meta data\n" );
            // double dNumericalAreaCov_loc =
                GenerateMetaData (
                    *m_meshInputCov,
                    nPin,
                    fBubble,
                    dataGLLNodesSrcCov,
                    dataGLLJacobian );

            double dNumericalArea_loc =
                GenerateMetaData (
                    *m_meshInput,
                    nPin,
                    fBubble,
                    dataGLLNodesSrc,
                    dataGLLJacobianSrc );

            // if ( is_root ) dbgprint.printf ( 0, "Input Mesh: Coverage Area: %1.15e, Output Area: %1.15e\n", dNumericalAreaCov_loc, dTotalAreaOutput_loc );
            // assert(dNumericalAreaCov_loc >= dTotalAreaOutput_loc);

            double dNumericalArea = dNumericalArea_loc;
#ifdef MOAB_HAVE_MPI
            if (m_pcomm) MPI_Allreduce ( &dNumericalArea_loc, &dNumericalArea, 1, MPI_DOUBLE, MPI_SUM, m_pcomm->comm() );
#endif
            if ( is_root ) dbgprint.printf ( 0, "Input Mesh Numerical Area: %1.15e\n", dNumericalArea );

            if ( fabs ( dNumericalArea - dTotalAreaInput ) > 1.0e-12 )
            {
                dbgprint.printf ( 0, "WARNING: Significant mismatch between input mesh "
                                  "numerical area and geometric area\n" );
            }

            if ( dataGLLNodesSrcCov.GetSubColumns() != m_meshInputCov->faces.size() )
            {
                _EXCEPTIONT ( "Number of element does not match between metadata and "
                              "input mesh" );
            }

            // Initialize coordinates for map
            this->InitializeSourceCoordinatesFromMeshFE (
                *m_meshInputCov, nPin, dataGLLNodesSrcCov );
            this->InitializeSourceCoordinatesFromMeshFE (
                *m_meshInput, nPin, dataGLLNodesSrc );
            this->InitializeTargetCoordinatesFromMeshFV ( *m_meshOutput );

            // Generate the continuous Jacobian for input mesh
            bool fContinuousIn = ( eInputType == DiscretizationType_CGLL );

            if ( eInputType == DiscretizationType_CGLL )
            {
                GenerateUniqueJacobian (
                    dataGLLNodesSrcCov,
                    dataGLLJacobian,
                    this->GetSourceAreas() );
            }
            else
            {
                GenerateDiscontinuousJacobian (
                    dataGLLJacobian,
                    this->GetSourceAreas() );
            }

            // Finite element input / Finite volume output
            rval = this->SetDOFmapAssociation(eInputType, (eInputType == DiscretizationType_CGLL), &dataGLLNodesSrcCov, &dataGLLNodesSrc, 
                eOutputType, false, NULL);MB_CHK_ERR(rval);

            // Generate remap
            if ( is_root ) dbgprint.printf ( 0, "Calculating remap weights\n" );

            if ( fVolumetric )
            {
                _EXCEPTIONT ( "Unimplemented: Volumetric currently unavailable for"
                              "GLL input mesh" );
            }

            LinearRemapSE4_Tempest_MOAB (
                dataGLLNodesSrcCov,
                dataGLLJacobian,
                nMonotoneType,
                fContinuousIn,
                fNoConservation
            );

        }
        else if (
            ( eInputType  != DiscretizationType_FV ) &&
            ( eOutputType != DiscretizationType_FV ) )
        {
            DataArray3D<double> dataGLLJacobianIn, dataGLLJacobianSrc;
            DataArray3D<double> dataGLLJacobianOut;

            // Input metadata
            if ( is_root ) dbgprint.printf ( 0, "Generating input mesh meta data\n" );
            double dNumericalAreaIn_loc =
                GenerateMetaData (
                    *m_meshInputCov,
                    nPin,
                    fBubble,
                    dataGLLNodesSrcCov,
                    dataGLLJacobianIn );

            double dNumericalAreaSrc_loc =
                GenerateMetaData (
                    *m_meshInput,
                    nPin,
                    fBubble,
                    dataGLLNodesSrc,
                    dataGLLJacobianSrc );

            assert(dNumericalAreaIn_loc >= dNumericalAreaSrc_loc);

            double dNumericalAreaIn = dNumericalAreaSrc_loc;
#ifdef MOAB_HAVE_MPI
            if (m_pcomm) MPI_Allreduce ( &dNumericalAreaSrc_loc, &dNumericalAreaIn, 1, MPI_DOUBLE, MPI_SUM, m_pcomm->comm() );
#endif
            if ( is_root ) dbgprint.printf ( 0, "Input Mesh Numerical Area: %1.15e\n", dNumericalAreaIn );

            if ( fabs ( dNumericalAreaIn - dTotalAreaInput ) > 1.0e-12 )
            {
                dbgprint.printf ( 0, "WARNING: Significant mismatch between input mesh "
                                  "numerical area and geometric area\n" );
            }

            // Output metadata
            if ( is_root ) dbgprint.printf ( 0, "Generating output mesh meta data\n" );
            double dNumericalAreaOut_loc =
                GenerateMetaData (
                    *m_meshOutput,
                    nPout,
                    fBubble,
                    dataGLLNodesDest,
                    dataGLLJacobianOut );

            double dNumericalAreaOut = dNumericalAreaOut_loc;
#ifdef MOAB_HAVE_MPI
            if (m_pcomm) MPI_Allreduce ( &dNumericalAreaOut_loc, &dNumericalAreaOut, 1, MPI_DOUBLE, MPI_SUM, m_pcomm->comm() );
#endif
            if ( is_root ) dbgprint.printf ( 0, "Output Mesh Numerical Area: %1.15e\n", dNumericalAreaOut );

            if ( fabs ( dNumericalAreaOut - dTotalAreaOutput ) > 1.0e-12 )
            {
                if ( is_root ) dbgprint.printf ( 0, "WARNING: Significant mismatch between output mesh "
                                                            "numerical area and geometric area\n" );
            }

            // Initialize coordinates for map
            this->InitializeSourceCoordinatesFromMeshFE (
                *m_meshInputCov, nPin, dataGLLNodesSrcCov );
            this->InitializeSourceCoordinatesFromMeshFE (
                *m_meshInput, nPin, dataGLLNodesSrc );
            this->InitializeTargetCoordinatesFromMeshFE (
                *m_meshOutput, nPout, dataGLLNodesDest );

            // Generate the continuous Jacobian for input mesh
            bool fContinuousIn = ( eInputType == DiscretizationType_CGLL );

            if ( eInputType == DiscretizationType_CGLL )
            {
                GenerateUniqueJacobian (
                    dataGLLNodesSrcCov,
                    dataGLLJacobianIn,
                    this->GetSourceAreas() );
            }
            else
            {
                GenerateDiscontinuousJacobian (
                    dataGLLJacobianIn,
                    this->GetSourceAreas() );
            }

            // Generate the continuous Jacobian for output mesh
            bool fContinuousOut = ( eOutputType == DiscretizationType_CGLL );

            if ( eOutputType == DiscretizationType_CGLL )
            {
                GenerateUniqueJacobian (
                    dataGLLNodesDest,
                    dataGLLJacobianOut,
                    this->GetTargetAreas() );
            }
            else
            {
                GenerateDiscontinuousJacobian (
                    dataGLLJacobianOut,
                    this->GetTargetAreas() );
            }

            // Input Finite Element to Output Finite Element
            rval = this->SetDOFmapAssociation(eInputType, (eInputType == DiscretizationType_CGLL), &dataGLLNodesSrcCov, &dataGLLNodesSrc, 
                eOutputType, (eOutputType == DiscretizationType_CGLL), &dataGLLNodesDest);MB_CHK_ERR(rval);

            // Generate remap
            if ( is_root ) dbgprint.printf ( 0, "Calculating remap weights\n" );

            LinearRemapGLLtoGLL2_MOAB (
                dataGLLNodesSrcCov,
                dataGLLJacobianIn,
                dataGLLNodesDest,
                dataGLLJacobianOut,
                this->GetTargetAreas(),
                nPin,
                nPout,
                nMonotoneType,
                fContinuousIn,
                fContinuousOut,
                fNoConservation
            );

        }
        else
        {
            _EXCEPTIONT ( "Not implemented" );
        }

#ifdef MOAB_HAVE_EIGEN
        copy_tempest_sparsemat_to_eigen3();
#endif

        // // Let us alos write out the TempestRemap equivalent so that we can do some verification checks
        // if ( is_root && size == 1)
        // {
        //     dbgprint.printf ( 0, "NOTE: Writing out moab_intersection mesh in TempestRemap format\n" );
        //     m_meshOverlap->Write ( "moab_intersection_tempest.g");
        // }

#ifdef MOAB_HAVE_MPI
        moab::Range sharedGhostEntities;
        rval = this->remove_ghosted_overlap_entities(sharedGhostEntities); MB_CHK_ERR ( rval );
        /* Use the following call to re-add them back */
        // moab::EntityHandle m_meshOverlapSet = m_remapper->GetMeshSet ( moab::Remapper::IntersectedMesh );
        // rval = m_interface->add_entities(m_meshOverlapSet, sharedGhostEntities);MB_CHK_SET_ERR(rval, "Adding entities dim 2 failed");
#endif

        // Verify consistency, conservation and monotonicity, globally
#ifdef MOAB_HAVE_MPI
        // first, we have to agree if checks are needed globally
        // if there is at least one that does not want checks, no-one should do checks
        int fck_int_loc = fNoCheck ? 1 : 0;
        int fck_int_glob = fck_int_loc;
        if (m_pcomm) MPI_Allreduce ( &fck_int_loc, &fck_int_glob, 1, MPI_INT, MPI_MAX, m_pcomm->comm() );
        fNoCheck = (0==fck_int_glob)? false : true;
#endif
        if ( !fNoCheck )
        {
            if ( is_root ) dbgprint.printf ( 0, "Verifying map" );
            this->IsConsistent ( 1.0e-8 );
            if ( !fNoConservation ) this->IsConservative ( 1.0e-8 );

            if ( nMonotoneType != 0 )
            {
                this->IsMonotone ( 1.0e-12 );
            }
        }

        // Apply Remapping Weights to data
        if ( strInputData != "" )
        {
            if ( is_root ) dbgprint.printf ( 0, "Applying remap weights to data\n" );

            this->SetFillValueOverride ( static_cast<float> ( dFillValueOverride ) );
            this->Apply (
                strInputData,
                strOutputData,
                vecVariableStrings,
                strNColName,
                fOutputDouble,
                false );
        }

        // Copy variables from input file to output file
        if ( ( strInputData != "" ) && ( strOutputData != "" ) )
        {
            if ( fPreserveAll )
            {
                if ( is_root ) dbgprint.printf ( 0, "Preserving variables" );
                this->PreserveAllVariables ( strInputData, strOutputData );

            }
            else if ( vecPreserveVariableStrings.size() != 0 )
            {
                if ( is_root ) dbgprint.printf ( 0, "Preserving variables" );
                this->PreserveVariables (
                    strInputData,
                    strOutputData,
                    vecPreserveVariableStrings );
            }
        }

    }
    catch ( Exception & e )
    {
        dbgprint.printf ( 0, "%s", e.ToString().c_str() );
        return ( moab::MB_FAILURE );

    }
    catch ( ... )
    {
        return ( moab::MB_FAILURE );
    }
    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////

int moab::TempestOnlineMap::IsConsistent (double dTolerance)
{
#ifndef MOAB_HAVE_MPI
    
    return OfflineMap::IsConsistent(dTolerance);

#else

    // Get map entries
    DataArray1D<int> dataRows;
    DataArray1D<int> dataCols;
    DataArray1D<double> dataEntries;

    // Calculate row sums
    DataArray1D<double> dRowSums;
    m_mapRemap.GetEntries ( dataRows, dataCols, dataEntries );
    dRowSums.Allocate ( m_mapRemap.GetRows() );

    for ( unsigned i = 0; i < dataRows.GetRows(); i++ )
    {
        dRowSums[dataRows[i]] += dataEntries[i];
    }

    // Verify all row sums are equal to 1
    int fConsistent = 0;
    for ( unsigned i = 0; i < dRowSums.GetRows(); i++ )
    {
        if ( fabs ( dRowSums[i] - 1.0 ) > dTolerance )
        {
            fConsistent++;
            Announce ( "TempestOnlineMap is not consistent in row %i (%1.15e)",
                       i, dRowSums[i] );
        }
    }

    int ierr;
    int fConsistentGlobal = 0;
    ierr = MPI_Allreduce(&fConsistent, &fConsistentGlobal, 1, MPI_INT, MPI_SUM, m_pcomm->comm());
    if ( ierr != MPI_SUCCESS ) return -1;

    return fConsistentGlobal;
#endif
}

///////////////////////////////////////////////////////////////////////////////

int moab::TempestOnlineMap::IsConservative (double dTolerance)
{
#ifndef MOAB_HAVE_MPI

    return OfflineMap::IsConservative(dTolerance);

#else
    // return OfflineMap::IsConservative(dTolerance);

    int ierr;
    // Get map entries
    DataArray1D<int> dataRows;
    DataArray1D<int> dataCols;
    DataArray1D<double> dataEntries;
    const DataArray1D<double>& dTargetAreas = this->GetTargetAreas();
    const DataArray1D<double>& dSourceAreas = this->GetSourceAreas();

    // Calculate column sums
    std::vector<int> dColumnsUnique;
    std::vector<double> dColumnSums;

    int nColumns = m_mapRemap.GetColumns();
    m_mapRemap.GetEntries ( dataRows, dataCols, dataEntries );
    dColumnSums.resize ( m_nTotDofs_SrcCov, 0.0 );
    dColumnsUnique.resize ( m_nTotDofs_SrcCov, -1 );

    for ( unsigned i = 0; i < dataEntries.GetRows(); i++ )
    {
        dColumnSums[dataCols[i]] += dataEntries[i] * dTargetAreas[dataRows[i]];//  / dSourceAreas[dataCols[i]];

        assert(dataCols[i] < m_nTotDofs_SrcCov);

        // GID for column DoFs: col_gdofmap[ col_ldofmap [ dataCols[i] ] ]
        int colGID = col_gdofmap[ col_ldofmap [ dataCols[i] ] ];
        // int colGID = col_dofmap [ col_ldofmap [ dataCols[i] ] ];
        dColumnsUnique[ dataCols[i] ] = colGID;

        // std::cout << "Column dataCols[i]=" << dataCols[i] << " with GID = " << colGID << std::endl;
    }

    int rootProc = 0;
    std::vector<int> nElementsInProc;
    const int nDATA = 3;
    if (!rank) nElementsInProc.resize(size*nDATA);
    int senddata[nDATA] = {nColumns, m_nTotDofs_SrcCov, m_nTotDofs_Src};
    ierr = MPI_Gather ( senddata, nDATA, MPI_INTEGER, nElementsInProc.data(), nDATA, MPI_INTEGER, rootProc, m_pcomm->comm() );
    if ( ierr != MPI_SUCCESS ) return -1;

    int nTotVals = 0, nTotColumns = 0, nTotColumnsUnq = 0;
    std::vector<int> dColumnIndices;
    std::vector<double> dColumnSourceAreas;
    std::vector<double> dColumnSumsTotal;
    std::vector<int> displs, rcount;
    if (rank == rootProc)
    {
        displs.resize ( size+1, 0 );
        rcount.resize ( size, 0 );
        int gsum = 0;
        for (int ir = 0; ir < size; ++ir)
        {
            nTotVals += nElementsInProc[ir*nDATA];
            nTotColumns += nElementsInProc[ir*nDATA+1];
            nTotColumnsUnq += nElementsInProc[ir*nDATA+2];

            displs[ir] = gsum;
            rcount[ir] = nElementsInProc[ir*nDATA+1];
            gsum += rcount[ir];
        }

        dColumnIndices.resize ( nTotColumns, -1 );
        dColumnSumsTotal.resize ( nTotColumns, 0.0 );
        dColumnSourceAreas.resize ( nTotColumns, 0.0 );

        // std::cout << "Total columns: " << nTotVals << " cols: " << nTotColumns << " and unique ones: " << nTotColumnsUnq << std::endl;
    }

    // Gather all ColumnSums to root process and accumulate
    // We expect that the sums of all columns equate to 1.0 within user specified tolerance
    // Need to do a gatherv here since different processes have different number of elements
    // MPI_Reduce(&dColumnSums[0], &dColumnSumsTotal[0], m_mapRemap.GetColumns(), MPI_DOUBLE, MPI_SUM, 0, m_pcomm->comm());
    ierr = MPI_Gatherv ( &dColumnsUnique[0], m_nTotDofs_SrcCov, MPI_INTEGER, &dColumnIndices[0], rcount.data(), displs.data(), MPI_INTEGER, rootProc, m_pcomm->comm() );
    if ( ierr != MPI_SUCCESS ) return -1;
    ierr = MPI_Gatherv ( &dColumnSums[0], m_nTotDofs_SrcCov, MPI_DOUBLE, &dColumnSumsTotal[0], rcount.data(), displs.data(), MPI_DOUBLE, rootProc, m_pcomm->comm() );
    if ( ierr != MPI_SUCCESS ) return -1;
    ierr = MPI_Gatherv ( &dSourceAreas[0], m_nTotDofs_SrcCov, MPI_DOUBLE, &dColumnSourceAreas[0], rcount.data(), displs.data(), MPI_DOUBLE, rootProc, m_pcomm->comm() );
    if ( ierr != MPI_SUCCESS ) return -1;

    // Clean out unwanted arrays now
    dColumnSums.clear();
    dColumnsUnique.clear();

    // Verify all column sums equal the input Jacobian
    int fConservative = 0;
    if (rank == rootProc)
    {
        displs[size] = (nTotColumns);
        // std::vector<double> dColumnSumsOnRoot(nTotColumnsUnq, 0.0);
        std::map<int, double> dColumnSumsOnRoot;
        std::map<int, double> dColumnSourceAreasOnRoot;
        for ( int ir = 0; ir < size; ir++ )
        {
            for ( int ips = displs[ir]; ips < displs[ir+1]; ips++ )
            {
                assert(dColumnIndices[ips] < nTotColumnsUnq);
                dColumnSumsOnRoot[ dColumnIndices[ips] ] += dColumnSumsTotal[ips];// / dColumnSourceAreas[ips];
                dColumnSourceAreasOnRoot[ dColumnIndices[ips] ] += dColumnSourceAreas[ips];
                // dColumnSourceAreas[ dColumnIndices[ips] ]
            }
        }

        for(std::map<int, double>::iterator it = dColumnSumsOnRoot.begin(); it != dColumnSumsOnRoot.end(); ++it)
        {
            if ( fabs ( it->second - dColumnSourceAreasOnRoot[it->first] ) > dTolerance )
            // if ( fabs ( it->second - 1.0 ) > dTolerance )
            {
                fConservative++;
                Announce ( "TempestOnlineMap is not conservative in column "
                        // "%i (%1.15e)", it->first, it->second );
                        "%i (%1.15e)", it->first, it->second / dColumnSourceAreasOnRoot[it->first] );
            }
        }
    }

    // TODO: Just do a broadcast from root instead of a reduction
    ierr = MPI_Bcast(&fConservative, 1, MPI_INT, rootProc, m_pcomm->comm());
    if ( ierr != MPI_SUCCESS ) return -1;

    return fConservative;
#endif
}

///////////////////////////////////////////////////////////////////////////////

int moab::TempestOnlineMap::IsMonotone (double dTolerance)
{
#ifndef MOAB_HAVE_MPI

    return OfflineMap::IsMonotone(dTolerance);

#else

    // Get map entries
    DataArray1D<int> dataRows;
    DataArray1D<int> dataCols;
    DataArray1D<double> dataEntries;

    m_mapRemap.GetEntries ( dataRows, dataCols, dataEntries );

    // Verify all entries are in the range [0,1]
    int fMonotone = 0;
    for ( unsigned i = 0; i < dataRows.GetRows(); i++ )
    {
        if ( ( dataEntries[i] < -dTolerance ) ||
                ( dataEntries[i] > 1.0 + dTolerance )
           )
        {
            fMonotone++;

            Announce ( "TempestOnlineMap is not monotone in entry (%i): %1.15e",
                       i, dataEntries[i] );
        }
    }

    int ierr;
    int fMonotoneGlobal = 0;
    ierr = MPI_Allreduce(&fMonotone, &fMonotoneGlobal, 1, MPI_INT, MPI_SUM, m_pcomm->comm());
    if ( ierr != MPI_SUCCESS ) return -1;

    return fMonotoneGlobal;
#endif
}

///////////////////////////////////////////////////////////////////////////////

#ifdef MOAB_HAVE_EIGEN
void moab::TempestOnlineMap::InitVectors()
{
    //assert(m_weightMatrix.rows() != 0 && m_weightMatrix.cols() != 0);
    m_rowVector.resize( m_weightMatrix.rows() );
    m_colVector.resize( m_weightMatrix.cols() );
}
#endif

///////////////////////////////////////////////////////////////////////////////

moab::ErrorCode moab::TempestOnlineMap::remove_ghosted_overlap_entities (moab::Range& sharedGhostEntities)
{
    sharedGhostEntities.clear();
#ifdef MOAB_HAVE_MPI
    moab::ErrorCode rval;

    // Remove entities in the intersection mesh that are part of the ghosted overlap
    if (is_parallel && size > 1)
    {
        moab::Range allents;
        moab::EntityHandle m_meshOverlapSet = m_remapper->GetMeshSet ( moab::Remapper::IntersectedMesh );
        rval = m_interface->get_entities_by_dimension(m_meshOverlapSet, 2, allents);MB_CHK_SET_ERR(rval, "Getting entities dim 2 failed");

        moab::Range sharedents;
        moab::Tag ghostTag;
        std::vector<int> ghFlags(allents.size());
        rval = m_interface->tag_get_handle ( "ORIG_PROC", ghostTag ); MB_CHK_ERR ( rval );
        rval = m_interface->tag_get_data ( ghostTag,  allents, &ghFlags[0] ); MB_CHK_ERR ( rval );
        for (unsigned i=0; i < allents.size(); ++i)
            if (ghFlags[i]>=0) // it means it is a ghost overlap element
                sharedents.insert(allents[i]); // this should not participate in smat!

        allents = subtract(allents,sharedents);

        // Get connectivity from all ghosted elements and filter out
        // the vertices that are not owned
        moab::Range ownedverts, sharedverts;
        rval = m_interface->get_connectivity(allents, ownedverts);MB_CHK_SET_ERR(rval, "Deleting entities dim 0 failed");
        rval = m_interface->get_connectivity(sharedents, sharedverts);MB_CHK_SET_ERR(rval, "Deleting entities dim 0 failed");
        sharedverts = subtract(sharedverts,ownedverts);                    
        rval = m_interface->remove_entities(m_meshOverlapSet, sharedents);MB_CHK_SET_ERR(rval, "Deleting entities dim 2 failed");
        rval = m_interface->remove_entities(m_meshOverlapSet, sharedverts);MB_CHK_SET_ERR(rval, "Deleting entities dim 0 failed");

        sharedGhostEntities.merge(sharedents);
        sharedGhostEntities.merge(sharedverts);
    }
#endif
    return moab::MB_SUCCESS;
}

moab::ErrorCode moab::TempestOnlineMap::WriteParallelMap (std::string strOutputFile)  const
{
    moab::ErrorCode rval;

    moab::EntityHandle& m_meshOverlapSet = m_remapper->m_overlap_set;
    int tot_src_ents = m_remapper->m_source_entities.size();
    int tot_tgt_ents = m_remapper->m_target_entities.size();

    const int weightMatNNZ = m_weightMatrix.nonZeros();
    moab::Tag tagMapMetaData, tagMapIndexRow, tagMapIndexCol, tagMapValues, srcEleIDs, tgtEleIDs, srcAreaValues, tgtAreaValues, srcMaskValues, tgtMaskValues;
    rval = m_interface->tag_get_handle("SMAT_DATA", 13, moab::MB_TYPE_INTEGER, tagMapMetaData, moab::MB_TAG_CREAT|moab::MB_TAG_SPARSE);MB_CHK_SET_ERR(rval, "Retrieving tag handles failed");
    rval = m_interface->tag_get_handle("SMAT_ROWS", weightMatNNZ, moab::MB_TYPE_INTEGER, tagMapIndexRow, moab::MB_TAG_CREAT|moab::MB_TAG_SPARSE|moab::MB_TAG_VARLEN);MB_CHK_SET_ERR(rval, "Retrieving tag handles failed");
    rval = m_interface->tag_get_handle("SMAT_COLS", weightMatNNZ, moab::MB_TYPE_INTEGER, tagMapIndexCol, moab::MB_TAG_CREAT|moab::MB_TAG_SPARSE|moab::MB_TAG_VARLEN);MB_CHK_SET_ERR(rval, "Retrieving tag handles failed");
    rval = m_interface->tag_get_handle("SMAT_VALS", weightMatNNZ, moab::MB_TYPE_DOUBLE, tagMapValues, moab::MB_TAG_CREAT|moab::MB_TAG_SPARSE|moab::MB_TAG_VARLEN);MB_CHK_SET_ERR(rval, "Retrieving tag handles failed");
    rval = m_interface->tag_get_handle("SourceGIDS", this->m_dSourceAreas.GetRows(), moab::MB_TYPE_INTEGER, srcEleIDs, moab::MB_TAG_CREAT|moab::MB_TAG_SPARSE|moab::MB_TAG_VARLEN);MB_CHK_SET_ERR(rval, "Retrieving tag handles failed");
    rval = m_interface->tag_get_handle("TargetGIDS", this->m_dTargetAreas.GetRows(), moab::MB_TYPE_INTEGER, tgtEleIDs, moab::MB_TAG_CREAT|moab::MB_TAG_SPARSE|moab::MB_TAG_VARLEN);MB_CHK_SET_ERR(rval, "Retrieving tag handles failed");
    rval = m_interface->tag_get_handle("SourceAreas", this->m_dSourceAreas.GetRows(), moab::MB_TYPE_DOUBLE, srcAreaValues, moab::MB_TAG_CREAT|moab::MB_TAG_SPARSE|moab::MB_TAG_VARLEN);MB_CHK_SET_ERR(rval, "Retrieving tag handles failed");
    rval = m_interface->tag_get_handle("TargetAreas", this->m_dTargetAreas.GetRows(), moab::MB_TYPE_DOUBLE, tgtAreaValues, moab::MB_TAG_CREAT|moab::MB_TAG_SPARSE|moab::MB_TAG_VARLEN);MB_CHK_SET_ERR(rval, "Retrieving tag handles failed");
    if (m_iSourceMask.IsAttached()) {
        rval = m_interface->tag_get_handle("SourceMask", m_iSourceMask.GetRows(), moab::MB_TYPE_INTEGER, srcMaskValues, moab::MB_TAG_CREAT|moab::MB_TAG_SPARSE|moab::MB_TAG_VARLEN);MB_CHK_SET_ERR(rval, "Retrieving tag handles failed");
    }
    if (m_iTargetMask.IsAttached()) {
        rval = m_interface->tag_get_handle("TargetMask", m_iTargetMask.GetRows(), moab::MB_TYPE_INTEGER, tgtMaskValues, moab::MB_TAG_CREAT|moab::MB_TAG_SPARSE|moab::MB_TAG_VARLEN);MB_CHK_SET_ERR(rval, "Retrieving tag handles failed");
    }

    std::vector<int> smatrowvals(weightMatNNZ),smatcolvals(weightMatNNZ);
    const double* smatvals = m_weightMatrix.valuePtr();
    int maxrow=0, maxcol=0, offset=0;

    // Loop over the matrix entries and find the max global ID for rows and columns
    for (int k=0; k < m_weightMatrix.outerSize(); ++k)
    {
        for (moab::TempestOnlineMap::WeightMatrix::InnerIterator it(m_weightMatrix,k); it; ++it)
        {
            smatrowvals[offset] = row_gdofmap [ row_ldofmap [ it.row() ] ]; // this->GetRowGlobalDoF ( it.row() );
            smatcolvals[offset] = col_gdofmap [ col_ldofmap [ it.col() ] ]; // this->GetColGlobalDoF ( it.col() );
            maxrow = (smatrowvals[offset] > maxrow) ? smatrowvals[offset] : maxrow;
            maxcol = (smatcolvals[offset] > maxcol) ? smatcolvals[offset] : maxcol;
            ++offset;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // The metadata in H5M file contains the following data:
    //
    //   1. n_a: Total source entities: (number of elements in source mesh)
    //   2. n_b: Total target entities: (number of elements in target mesh)
    //   3. nv_a: Max edge size of elements in source mesh
    //   4. nv_b: Max edge size of elements in target mesh
    //   5. maxrows: Number of rows in remap weight matrix
    //   6. maxcols: Number of cols in remap weight matrix
    //   7. nnz: Number of total nnz in sparse remap weight matrix
    //   8. np_a: The order of the field description on the source mesh: >= 1 
    //   9. np_b: The order of the field description on the target mesh: >= 1
    //   10. method_a: The type of discretization for field on source mesh: [0 = FV, 1 = cGLL, 2 = dGLL]
    //   11. method_b: The type of discretization for field on target mesh: [0 = FV, 1 = cGLL, 2 = dGLL]
    //   12. conserved: Flag to specify whether the remap operator has conservation constraints: [0, 1]
    //   13. monotonicity: Flags to specify whether the remap operator has monotonicity constraints: [0, 1, 2]
    //
    ///////////////////////////////////////////////////////////////////////////
    int map_disc_details[6];
    map_disc_details[0] = m_nDofsPEl_Src;
    map_disc_details[1] = m_nDofsPEl_Dest;
    map_disc_details[2] = (m_srcDiscType == DiscretizationType_FV ? 0 : (m_srcDiscType == DiscretizationType_CGLL ? 1 : 2));
    map_disc_details[3] = (m_destDiscType == DiscretizationType_FV ? 0 : (m_destDiscType == DiscretizationType_CGLL ? 1 : 2));
    map_disc_details[4] = (m_bConserved ? 1 : 0);
    map_disc_details[5] = m_iMonotonicity;

#ifdef MOAB_HAVE_MPI
    int loc_smatmetadata[13] = {tot_src_ents, tot_tgt_ents, m_remapper->max_source_edges, m_remapper->max_target_edges, maxrow+1, maxcol+1, weightMatNNZ,
                                map_disc_details[0], map_disc_details[1], map_disc_details[2], map_disc_details[3], map_disc_details[4], map_disc_details[5]};
    rval = m_interface->tag_set_data(tagMapMetaData, &m_meshOverlapSet, 1, &loc_smatmetadata[0]);MB_CHK_SET_ERR(rval, "Setting local tag data failed");
    int glb_smatmetadata[13] = {0, 0, 0, 0, 0, 0, 0,
                                map_disc_details[0], map_disc_details[1], map_disc_details[2], map_disc_details[3], map_disc_details[4], map_disc_details[5]};
    int loc_buf[7] = {tot_src_ents, tot_tgt_ents, weightMatNNZ, m_remapper->max_source_edges, m_remapper->max_target_edges, maxrow, maxcol};
    int glb_buf[4] = {0,0,0,0};
    MPI_Reduce(&loc_buf[0], &glb_buf[0], 3, MPI_INT, MPI_SUM, 0, m_pcomm->comm());
    glb_smatmetadata[0] = glb_buf[0];
    glb_smatmetadata[1] = glb_buf[1];
    glb_smatmetadata[6] = glb_buf[2];
    MPI_Reduce(&loc_buf[3], &glb_buf[0], 4, MPI_INT, MPI_MAX, 0, m_pcomm->comm());
    glb_smatmetadata[2] = glb_buf[0];
    glb_smatmetadata[3] = glb_buf[1];
    glb_smatmetadata[4] = glb_buf[2];
    glb_smatmetadata[5] = glb_buf[3];
#else
    int glb_smatmetadata[13] = {tot_src_ents, tot_tgt_ents, m_remapper->max_source_edges, m_remapper->max_target_edges, maxrow, maxcol, weightMatNNZ,
                                map_disc_details[0], map_disc_details[1], map_disc_details[2], map_disc_details[3], map_disc_details[4], map_disc_details[5]};
#endif
    // These values represent number of rows and columns. So should be 1-based.
    glb_smatmetadata[4]++;
    glb_smatmetadata[5]++;

    if (this->is_root)
    {
        std::cout << "Global data = " << glb_smatmetadata[0] << " " << glb_smatmetadata[1] << " " << glb_smatmetadata[2] << " " << glb_smatmetadata[3] << " " << glb_smatmetadata[4] << " " << glb_smatmetadata[5] << " " << glb_smatmetadata[6] << "\n"; 
        std::cout << "  " << this->rank << "  Writing remap weights with size [" << glb_smatmetadata[4] << " X " << glb_smatmetadata[5] << "] and NNZ = " << glb_smatmetadata[6] << std::endl;
        EntityHandle root_set = 0;
        rval = m_interface->tag_set_data(tagMapMetaData, &root_set, 1, &glb_smatmetadata[0]);MB_CHK_SET_ERR(rval, "Setting local tag data failed");
    }

    int dsize;
    const int numval = weightMatNNZ;
    const void* smatrowvals_d = smatrowvals.data();
    const void* smatcolvals_d = smatcolvals.data();
    const void* smatvals_d = smatvals;
    rval = m_interface->tag_set_by_ptr(tagMapIndexRow, &m_meshOverlapSet, 1, &smatrowvals_d, &numval);MB_CHK_SET_ERR(rval, "Setting local tag data failed");
    rval = m_interface->tag_set_by_ptr(tagMapIndexCol, &m_meshOverlapSet, 1, &smatcolvals_d, &numval);MB_CHK_SET_ERR(rval, "Setting local tag data failed");
    rval = m_interface->tag_set_by_ptr(tagMapValues, &m_meshOverlapSet, 1, &smatvals_d, &numval);MB_CHK_SET_ERR(rval, "Setting local tag data failed");

    /* Set the global IDs for the DoFs */
    ////
    // col_gdofmap [ col_ldofmap [ 0 : local_ndofs ] ] = GDOF
    // row_gdofmap [ row_ldofmap [ 0 : local_ndofs ] ] = GDOF
    ////
    std::vector<int> src_global_dofs(m_dSourceAreas.GetRows()), tgt_global_dofs(m_dTargetAreas.GetRows());
    for (unsigned i=0; i < m_dSourceAreas.GetRows(); ++i)
        src_global_dofs[i] = col_gdofmap [ col_ldofmap [ i ] ];
    for (unsigned i=0; i < m_dTargetAreas.GetRows(); ++i)
        tgt_global_dofs[i] = row_gdofmap [ row_ldofmap [ i ] ];
    const void* srceleidvals_d = src_global_dofs.data(); //this->col_gdofmap.data();
    const void* tgteleidvals_d = tgt_global_dofs.data(); //this->row_gdofmap.data();
    dsize = src_global_dofs.size();
    rval = m_interface->tag_set_by_ptr(srcEleIDs, &m_meshOverlapSet, 1, &srceleidvals_d, &dsize);MB_CHK_SET_ERR(rval, "Setting local tag data failed");
    dsize = tgt_global_dofs.size();
    rval = m_interface->tag_set_by_ptr(tgtEleIDs, &m_meshOverlapSet, 1, &tgteleidvals_d, &dsize);MB_CHK_SET_ERR(rval, "Setting local tag data failed");

    /* Set the source and target areas */
    const void* srcareavals_d = /*m_remapper->m_source->vecFaceArea*/ m_dSourceAreas;
    const void* tgtareavals_d = /*m_remapper->m_target->vecFaceArea*/ m_dTargetAreas;
    dsize = /*m_remapper->m_source->vecFaceArea.GetRows()*/m_dSourceAreas.GetRows();
    rval = m_interface->tag_set_by_ptr(srcAreaValues, &m_meshOverlapSet, 1, &srcareavals_d, &dsize);MB_CHK_SET_ERR(rval, "Setting local tag data failed");
    dsize = /*m_remapper->m_target->vecFaceArea.GetRows()*/m_dTargetAreas.GetRows();
    rval = m_interface->tag_set_by_ptr(tgtAreaValues, &m_meshOverlapSet, 1, &tgtareavals_d, &dsize);MB_CHK_SET_ERR(rval, "Setting local tag data failed");

    if (m_iSourceMask.IsAttached()) {
        const void* srcmaskvals_d = m_iSourceMask;
        dsize = m_iSourceMask.GetRows();
        rval = m_interface->tag_set_by_ptr(srcMaskValues, &m_meshOverlapSet, 1, &srcmaskvals_d, &dsize);MB_CHK_SET_ERR(rval, "Setting local tag data failed");
    }

    if (m_iTargetMask.IsAttached()) {
        const void* tgtmaskvals_d = m_iTargetMask;
        dsize = m_iTargetMask.GetRows();
        rval = m_interface->tag_set_by_ptr(tgtMaskValues, &m_meshOverlapSet, 1, &tgtmaskvals_d, &dsize);MB_CHK_SET_ERR(rval, "Setting local tag data failed");
    }


#ifdef MOAB_HAVE_MPI
    const char *writeOptions = (this->size > 1 ? "PARALLEL=WRITE_PART" : "");
#else
    const char *writeOptions = "";
#endif

    // EntityHandle sets[3] = {m_remapper->m_source_set, m_remapper->m_target_set, m_remapper->m_overlap_set};
    // rval = m_interface->write_file ( strOutputFile.c_str(), NULL, writeOptions, sets, 3 ); MB_CHK_ERR ( rval );
    rval = m_interface->write_file ( strOutputFile.c_str(), NULL, writeOptions ); MB_CHK_ERR ( rval );

#ifdef WRITE_SCRIP_FILE
    sstr.str("");
    sstr << ctx.outFilename.substr(0, lastindex) << "_" << proc_id << ".nc";
    std::map<std::string, std::string> mapAttributes;
    mapAttributes["Creator"] = "MOAB mbtempest workflow";
    if (!ctx.proc_id) std::cout << "Writing offline map to file: " << sstr.str() << std::endl;
    this->Write(strOutputFile.c_str(), mapAttributes, NcFile::Netcdf4);
    sstr.str("");
#endif

    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
