///////////////////////////////////////////////////////////////////////////////
///
///	\file    OfflineMap.h
///	\author  Paul Ullrich
///	\version August 14, 2014
///
///	<remarks>
///		Copyright 2000-2014 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _TEMPESTOFFLINEMAP_H_
#define _TEMPESTOFFLINEMAP_H_

#include "SparseMatrix.h"
#include "DataVector.h"
#include "DataMatrix.h"
#include "DataMatrix3D.h"
#include "OfflineMap.h"
#include <string>
#include <vector>

#include "moab/TempestRemapper.hpp"

class Mesh;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An offline map between two Meshes.
///	</summary>
class TempestOfflineMap : public OfflineMap {

public:

	///	<summary>
	///		Generate the metadata associated with the offline map.
	///	</summary>
	TempestOfflineMap(moab::TempestRemapper* remapper);

	///	<summary>
	///		Define a virtual destructor.
	///	</summary>
	virtual ~TempestOfflineMap();

	///	<summary>
	///		Gather the mapping matrix that was computed in different processors and accumulate the data
	///     on the root so that OfflineMap can be generated in parallel.
	///	</summary>
	virtual void GatherAllToRoot();

public:
	///	<summary>
	///		Generate the offline map, given the source and target mesh and discretization details.
	///     This method generates the mapping between the two meshes based on the overlap and stores 
	///     the result in the SparseMatrix.
	///	</summary>
	moab::ErrorCode GenerateOfflineMap( std::string strInputType, std::string strOutputType,
                                        int nPin=4, int nPout=4,
                                        bool fBubble=false, int fMonotoneTypeID=0,
                                        bool fVolumetric=false, bool fNoConservation=false, bool fNoCheck=false,
                                        std::string strVariables="", std::string strOutputMap="",
                                        std::string strInputData="", std::string strOutputData="",
                                        std::string strNColName="", bool fOutputDouble=false,
                                        std::string strPreserveVariables="", bool fPreserveAll=false, double dFillValueOverride=0.0 );

	///	<summary>
	///		Generate the metadata associated with the offline map.
	///	</summary>
	// moab::ErrorCode GenerateMetaData();

public:

	///	<summary>
	///		Read the OfflineMap from a NetCDF file.
	///	</summary>
	// virtual void Read(
	// 	const std::string & strSource
	// );

	///	<summary>
	///		Write the TempestOfflineMap to a parallel NetCDF file.
	///	</summary>
	virtual void Write(
		const std::string & strTarget
	);

public:
	///	<summary>
	///		Determine if the map is first-order accurate.
	///	</summary>
	virtual bool IsConsistent(
		double dTolerance
	);

	///	<summary>
	///		Determine if the map is conservative.
	///	</summary>
	virtual bool IsConservative(
		double dTolerance
	);

	///	<summary>
	///		Determine if the map is monotone.
	///	</summary>
	virtual bool IsMonotone(
		double dTolerance
	);

private:
	///	<summary>
	///		Compute the remapping weights for a FV field defined on the source to a 
	///     FV field defined on the target mesh.
	///	</summary>
	void LinearRemapFVtoFV_Tempest_MOAB( int nOrder );

protected:
	///	<summary>
	///		The fundamental remapping operator object.
	///	</summary>
	moab::TempestRemapper* m_remapper;

	///	<summary>
	///		The SparseMatrix representing this operator.
	///	</summary>
	SparseMatrix<double> m_mapRemapGlobal;

	///	<summary>
	///		The reference to the moab::Core object that contains source/target and overlap sets.
	///	</summary>
	moab::Interface* mbCore;

	///	<summary>
	///		The reference to the parallel communicator object used by the Core object.
	///	</summary>
	moab::ParallelComm* pcomm;

	Mesh* m_meshInput;
	Mesh* m_meshInputCov;
	Mesh* m_meshOutput;
	Mesh* m_meshOverlap;

};

///////////////////////////////////////////////////////////////////////////////

#endif

