// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Common/CFMap.hh"
#include "Common/PE.hh"

#include "Petsc/PetscIdxMapping.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

void PetscIdxMapping::computeMapping(std::valarray<CFuint>& localToPetscIDs,
				     std::valarray<bool>& isNonLocalRow)
{
  CFMap<CFuint, CFuint> globalToPetscIDs;

  cf_assert(_globalIDs.size() == _isGhost.size());
  const CFuint nbPoints = _globalIDs.size();
  globalToPetscIDs.reserve(nbPoints);

  // allocation and initialization of the flags telling if a row is non locally
  // owned
  isNonLocalRow.resize(nbPoints);
  isNonLocalRow = false;

  const CFuint currRank = Common::PE::GetPE().GetRank();
  const CFuint nbProcesses = Common::PE::GetPE().GetProcessorCount();
  cf_assert(currRank < nbProcesses);
  cf_assert(_updatablesPerRank.size() == nbProcesses);

  // compute the start global
  std::valarray<CFuint> startPetscID(nbProcesses);
  startPetscID = 0;

  for (CFuint ip = 0; ip < nbProcesses; ++ip) {
    for (CFuint ir = 0; ir < ip; ++ir) {
      startPetscID[ip] += _updatablesPerRank[ir].size();
    }
    CFLog(NOTICE, "startPetscID[" << ip << "] = " << startPetscID[ip] << "\n");
  }

  vector<CFuint> globalGhostIDs;
  // at first compute and store the Petsc IDs for locally owned points
  // store global IDs of ghost points
  CFuint startGlobalID = startPetscID[currRank];
  for (CFuint i = 0; i < nbPoints; ++i) {
    if (!_isGhost[i]) {
      globalToPetscIDs.insert(_globalIDs[i],startGlobalID);
      ++startGlobalID;
    }
    else {
      globalGhostIDs.push_back(_globalIDs[i]);
      isNonLocalRow[i] = true;
    }
  }

  // process global IDs of ghost states (owned by other processes)
  // to calculate the corresponding global Petsc IDs
  for (CFuint i = 0; i < globalGhostIDs.size(); ++i) {
    const CFuint globGhostID = globalGhostIDs[i];
    CFuint updatableID = 0;
    CFuint iProc = 0;
    bool found = false;

    for (CFuint ip = 0; ip < nbProcesses && !found; ++ip) {
      // skip the case ip == currRank, because you already know that
      // ghost points are not updatable in the current process
      if (ip != currRank) {
	const CFuint nbUpdatables = _updatablesPerRank[ip].size();
	for (CFuint j = 0; j < nbUpdatables; ++j) {
	  if (_updatablesPerRank[ip][j] == globGhostID) {
	    updatableID = j;
	    iProc = ip;
	    found = true;
	    break;
	  }
	}
      }
    }

    cf_assert(found);
    cf_assert(iProc != currRank);
    globalToPetscIDs.insert(globGhostID, startPetscID[iProc] + updatableID);
  }

  cf_assert(globalToPetscIDs.size() == nbPoints);

  // fundamental !!!!
  globalToPetscIDs.sortKeys();

  // finally the mapping from local IDs to Petsc global IDs is computed
  localToPetscIDs.resize(nbPoints);
  for (CFuint i = 0; i < nbPoints; ++i) {
    localToPetscIDs[i] = globalToPetscIDs.find(_globalIDs[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
