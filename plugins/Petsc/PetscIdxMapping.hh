// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_PetscIdxMapping_hh
#define COOLFluiD_Numerics_Petsc_PetscIdxMapping_hh

//////////////////////////////////////////////////////////////////////////////

#include "Petsc/PetscHeaders.hh" // must come before any header

#include <valarray>

#include "Common/COOLFluiD.hh"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the index mapping from global to Petsc IDs in parallel
 *
 * @author Andrea Lani
 *
 */
class PetscIdxMapping  {
public:

  PetscIdxMapping(const std::vector< std::vector<CFuint> >& updatablesPerRank,
		  const std::vector<CFuint>& globalIDs,
		  const std::vector<bool>& isGhost) :
    _updatablesPerRank(updatablesPerRank),
    _globalIDs(globalIDs),
    _isGhost(isGhost)
  {
  }

  /**
   * Destructor.
   */
  ~PetscIdxMapping()
  {
  }

  /**
   * Create a mapping from local IDs to global Petsc IDS
   * @param nbUpdatablesPerRank array storing the number of
   *                            updatable points per process
   * @param globalIDs array storing the globall Ids of all the points
   *                  (ghost and non) in this process
   * @param isGhost array of flags indicating the ghost points
   *
   * @pre isGhost.size() == globalIDs.size()
   *
   * @post for the moment a mapping global IDs to global Petsc IDs
   *       is created but could be replaced with a mappping localID to
   *       global Petsc IDs (much more efficient in speed and memory)
   */
  void computeMapping(std::valarray<CFuint>& localToPetscIDs, std::valarray<bool>& isNonLocalRow);

private:

  /// updatable global state IDs per process
  const std::vector< std::vector<CFuint> >& _updatablesPerRank;

  /// all global IDs of the points in this processor
  const std::vector<CFuint>& _globalIDs;

  /// array of flag to telling if a point is ghost (non updatable)
  /// in this process
  const std::vector<bool>& _isGhost;

}; // end of class PetscIdxMapping

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_PetscIdxMapping_hh
