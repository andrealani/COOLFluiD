// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Paralution_ParalutionIdxMapping_hh
#define COOLFluiD_Numerics_Paralution_ParalutionIdxMapping_hh

//////////////////////////////////////////////////////////////////////////////

#include <valarray>

#include "Paralution/Paralution.hh"
#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

/**
 * Class responsible for the mapping of coolfluid state IDs
 * to state IDs that are used by the trilinos packages.
 * This class holds the information that allows each processor
 * to calculate the trilinos ID for each states that it knowns
 * (= locally owned + ghost states).
 *
 * @author Tim Boonen
 *
 */
class ParalutionIdxMapping  {
public:

  /**
   *  Registering all coolfluid state IDs that 'live' on a processor
   *  (= locally owned + ghost states)
   *  @param   myGlobalStateIDs  global coolfluid state IDs
   *  @param   isNonLocalState     flags telling if state is updatable
   *
   *  @post    after computeMapping(id):
   *             forall 0<=i<nbMyGlobalStateIDs: id[i] == myGlobalStateIDs[i]
   */
  ParalutionIdxMapping(const std::vector<int>& myGlobalStateIDs,
                     const std::vector<bool>& isNonLocalState) :
    _myGlobalStateIDs(myGlobalStateIDs),
    _isNonLocalState(isNonLocalState)
  {
    // check size
    cf_assert(_myGlobalStateIDs.size() == isNonLocalState.size());
  }

  /**
   * Destructor.
   */
  ~ParalutionIdxMapping()
  {
  }

  /**
   * Creates a mapping from local state IDs to the global state IDs used by Paralution
   * @post      localToParalutionIDs[localStateID] == globalParalutionStateID
   * @remarks   Also the local ghost state IDs are be mapped to their
   *            global trilinos state IDs
   */
  void computeMapping(std::valarray<CFuint>& localToParalutionIDs,
                      std::valarray<bool>& isNonLocalState);

private:

  /// The global state IDs seen by this processor (locally owned + ghost)
  std::vector<int> _myGlobalStateIDs;

  /// The global state IDs seen by this processor (locally owned + ghost)
  std::vector<bool> _isNonLocalState;

}; // end of class ParalutionIdxMapping

//////////////////////////////////////////////////////////////////////////////

    } // namespace Paralution

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Paralution_ParalutionIdxMapping_hh
