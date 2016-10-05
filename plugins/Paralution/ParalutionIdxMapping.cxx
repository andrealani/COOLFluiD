// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Paralution/ParalutionIdxMapping.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
      
    namespace Paralution {
      
//////////////////////////////////////////////////////////////////////////////
      
void ParalutionIdxMapping::computeMapping
(std::valarray<CFuint>& localToParalutionStateIDs,
 std::valarray<bool>& isNonLocalState)
{
  // principle: since Epetra can deal with non-contiguous indices,
  //            the global coolfluid state IDs will do the job !!!
  const CFuint nbStates = _myGlobalStateIDs.size();
  localToParalutionStateIDs.resize(nbStates);
  isNonLocalState.resize(nbStates);
  
  for (CFuint i = 0; i< nbStates; i++) {
    localToParalutionStateIDs[i] = _myGlobalStateIDs[i];
    isNonLocalState[i] = _isNonLocalState[i];
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace Paralution

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
