// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "TrilinosIdxMapping.hh"
#include "Common/CFMap.hh"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

#ifdef linefile
#undef linefile
#endif
#define linefile printf("Line %u in file %s.\n",__LINE__,__FILE__);fflush(stdout); //MPI_Barrier(MPI_COMM_WORLD);


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace Trilinos {
      
//////////////////////////////////////////////////////////////////////////////
      
void TrilinosIdxMapping::computeMapping(std::valarray<CFuint>& localToTrilinosStateIDs,
					std::valarray<bool>& isNonLocalState)
{
  // principle: since Epetra can deal with non-contiguous indices,
  //            the global coolfluid state IDs will do the job !!!
  const CFuint nbStates = _myGlobalStateIDs.size();
  localToTrilinosStateIDs.resize(nbStates);
  isNonLocalState.resize(nbStates);
  
  for (CFuint i = 0; i< nbStates; i++) {
    localToTrilinosStateIDs[i] = _myGlobalStateIDs[i];
    isNonLocalState[i] = _isNonLocalState[i];
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
