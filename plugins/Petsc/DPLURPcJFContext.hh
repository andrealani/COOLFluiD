// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_DPLURPcJFContext_hh
#define COOLFluiD_Numerics_Petsc_DPLURPcJFContext_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataStorage.hh"
#include "MathTools/RealMatrix.hh"
#include "Common/ConnectivityTable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Petsc {
    class JFContext;
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a tuple of data to be passed to the matrix free solver DP-LUR
 * user defined preconditioner of Petsc
 *
 * @author Jiri Simonek 
 *
 */
class DPLURPcJFContext {
public: // functions
  
  /// Constructor
  DPLURPcJFContext() : diagMatrices(CFNULL), upLocalIDsAll(CFNULL) {}
  
  /// handle of diagonal inverted matrices
  Common::SafePtr<Framework::DataSocketSink <CFreal> > diagMatrices;
  
  /// handle of local updatable IDs or -1 (ghost) for all local states
  Common::SafePtr<Framework::DataSocketSink <CFint> > upLocalIDsAll;
  
  /// pointer to JFContext - we will use bkpStates from this object during the DP-LUR preconditioning
  JFContext* pJFC;
  
	/// Omega - relaxation factor (for underrelaxation/overrelaxation)
	CFreal omega;
	
	/// Number of sweeps in the DP-LUR method
	CFuint nbSweeps;
	
  /// state neighbors connectivity
  Common::ConnectivityTable<CFuint> stateNeighbors;
    
}; // end of class DPLURPcJFContext

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_DPLURPcJFContext_hh
