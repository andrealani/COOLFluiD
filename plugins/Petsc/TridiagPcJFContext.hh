// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_TridiagPcJFContext_hh
#define COOLFluiD_Numerics_Petsc_TridiagPcJFContext_hh

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
 * This class represents a tuple of data to be passed to the matrix free solver
 * user defined tridiagonal preconditioner of Petsc
 *
 * @author Jiri Simonek
 *
 */
class TridiagPcJFContext {
public: // functions

  /// Constructor
  TridiagPcJFContext() : diagMatrices(CFNULL),
                         underDiagMatrices(CFNULL),
                         aboveDiagMatrices(CFNULL),
                         upLocalIDsAll(CFNULL),
                         dplrIsFirstInLine(CFNULL),
                         dplrToLocalIDs(CFNULL),
                         localToDplrIDs(CFNULL),
                         dplrCellInLine(CFNULL) {}

  /// pointer to the Petsc method data
  Common::SafePtr<PetscLSSData> petscData; /// debug

  /// handle of states
  Common::SafePtr<Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > > states; /// debug

  /// handle of diagonal jacobians
  Common::SafePtr<Framework::DataSocketSink<CFreal> > diagMatrices;

  /// handle of under-diagonal jacobians
  Common::SafePtr<Framework::DataSocketSink<CFreal> > underDiagMatrices;

  /// handle of above-diagonal jacobians
  Common::SafePtr<Framework::DataSocketSink<CFreal> > aboveDiagMatrices;

  /// handle of local updatable IDs or -1 (ghost) for all local states
  Common::SafePtr<Framework::DataSocketSink <CFint> > upLocalIDsAll;

  /// handle of bool array, which indicates if the cell is first in the line (DPLR)
  Common::SafePtr<Framework::DataSocketSink <bool> > dplrIsFirstInLine;

  /// handle of DPLR IDs -> local IDs conversion array
  Common::SafePtr<Framework::DataSocketSink <CFint> > dplrToLocalIDs;

  /// handle of local IDs -> DPLR IDs conversion array
  Common::SafePtr<Framework::DataSocketSink <CFint> > localToDplrIDs;

  /// handle of array which indicates to which line the cell belongs (DPLR) - may be it is useless
  Common::SafePtr<Framework::DataSocketSink <CFint> > dplrCellInLine;

  // pointer to JFContext - we will use bkpStates from this object during the LU-SGS preconditioning
  // JFContext* pJFC;

}; // end of class TridiagPcJFContext

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_TridiagPcJFContext_hh
