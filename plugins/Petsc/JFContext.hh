// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_JFContext_hh
#define COOLFluiD_Numerics_Petsc_JFContext_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/NumericalCommand.hh"
#include "Framework/DataStorage.hh"

#include "Petsc/PetscHeaders.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State;
    class SpaceMethod;
  }

  namespace Petsc {
    class PetscLSSData;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a tuple of data to be passed to the matrix free solver of
 * Petsc
 *
 * @author Andrea Lani
 * @author Jiri Simonek
 *
 */
class JFContext {
public: // functions

  /// Constructor
  JFContext() :
    states(CFNULL),
    rhs(CFNULL),
    rhsVec(CFNULL),
    stateVec(CFNULL),
    eps(1e-6),
    jfApprox2ndOrder(false),
    differentPreconditionerMatrix(false),
    useEisenstatWalker(true),
    prevNonlinResNorm(-1.0),
    prevEta(0.5),
    ewGamma(0.9),
    ewAlpha(2.0),
    ewMaxEta(0.9),
    ewMinEta(1e-4),
    ewLastTimeStep(0)
  {}

  /// pointer to the Petsc method data
  Common::SafePtr<PetscLSSData> petscData;

  /// handle of states
  Common::SafePtr<Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > > states;

  /// handle of rhs
  Common::SafePtr<Framework::DataSocketSink < CFreal > > rhs;

  /// handle of updateCoeff
  Common::SafePtr<Framework::DataSocketSink < CFreal > > updateCoeff;

  /// Petsc RHS vector
  Common::SafePtr<PetscVector> rhsVec;

  /// pointer to the SpaceMethod
  Common::SafePtr<Framework::SpaceMethod> spaceMethod;

  /// Epsilon for numerical derivative.
  /// DEPRECATED for ParJFSolveSys (MatMFFD uses adaptive epsilon).
  /// Still used by GMRESR variant and LUSGS/DPLUR preconditioners.
  CFreal eps;

  /// DEPRECATED: 2nd-order FD not supported with MatMFFD.
  /// Kept for backward compatibility — logged as deprecation warning if true.
  bool jfApprox2ndOrder;

  /// Enable/Disable usage of different preconditioner matrix
  bool differentPreconditionerMatrix;

  /// PETSc Vec for packing current states U (used by MatMFFDSetBase)
  Vec stateVec;

  // --- Eisenstat-Walker adaptive KSP tolerance ---

  /// Enable/Disable Eisenstat-Walker forcing term
  bool useEisenstatWalker;

  /// ||F(x_{k-1})|| from previous Newton step (-1 signals first iteration)
  CFreal prevNonlinResNorm;

  /// Previous forcing term (for safeguard)
  CFreal prevEta;

  /// EW parameter gamma (default 0.9)
  CFreal ewGamma;

  /// EW exponent alpha (default 2.0, Choice 2)
  CFreal ewAlpha;

  /// Maximum forcing term (default 0.9)
  CFreal ewMaxEta;

  /// Minimum forcing term (default 1e-4)
  CFreal ewMinEta;

  /// Last outer iteration number (for detecting new time step in EW reset)
  CFuint ewLastTimeStep;

  // --- End Eisenstat-Walker ---

  /// backup of the states array
  RealVector bkpStates;

  /// backup of the update coefficients array
  RealVector bkpUpdateCoeff;

  /// vector of the local ids of only the updatable states
  std::vector<CFint> upLocalIDs;

  /// indexes for the insertion of elements in a PetscVector
  std::vector<CFint> upStatesGlobalIDs;

}; // end of class JFContext

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_JFContext_hh
