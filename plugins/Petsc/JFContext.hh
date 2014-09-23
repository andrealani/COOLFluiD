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
  JFContext() : states(CFNULL), rhs(CFNULL), rhsVec(CFNULL) {}
  
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
  
  /// epsilon for the numerical derivative
  CFreal eps;
  
  /// Order of the Jacobian-free matrix vector product approximation (1 or 2)
  bool jfApprox2ndOrder;
  
  /// Enable/Disable usage of different preconditioner matrix
	bool differentPreconditionerMatrix;
  
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
