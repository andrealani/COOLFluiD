// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_ShellPreconditioner_hh
#define COOLFluiD_Numerics_Petsc_ShellPreconditioner_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Common/SafePtr.hh"
#include "Common/CFMap.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/MethodStrategy.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

      class PetscLSSData;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a shell preconditioner object
 *
 * @author Andrea Lani
 *
 */

class ShellPreconditioner : public Framework::MethodStrategy<PetscLSSData> {
public:
  
  typedef Framework::BaseMethodStrategyProvider<PetscLSSData,ShellPreconditioner> PROVIDER;
  
  /**
   * Constructor
   */
  ShellPreconditioner(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ShellPreconditioner();

  /**
   * Set the preconditioner
   */
  virtual void setPreconditioner() = 0;
  
  /**
   * Compute before solving the system
   */
  virtual void computeBeforeSolving() {};
  
  /**
   * Compute after solving the system
   */
  virtual void computeAfterSolving() = 0;
  
  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ShellPreconditioner";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
}; // end of class ShellPreconditioner

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Petsc_ShellPreconditioner_hh
