// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Petsc_hh
#define COOLFluiD_Numerics_Petsc_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement the interface to Petsc Linear System Solver
  namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Petsc
 */
class PetscModule : public Environment::ModuleRegister<PetscModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Petsc";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an interface to the Petsc linear system solver.";
  }

  /**
   * Initiates the Petsc environment
   */
  virtual void initiate();

  /**
   * Terminates the Petsc environment
   */
  virtual void terminate();

}; // end PetscModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_Petsc_hh

