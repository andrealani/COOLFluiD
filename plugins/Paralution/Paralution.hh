// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Paralution_hh
#define COOLFluiD_Numerics_Paralution_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    /// The classes that implement the interface to Paralution Linear System Solver
  namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Paralution
 */
class ParalutionModule : public Environment::ModuleRegister<ParalutionModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Paralution";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements an interface to the Paralution linear system solver.";
  }

  /**
   * Initiates the Paralution environment
   */
  //virtual void initiate();

  /**
   * Terminates the Paralution environment
   */
  //virtual void terminate();

}; // end ParalutionModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace Paralution

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_Paralution_hh

