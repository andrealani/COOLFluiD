// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_RungeKuttaLS_hh
#define COOLFluiD_Numerics_RungeKuttaLS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

  /// The classes that implement a low storage Runge-Kutta time stepper.
  namespace RungeKuttaLS {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module RungeKuttaLS
 */
class RungeKuttaLSModule : public Environment::ModuleRegister<RungeKuttaLSModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "RungeKuttaLS";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a low storage Runge-Kutta time stepper.";
  }

}; // end RungeKuttaLSModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKuttaLS

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_RungeKuttaLS_hh
