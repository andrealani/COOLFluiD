// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdvSys_hh
#define COOLFluiD_Physics_LinearAdvSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    
    namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

/**
 * @author Tiago Quintino
 * @author Tomas Kopacek
*/

/**
 * This class defines the Module System of Linear Advction Equations
 */
class LinearAdvSysModule : public Environment::ModuleRegister<LinearAdvSysModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "LinearAdvSys";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a model of a set of linear advection equations.";
  }

}; // end LinearAdvSysModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearAdvSys

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_LinEuler_hh
