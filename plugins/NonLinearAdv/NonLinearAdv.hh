// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_NonLinearAdv_hh
#define COOLFluiD_Physics_NonLinearAdv_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Physics {

    /// The classes that implement a scalar rotation advection physical model.
  namespace NonLinearAdv {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module NonLinearAdv
 */
class NonLinearAdvModule : public Environment::ModuleRegister<NonLinearAdvModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "NonLinearAdv";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a scalar rotation advection physical model.";
  }

}; // end NonLinearAdvModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace NonLinearAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_NonLinearAdv_hh
