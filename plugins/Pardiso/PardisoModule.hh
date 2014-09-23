// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Pardiso_PardisoModule_hh
#define COOLFluiD_Pardiso_PardisoModule_hh

#include "Environment/ModuleRegister.hh"

namespace COOLFluiD {
  namespace Pardiso {

/// This class defines the Module Pardiso
class PardisoModule : public Environment::ModuleRegister< PardisoModule > {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName() {
    return "Pardiso";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription() {
    return "This module interfaces the PARDISO linear system solver.";
  }

}; // end PardisoModule

  } // namespace Pardiso
} // namespace COOLFluiD

#endif // COOLFluiD_Pardiso_PardisoModule_hh

