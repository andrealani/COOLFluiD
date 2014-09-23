// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_SAMGLSS_SAMGLSSModule_hh
#define COOLFluiD_SAMGLSS_SAMGLSSModule_hh

#include "Environment/ModuleRegister.hh"

namespace COOLFluiD {
  namespace SAMGLSS {

/// This class defines the Module SAMGLSS
class SAMGLSSModule : public Environment::ModuleRegister< SAMGLSSModule > {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName() {
    return "SAMGLSS";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription() {
    return "This module interfaces the SAMG linear system solver.";
  }

}; // end SAMGLSSModule

  }  // namespace SAMGLSS
}  // namespace COOLFluiD

#endif // COOLFluiD_SAMGLSS_SAMGLSSModule_hh

