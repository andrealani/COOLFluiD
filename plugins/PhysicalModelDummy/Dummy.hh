// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_PhysicalModelDummy_hh
#define COOLFluiD_PhysicalModelDummy_hh

#include "Environment/ModuleRegister.hh"

namespace COOLFluiD {
  namespace PhysicalModelDummy {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module Dummy
class DummyModule : public Environment::ModuleRegister< DummyModule > {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "PhysicalModelDummy";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a dummy physical model.";
  }

}; // end DummyModule

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PhysicalModelDummy
}  // namespace COOLFluiD

#endif // COOLFLUID_PhysicalModelDummy_hh

