// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MarcoTest_MarcoTestModule_hh
#define COOLFluiD_MarcoTest_MarcoTestModule_hh

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace MarcoTest {

//////////////////////////////////////////////////////////////////////////////

/// This class defines an empty module
class MarcoTestModule : public Environment::ModuleRegister< MarcoTestModule > {
public:

  /// Static function that returns the module name. Must be implemented for the
  /// ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {
    return "MarcoTest";
  }

  /// Static function that returns the description of the module. Must be
  /// implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription() {
    return "This module interfaces and empty solver.";
  }

}; // end MarcoTestModule

//////////////////////////////////////////////////////////////////////////////

  }  // namespace MarcoTest
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MarcoTestSpcaeMethod_MarcoTestModule_hh

