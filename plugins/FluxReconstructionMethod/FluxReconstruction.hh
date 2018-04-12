// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_FluxReconstructionModule_hh
#define COOLFluiD_FluxReconstructionMethod_FluxReconstructionModule_hh

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a flux reconstruction module
class FluxReconstructionModule : public Environment::ModuleRegister< FluxReconstructionModule > {
public:

  /// Static function that returns the module name. Must be implemented for the
  /// ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() 
  {
    return "FluxReconstruction";
  }

  /// Static function that returns the description of the module. Must be
  /// implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription() 
  {
    return "This module implements a flux reconstruction solver.";
  }

}; // end FluxReconstructionModule

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_FluxReconstructionModule_hh

