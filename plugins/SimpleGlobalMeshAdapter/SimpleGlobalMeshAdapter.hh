// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace Numerics {

    /// The classes that implement SimpleGlobalMeshAdapter mesh adapter.
  namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module SimpleGlobalMeshAdapter
 */
class SimpleGlobalMeshAdapterModule : public Environment::ModuleRegister<SimpleGlobalMeshAdapterModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "SimpleGlobalMeshAdapter";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a simple global mesh adapter.";
  }

}; // end SimpleGlobalMeshAdapterModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_SimpleGlobalMeshAdapter_hh
