// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_FAST2CFmesh_hh
#define COOLFluiD_IO_FAST2CFmesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement a converter from FAST format to CFmesh.
  namespace FAST2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module FAST2CFmesh
 */
class FAST2CFmeshModule : public Environment::ModuleRegister<FAST2CFmeshModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "FAST2CFmesh";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a converter from FAST format to CFmesh.";
  }

}; // end FAST2CFmeshModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace FAST2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_IO_FAST2CFmesh_hh
