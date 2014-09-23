// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_THOR2CFmesh_hh
#define COOLFluiD_IO_THOR2CFmesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement a converter from THOR format to CFmesh.
  namespace THOR2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module THOR2CFmesh
 */
class THOR2CFmeshModule : public Environment::ModuleRegister<THOR2CFmeshModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "THOR2CFmesh";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a converter from THOR format to CFmesh.";
  }

}; // end THOR2CFmeshModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace THOR2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_IO_THOR2CFmesh_hh
