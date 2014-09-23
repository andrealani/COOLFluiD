// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS2CFmesh_hh
#define COOLFluiD_CGNS2CFmesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement a converter from CGNS format to CFmesh.
  namespace CGNS2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module CGNS2CFmesh
 */
class CGNS2CFmeshModule : public Environment::ModuleRegister<CGNS2CFmeshModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "CGNS2CFmesh";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a converter from CGNS format to CFmesh.";
  }

}; // end CGNS2CFmeshModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNS2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_CGNS2CFmesh_hh
