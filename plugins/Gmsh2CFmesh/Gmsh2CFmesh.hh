// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_Gmsh2CFmesh_hh
#define COOLFluiD_IO_Gmsh2CFmesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace IO {

  /// The classes that implement a converter from Gmsh format to CFmesh.
  namespace Gmsh2CFmesh {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module Gmsh2CFmesh
 */
class Gmsh2CFmeshModule : public Environment::ModuleRegister<Gmsh2CFmeshModule> {
public:
  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "Gmsh2CFmesh";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a converter from Gmsh format to CFmesh.";
  }

}; // end Gmsh2CFmeshModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace Gmsh2CFmesh

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_IO_Gmsh2CFmesh_hh

#ifndef COOLFluiD_IO_Gmsh2CFmesh_hh
#define COOLFluiD_IO_Gmsh2CFmesh_hh

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    /// The classes that implement the Gmsh2CFmesh file converter.
    namespace Gmsh2CFmesh {}

  }  // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_Gmsh2CFmesh_hh
