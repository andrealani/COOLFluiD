// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFLUID_CFmeshExtruder_hh
#define COOLFLUID_CFmeshExtruder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace IO {

    /// The classes that implement a MeshFormatConverter
    /// that reads 2D CFmesh files, and writes back another 3D CFmesh file,
    /// by extruding the mesh in the z coordinate.
  namespace CFmeshExtruder {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module CFmeshExtruder
 */
class CFmeshExtruderModule : public Environment::ModuleRegister<CFmeshExtruderModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "CFmeshExtruder";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a MeshFormatConverter that reads 2D CFmesh files, and writes back another 3D CFmesh file.";
  }

}; // end CFmeshExtruderModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshExtruder

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_CFmeshExtruder_hh
