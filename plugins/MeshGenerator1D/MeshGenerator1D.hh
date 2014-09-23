// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MeshGenerator1D_MeshGenerator1D_hh
#define COOLFluiD_MeshGenerator1D_MeshGenerator1D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement a converter from Gambit neutral format to CFmesh.
  namespace MeshGenerator1D {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module MeshGenerator1D
 */
class MeshGenerator1DModule : public Environment::ModuleRegister<MeshGenerator1DModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "MeshGenerator1D";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a 1D mesh generator that writes a CFmesh file.";
  }

}; // end MeshGenerator1DModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshGenerator1D

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLU1D_MeshGenerator1D_MeshGenerator1D_hh
