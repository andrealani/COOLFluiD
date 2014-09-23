// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFLUID_CFmeshCellSplitter_hh
#define COOLFLUID_CFmeshCellSplitter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace IO {

    /// The classes that implement a MeshFormatConverter
    /// that reads a CFmesh files, and writes back another CFmesh file,
    /// with quads split into triangles or hexa/pyram split into tetrahedras.
  namespace CFmeshCellSplitter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines the Module CFmeshCellSplitter
 */
class CFmeshCellSplitterModule : public Environment::ModuleRegister<CFmeshCellSplitterModule> {
public:

  /**
   * Static function that returns the module name.
   * Must be implemented for the ModuleRegister template
   * @return name of the module
   */
  static std::string getModuleName()
  {
    return "CFmeshCellSplitter";
  }

  /**
   * Static function that returns the description of the module.
   * Must be implemented for the ModuleRegister template
   * @return descripton of the module
   */
  static std::string getModuleDescription()
  {
    return "This module implements a MeshFormatConverter that reads a CFmesh files, and writes back another CFmesh file, with quads split into triangles or hexa/pyram split into tetrahedras.";
  }

}; // end CFmeshCellSplitterModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshCellSplitter

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_CFmeshCellSplitter_hh
