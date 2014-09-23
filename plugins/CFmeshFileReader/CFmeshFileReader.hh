// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshFileReader_hh
#define COOLFluiD_CFmeshFileReader_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "CFmeshFileReader/CFmeshFileReaderAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement MeshCreator
  /// that reads CFmesh files.
  namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module CFmeshFileReader
class CFmeshFileReader_API CFmeshFileReaderPlugin : public Environment::ModuleRegister<CFmeshFileReaderPlugin> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {  return "CFmeshFileReader";  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()  { return "This module implements the MeshCreator that read CFmesh files."; }

}; // end CFmeshFileReaderModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_CFmeshFileReader_hh
