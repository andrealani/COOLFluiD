// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CFmeshFileWriter_hh
#define COOLFluiD_IO_CFmeshFileWriter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "CFmeshFileWriter/CFmeshFileWriterAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the CFmesh file OutputFormatter
  namespace CFmeshFileWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module CFmeshFileWriter
class CFmeshFileWriter_API CFmeshFileWriterModule : public Environment::ModuleRegister<CFmeshFileWriterModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {  return "CFmeshFileWriter";  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription() { return "This module implements the OutputFormatter that write CFmesh files.";   }

}; // end CFmeshFileWriterModule

//////////////////////////////////////////////////////////////////////////////

} // namespace CFmeshFileWriter
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_IO_CFmeshFileWriter_hh
