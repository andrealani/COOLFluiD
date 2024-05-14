// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CGNSWriter_hh
#define COOLFluiD_IO_CGNSWriter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

#include "CGNSWriter/CGNSWriterAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the CGNS writer
  namespace CGNSWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module CGNSWriter
class CGNSWriter_API CGNSWriterModule : public Environment::ModuleRegister<CGNSWriterModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "CGNSWriter";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implements the CGNS writer";
  }

}; // end CGNSWriterModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNSWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_IO_CGNSWriter_hh
