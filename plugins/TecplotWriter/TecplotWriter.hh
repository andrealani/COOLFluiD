// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_TecplotWriter_hh
#define COOLFluiD_IO_TecplotWriter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

#include "TecplotWriter/TecplotWriterAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the Tecplot writer
  namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module TecplotWriter
class TecplotWriter_API TecplotWriterModule : public Environment::ModuleRegister<TecplotWriterModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()
  {
    return "TecplotWriter";
  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implements the Tecplot writer";
  }

}; // end TecplotWriterModule

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_IO_TecplotWriter_hh
