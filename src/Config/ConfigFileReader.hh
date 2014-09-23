// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ConfigFileReader_hh
#define COOLFluiD_ConfigFileReader_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/path.hpp>

#include "Common/COOLFluiD.hh"
#include "Config/ConfigArgs.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

/// Parses a Config file
/// @author Tiago Quintino
class Config_API ConfigFileReader {
public: // functions

  /// Empty constructor
  ConfigFileReader();

  /// Destructor
  ~ConfigFileReader();

  /// Parses the file but first discovers which is the concrete type of the file
  /// and dispatches the parsing to the correct parser
  /// @param name path to the file
  /// @param args the argument map to fill in with the contents
  void parse (const boost::filesystem::path& name, ConfigArgs& args);

}; // Class ConfigFileReader

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ConfigFileReader_hh

