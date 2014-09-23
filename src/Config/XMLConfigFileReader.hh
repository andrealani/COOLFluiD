// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_XMLConfigFileReader_hh
#define COOLFluiD_XMLConfigFileReader_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/xmlParser.h"
#include "Config/ConfigArgs.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

/// Parses an XML config file
/// @author Tiago Quintino
class Config_API XMLConfigFileReader {
public: // functions

  /// Empty constructor
  XMLConfigFileReader();

  /// Destructor
  ~XMLConfigFileReader();

  /// Parses the xml file
  /// @param xcfcase is the file name
  /// @param args the argument map to fill in with the contents
  void parse(const std::string& xcfcase, ConfigArgs& args);

  /// Parses the xml string
  /// @param xcfcase is the file name
  /// @param args the argument map to fill in with the contents
  void parseString(const std::string& xml_string, ConfigArgs& args);

private: // functions

  /// Processes the informatio of an XML node
  /// This function is recursive.
  /// @param xNode is the XML node to process
  /// @param args the argument map to fill in with the contents
  /// @param nest the nested name derived from the parent XML nodes
  void process_xml_node (XMLNode xNode, ConfigArgs& args, std::string nest = std::string());

}; // Class XMLConfigFileReader

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_XMLConfigFileReader_hh

