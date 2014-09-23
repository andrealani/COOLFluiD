// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_NestedConfigFileReader_hh
#define COOLFluiD_NestedConfigFileReader_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Config {

    class ConfigArgs;

//////////////////////////////////////////////////////////////////////////////

/// Parses a config file in ASCII format
/// @author Tiago Quintino
class Config_API NestedConfigFileReader {
public:

  /// Reads the given file into a set of name/value pairs, both as std::strings.
  /// @param filename name of the file to parse
  explicit NestedConfigFileReader();

  /// Virtual destructor.
  virtual ~NestedConfigFileReader();

  /// Copys all of the name/value pairs into the supplied map of
  /// the labels and their values.
  /// @param m map where to place the labels and values of the options
  virtual void parse (const std::string& filename, ConfigArgs& args);

protected:

  /// The name of the file that was used to get these contents.
  /// @param fname missing documentation
  /// @return missing documentation
  std::string expandFileName(const std::string& fname) const;

private:

  /// Indicates to continue in nextline.
  static const CFchar CONTINUATOR;

  /// Separates tags and values.
  static const std::string SEPARATOR;

  /// Initiates comments. Rest of the line is ignored.
  static const std::string COMMENTOR;

  /// Initiates meta-comments
  static const std::string METACOMMENTOR;

  /// Initiates block comments. Starting from that point
  /// all the lines until BLOCKCOMMENTOREND is ignored.
  static const std::string BLOCKCOMMENTORBEGIN;
  static const std::string BLOCKCOMMENTOREND;

}; // Class NestedConfigFileReader

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_NestedConfigFileReader_hh

