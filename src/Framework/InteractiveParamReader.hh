// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_InteractiveParamReader_hh
#define COOLFluiD_Framework_InteractiveParamReader_hh

//////////////////////////////////////////////////////////////////////////////


#include <fstream>

#include "Common/StringOps.hh"
#include "Common/SafePtr.hh"
#include "Config/ConfigObject.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a reader for interactive parameters
/// @author Andrea Lani
/// @author Mehmet Sarp Yalim
/// @author Tiago Quintino
class Framework_API InteractiveParamReader : public Config::ConfigObject {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  ///  Constructor
  explicit InteractiveParamReader(const std::string& name);

  /// Default destructor
  virtual ~InteractiveParamReader();

  /// Configures the options for this object.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure ( Config::ConfigArgs& args );

  /// Reads the interactive configuration file
  void readFile();

private: // helper functions

  /// Read the file sequentially
  void readFileSequentially();

private:

  /// frequency of reading from the file containing
  /// interactive configuration parameters
  CFuint m_read_rate;

  /// name of the file from which reading all interactive
  /// parameters
  std::string m_filename;

}; // end of class InteractiveParamReader

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_InteractiveParamReader_hh
