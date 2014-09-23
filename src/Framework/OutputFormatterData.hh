// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_OutputFormatterData_hh
#define COOLFluiD_Framework_OutputFormatterData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class DataHandleOutput;

//////////////////////////////////////////////////////////////////////////////

/// This class collects some data shared by all OutputFormatter's
/// @author Tiago Quintino
class Framework_API OutputFormatterData : public Framework::MethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  OutputFormatterData(Common::SafePtr<Method> owner);

  /// Default destructor
  virtual ~OutputFormatterData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// ???
  virtual Common::SafePtr<DataHandleOutput> getDataHOutput();

protected: // data

  DataHandleOutput * m_datah_out;

}; // end of class OutputFormatterData

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_OutputFormatterData_hh
