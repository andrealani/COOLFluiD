// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MeshAdapterData_hh
#define COOLFluiD_Framework_MeshAdapterData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class collects some data shared by all MeshAdapter's
/// @author Tiago Quintino
class Framework_API MeshAdapterData : public Framework::MethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  MeshAdapterData(Common::SafePtr<Framework::Method> owner);

  /// Default destructor
  virtual ~MeshAdapterData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

protected: // data

}; // end of class MeshAdapterData

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MeshAdapterData_hh
