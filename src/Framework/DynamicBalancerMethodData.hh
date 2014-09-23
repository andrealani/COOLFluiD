// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DynamicBalancerMethodData_hh
#define COOLFluiD_Framework_DynamicBalancerMethodData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FilterState.hh"
#include "Framework/MethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Base class for the Data of the DynamicBalancerMethod's.
/// @author
class DynamicBalancerMethodData: public Framework::MethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  DynamicBalancerMethodData(Common::SafePtr<Method> owner);

  /// Destructor.
  ~DynamicBalancerMethodData();

  /// Configure the data from the supplied arguments.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

private: // data

}; // end of class DynamicBalancerMethodData

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DynamicBalancerMethodData_hh
