// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ErrorEstimatorData_hh
#define COOLFluiD_Framework_ErrorEstimatorData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class collects some data shared by all ErrorEstimator's
/// @author Tiago Quintino
class Framework_API ErrorEstimatorData : public Framework::MethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  ErrorEstimatorData(Common::SafePtr<Method> owner);

  /// Default destructor
  virtual ~ErrorEstimatorData();

  /// Configure this object and the nested objects
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

protected: // data

}; // end of class ErrorEstimatorData

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ErrorEstimatorData_hh
