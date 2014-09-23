// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullNullErrorEstimatorMethod_hh
#define COOLFluiD_Framework_NullNullErrorEstimatorMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "ErrorEstimatorMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class LinearSystemSolver;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NullNullErrorEstimatorMethod.
/// @author Jurek Majewski
class Framework_API NullErrorEstimatorMethod :
    public ErrorEstimatorMethod {

public:

  /// Constructor.
  explicit NullErrorEstimatorMethod(const std::string& name);

  /// Destructor
  virtual ~NullErrorEstimatorMethod();

  /// Configures this Method.
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Checks if this object is a Null object.
  /// @return true since this is NullErrorEstimatorMethod
  virtual bool isNull() const {  return true;  }

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

protected: // abstract interface implementations

  /// Estimate the Error
  /// @see ErrorEstimatorMethod::estimate()
  virtual void estimateImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

}; // class NullErrorEstimatorMethod

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullErrorEstimatorMethod_hh
