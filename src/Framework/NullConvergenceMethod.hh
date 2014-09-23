// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullConvergenceMethod_hh
#define COOLFluiD_Framework_NullConvergenceMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "ConvergenceMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NullConvergenceMethod.
/// @author Tiago Quintino
class Framework_API NullConvergenceMethod : public ConvergenceMethod {
public:

  /// Constructor.
  explicit NullConvergenceMethod(const std::string& name);

  /// Destructor
  virtual ~NullConvergenceMethod();

  /// Checks if this object is a Null object.
  /// Since this is NullSpaceMethod
  /// @return true
  virtual bool isNull() const {  return true;  }

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

protected: // abstract interface implementations

  /// Gets the Data aggregator of this Convergence method
  /// @return SafePtr to the ConvergenceMethodData
  virtual Common::SafePtr<Framework::ConvergenceMethodData> getConvergenceMethodData();

  /// Take one timestep
  /// @see ConvergenceMethod::takeStep()
  virtual void takeStepImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

}; // end NullConvergenceMethod

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ConvergenceMethod_hh
