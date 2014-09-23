// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullDataProcessing_hh
#define COOLFluiD_Framework_NullDataProcessing_hh

//////////////////////////////////////////////////////////////////////////////

#include "DataProcessingMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NullDataProcessing.
/// @author Thomas Wuilbaut
class Framework_API NullDataProcessing : public DataProcessingMethod {
public:

  /// Constructor.
  explicit NullDataProcessing(const std::string& name);

  /// Destructor
  virtual ~NullDataProcessing();

  /// Sets the SpaceMethod
  void setCollaborator(MultiMethodHandle<SpaceMethod> spaceMtd);
  /// Sets the ConvergenceMethod
  void setCollaborator(MultiMethodHandle<ConvergenceMethod> convergenceMtd);

  /// Checks if this object is a Null object.
  /// @return true since this is NullDataProcessing
  bool isNull() const {  return true; }

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

protected:

  /// Process the data
  void processDataImpl();

  /// Sets up the data, commands and strategies of this Method
  /// @see Method::setMethod()
  virtual void setMethodImpl();

  /// Unsets the data, commands and strategies of this Method
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

}; // class NullDataProcessing

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullDataProcessing_hh
