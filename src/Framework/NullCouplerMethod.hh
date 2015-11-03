// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullCouplerMethod_hh
#define COOLFluiD_Framework_NullCouplerMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "CouplerMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NullCouplerMethod.
/// @author Thomas Wuilbaut
class Framework_API NullCouplerMethod : public CouplerMethod {
public:

  /// Constructor.
  explicit NullCouplerMethod(const std::string& name);

  /// Destructor
  virtual ~NullCouplerMethod();

  /// Setting up communication between SubSystems
  /// by sending the info about the datahandles that will be
  /// transfered
  /// Since this is a Null Method it doesn't do anything but
  /// issue a warning.
  void setupCommunicationChannels();

  /// Checks if this object is a Null object.
  /// Since this is NullCouplerMethod
  /// @return true
  virtual bool isNull() const {  return true; }

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

protected: // abstract interface implementations

  /// Writes the data necessary for the preProcessing
  virtual void preProcessWriteImpl();

  /// Reads the data necessary for the preProcessing
  virtual void preProcessReadImpl();

  /// Writes the data necessary for the mesh matching
  virtual void meshMatchingWriteImpl();

  /// Reads the data necessary for the mesh matching
  virtual void meshMatchingReadImpl();

  /// Reads the data from the Coupled SubSystems
  virtual void dataTransferReadImpl();

  /// Writes the data for the Coupled SubSystems
  virtual void dataTransferWriteImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

  /// Finalize the coupling
  virtual void finalizeImpl();
  
}; // end NullCouplerMethod

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CouplerMethod_hh
