// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullMeshAdapterMethod_hh
#define COOLFluiD_Framework_NullMeshAdapterMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "MeshAdapterMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NullMeshAdapterMethod.
/// @author Thomas Wuilbaut
class Framework_API NullMeshAdapterMethod : public MeshAdapterMethod {
public:

  /// Constructor.
  explicit NullMeshAdapterMethod(const std::string& name);

  /// Destructor
  virtual ~NullMeshAdapterMethod();

  /// Checks if this object is a Null object.
  /// Since this is NullMeshAdapterMethod
  /// @return true
  virtual bool isNull() const {  return true;  }

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

protected: // abstract interface implementations

  /// Adapts of the mesh
  /// @see MeshAdapterMethod::adaptMesh()
  virtual void adaptMeshImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

}; // end NullMeshAdapterMethod

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MeshAdapterMethod_hh
