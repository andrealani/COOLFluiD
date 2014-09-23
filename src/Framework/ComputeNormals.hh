// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ComputeNormals_hh
#define COOLFluiD_Framework_ComputeNormals_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NotImplementedException.hh"

#include "Common/StringOps.hh"
#include "Common/SafePtr.hh"
#include "Common/OwnedObject.hh"

#include "Environment/ConcreteProvider.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Functor that computes the face (outward) normals
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API ComputeNormals : public Common::OwnedObject {
public:

  typedef Environment::ConcreteProvider<ComputeNormals> PROVIDER;

  /// Constructor
  ComputeNormals() : Common::OwnedObject() {}

  /// Virtual destructor
  virtual ~ComputeNormals() {}

  /// Overloading of the operator () to make this class act as a functor
  /// @param iFirstElem index of the first element of this kind
  /// @param iLastElem  index of the last element of this kind
  /// @param iType      type ID of the element
  virtual void operator() (const CFuint& iFirstElem, const CFuint& iLastElem, const CFuint& iType = 0) = 0;

  /// Update the Normals
  /// @param iFirstElem index of the first element of this kind
  /// @param iLastElem  index of the last element of this kind
  /// @param iType      type ID of the element
  virtual void update ( const CFuint& iFirstElem, const CFuint& iLastElem, const CFuint& iType = 0)
  {
    throw Common::NotImplementedException (FromHere(),getClassName() + "::update()");
  }

  /// Compute the averaged Normals between pastStates and Current States
  /// @todo make this more general !!!
  virtual void average (const CFuint& iFirstElem, const CFuint& iLastElem, const CFuint& iType = 0)
  {
    throw Common::NotImplementedException (FromHere(),getClassName() + "::average()");
  }

  /// Gets the Class name
  static std::string getClassName ()  {  return "ComputeNormals"; }

}; // end of class ComputeNormals

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(ComputeNormals) // define the factory

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeNormals_hh
