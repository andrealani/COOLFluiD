// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ProxyDofIterator_hh
#define COOLFluiD_Framework_ProxyDofIterator_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a proxy for a degree of freedom iterator.
/// It offers a uniform interface to use underlying different
/// storages.
/// @author Andrea Lani
template <class RETURNTYPE>
class ProxyDofIterator {
public:

  /// Default constructor
  ProxyDofIterator()
  {
  }

  /// Default destructor
  virtual ~ProxyDofIterator()
  {
  }

  /// Gets the dof corresponding to the nodes in the element
  /// Node* and State* can both be converted to RealVector*
  virtual RETURNTYPE* getState(const CFuint nodeID) = 0;
  
  /// Gets the dof corresponding to the nodes in the element
  /// Node* and State* can both be converted to RealVector*
  virtual RETURNTYPE* getNode(const CFuint nodeID) = 0;
  
  /// Gets the ID of the node in this dof
  virtual CFuint getNodeLocalID(const CFuint nodeID) = 0;

  /// Gets the ID of the node in this dof
  virtual CFuint getStateLocalID(const CFuint nodeID) = 0;

  /// Gets the dof corresponding to thenodes in the element
  virtual CFuint getSize() const = 0;

}; // end of class ProxyDofIterator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ProxyDofIterator_hh
