// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullLimiter_hh
#define COOLFluiD_Framework_NullLimiter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Limiter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class implments a NULL limiter
/// @author Andrea Lani

template < typename METHODDATA >
class NullLimiter : public Limiter<METHODDATA> {
public:

  /// Constructor
  NullLimiter(const std::string& name) : Limiter<METHODDATA>(name)
  {
  }

  /// Default destructor
  ~NullLimiter()
  {
  }
  
  /// Checks if this object is a Null object.
  /// By default is always false, except for
  /// the concrete Null types that should implement
  /// it as returning true.
  /// @return true if Null and false otherwise
  bool isNull() const {return true;}
  
  /// Set private data that will be used during the computation
  void setup();

  /// Unsetup private data that will be used during the computation
  void unsetup();

  /// Compute the flux in the current face
  /// @post in limiterValue you put the value of the limiter for
  ///       each variable
  void limit(const std::vector<std::vector<Node*> >& coord,
       Framework::GeometricEntity* const cell,
       CFreal* limiterValue);

  /// Apply the face limiter
  void limitOnFace(const RealVector& rLeft,
    const RealVector& rRight,
    CFreal* limiterValue);

  /// Apply the limiter to a scalar quantity
  void limitScalar(CFreal r, CFreal& limiterValue);

}; // end of class NullLimiter

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NullLimiter.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullLimiter_hh
