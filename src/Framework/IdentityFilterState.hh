// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_IdentityFilterState_hh
#define COOLFluiD_Framework_IdentityFilterState_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FilterState.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a configurable filter for residual values
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API IdentityFilterState : public FilterState {
public:

  /// Constructor
  IdentityFilterState(const std::string& name);

  /// Constructor
  virtual ~IdentityFilterState();

  /// Filter the given value of the residual
  /// @param state state to be filtered
  virtual void filter (RealVector& state) const;

}; // end of class IdentityFilterState

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_IdentityFilterState_hh
