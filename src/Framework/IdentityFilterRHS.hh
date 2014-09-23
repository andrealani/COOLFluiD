// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_IdentityFilterRHS_hh
#define COOLFluiD_Framework_IdentityFilterRHS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FilterRHS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a configurable filter for residual values
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API IdentityFilterRHS : public FilterRHS {
public:

  /// Constructor
  IdentityFilterRHS(const std::string& name);

  /// Constructor
  virtual ~IdentityFilterRHS();

  /// Filter the given value of the residual
  /// @param state state to be filtered
  virtual void filter (CFuint iVar, CFreal& dU);

}; // end of class IdentityFilterRHS

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_IdentityFilterRHS_hh
