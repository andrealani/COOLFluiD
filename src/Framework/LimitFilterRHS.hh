// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_LimitFilterRHS_hh
#define COOLFluiD_Framework_LimitFilterRHS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FilterRHS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a configurable filter for residual values
/// @author Andrea Lani
class Framework_API LimitFilterRHS : public FilterRHS {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  LimitFilterRHS(const std::string& name);

  /// Destructor
  virtual ~LimitFilterRHS();

  /// Filter the given value of the residual
  /// @param iEq    equation ID
  /// @param value  value to filter
  virtual void filter (CFuint iVar, CFreal& dU);
  
private: // data
  
  /// Maximum allowable value to impose to the RHS component
  std::vector<CFreal> m_maxDelta;
  
}; // end of class LimitFilterRHS

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_LimitFilterRHS_hh
