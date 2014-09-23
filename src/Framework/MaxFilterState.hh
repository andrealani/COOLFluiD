// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MaxFilterState_hh
#define COOLFluiD_Framework_MaxFilterState_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FilterState.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a configurable filter for residual values
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API MaxFilterState : public FilterState {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  MaxFilterState(const std::string& name);

  /// Destructor
  virtual ~MaxFilterState();

  /// Filter the given value of the residual
  /// @param state state to be filtered
  virtual void filter (RealVector& state) const;

private: // data

  /// array of minimum values allowed for the variable to filter
  std::vector<CFreal> m_minValues;

}; // end of class MaxFilterState

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MaxFilterState_hh
