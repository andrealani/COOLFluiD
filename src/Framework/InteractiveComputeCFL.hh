// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_InteractiveComputeCFL_hh
#define COOLFluiD_Framework_InteractiveComputeCFL_hh

//////////////////////////////////////////////////////////////////////////////

#include "Config/ConfigObject.hh"
#include "Framework/ComputeCFL.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class allows to change the CFL interactively
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API InteractiveComputeCFL : public ComputeCFL {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  InteractiveComputeCFL(const std::string& name);

  /// Default destructor
  virtual ~InteractiveComputeCFL();

  /// Computes the CFL value
  void operator() (const ConvergenceStatus& m_cstatus);

  /// Configures the options for this object.
  /// @param args is the ConfigArgs with the arguments to be parsed.
  virtual void configure ( Config::ConfigArgs& args );

private:

  /// current interactive value of the CFL
  CFreal m_curr_cfl;

  /// previous value that has been applied
  CFreal m_prev_cfl;

}; // end of class InteractiveComputeCFL

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_InteractiveComputeCFL_hh
