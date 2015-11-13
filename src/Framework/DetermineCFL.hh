// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DetermineCFL_hh
#define COOLFluiD_Framework_DetermineCFL_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeCFL.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class allows the computation of the CFL using a vectorial function
/// @author Thomas Wuilbaut
class Framework_API DetermineCFL : public ComputeCFL {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  DetermineCFL(const std::string& name);

  /// Default destructor
  ~DetermineCFL();

  /// Check if the stop condition has been achieved
  void operator() (const ConvergenceStatus& m_cstatus);

  /// Configures this object with supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

private:

  /// a string holding the function definition
  std::string _function;
  
  /// ID of the subsystem to consider for the CFL calculation
  CFuint _subSystemID;
  
}; // end of class DetermineCFL
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DetermineCFL_hh
