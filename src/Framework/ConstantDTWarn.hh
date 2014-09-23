// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ConstantDTWarn_hh
#define COOLFluiD_Framework_ConstantDTWarn_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ComputeDT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class allows the computation of the maximal allowable time step
/// @author Tiago Quintino

class Framework_API ConstantDTWarn : public ComputeDT {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  ConstantDTWarn(const std::string& name);

  /// Default destructor
  ~ConstantDTWarn();

  /// Check if the stop condition has been achieved
  void operator() ();

  /// Configures this object with supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

}; // end of class ConstantDTWarn

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ConstantDTWarn_hh
