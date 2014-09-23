// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MaxComputeDT_hh
#define COOLFluiD_Framework_MaxComputeDT_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeDT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class allows the computation of the maximal allowable time step
/// @author Thomas Wuilbaut
/// @author Tiago Quintino
class Framework_API MaxComputeDT : public ComputeDT {
public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  MaxComputeDT(const std::string& name);

  /// Default destructor
  ~MaxComputeDT();

  /// Check if the stop condition has been achieved
  void operator() ();

  /// Configures this object with supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

private: // data

  /// Ratio between DT and the Maximum DT
  CFreal m_dt_ratio;

  /// reset the last DT to finish exactly
  bool m_exact_finish;

}; // end of class MaxComputeDT

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MaxComputeDT_hh
