// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullComputeCFL_hh
#define COOLFluiD_Framework_NullComputeCFL_hh

//////////////////////////////////////////////////////////////////////////////

#include "Config/ConfigObject.hh"
#include "ComputeCFL.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the condition to be satisfied to stop
/// the computation
/// @author Andrea Lani
class Framework_API NullComputeCFL : public ComputeCFL {
public:

  /// Constructor
  NullComputeCFL(const std::string& name);

  /// Default destructor
  ~NullComputeCFL();

  /// Check if the stop condition has been achieved
  void operator() (const ConvergenceStatus& m_cstatus);

}; // end of class NullComputeCFL

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullComputeCFL_hh
