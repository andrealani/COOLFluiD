// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullComputeDT_hh
#define COOLFluiD_Framework_NullComputeDT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Config/ConfigObject.hh"
#include "ComputeDT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is the null class for the computation of the TimeStep
/// @author Andrea Lani

class Framework_API NullComputeDT : public ComputeDT {
public:

  /// Constructor
  NullComputeDT(const std::string& name);

  /// Default destructor
  ~NullComputeDT();

  /// Check if the stop condition has been achieved
  void operator() ();

}; // end of class NullComputeDT

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullComputeDT_hh
