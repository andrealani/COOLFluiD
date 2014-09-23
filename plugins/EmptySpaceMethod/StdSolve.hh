// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_EmptySpaceMethod_StdSolve_hh
#define COOLFluiD_EmptySpaceMethod_StdSolve_hh

//////////////////////////////////////////////////////////////////////////////

#include "EmptySpaceMethod/EmptySolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace EmptySpaceMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to assemble the system using Empty solver
/// @author Tiago Quintino
/// @author Pedro Maciel
class StdSolve : public EmptySolverCom {

public: // functions

  /// Constructor
  explicit StdSolve(const std::string& name) : EmptySolverCom(name) {}

  /// Destructor
  virtual ~StdSolve() {}

  /// Execute processing actions
  void execute();

}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace EmptySpaceMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_EmptySpaceMethod_StdSolve_hh

