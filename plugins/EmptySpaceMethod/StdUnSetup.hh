// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_EmptySpaceMethod_StdUnSetup_hh
#define COOLFluiD_EmptySpaceMethod_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "EmptySpaceMethod/EmptySolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace EmptySpaceMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to deallocate data specific to an empty method
/// @author Tiago Quintino
/// @author Pedro Maciel
class StdUnSetup : public EmptySolverCom {

public: // functions

  /// Constructor
  explicit StdUnSetup(const std::string& name) : EmptySolverCom(name) {}

  /// Destructor
  ~StdUnSetup() {}

  /// Execute processing actions
  void execute();

}; // class StdUnSetup

//////////////////////////////////////////////////////////////////////////////

  }  // namespace EmptySpaceMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_EmptySpaceMethod_StdUnSetup_hh

