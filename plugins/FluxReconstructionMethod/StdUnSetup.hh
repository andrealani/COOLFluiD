// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_StdUnSetup_hh
#define COOLFluiD_FluxReconstructionMethod_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to deallocate data specific to a FluxReconstruction method
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class StdUnSetup : public FluxReconstructionSolverCom {

public: // functions

  /// Constructor
  explicit StdUnSetup(const std::string& name) : FluxReconstructionSolverCom(name) {}

  /// Destructor
  ~StdUnSetup() {}

  /// Execute processing actions
  void execute();

}; // class StdUnSetup

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_StdUnSetup_hh

