// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_StdSolve_hh
#define COOLFluiD_FluxReconstructionMethod_StdSolve_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to assemble the system using FluxReconstruction solver
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class StdSolve : public FluxReconstructionSolverCom {

public: // functions

  /// Constructor
  explicit StdSolve(const std::string& name);

  /// Destructor
  virtual ~StdSolve() {}

  /// Execute processing actions
  void execute();
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected: //data
  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;
  
  /// socket for normals
  Framework::DataSocketSink< CFreal > socket_normals;
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_StdSolve_hh

