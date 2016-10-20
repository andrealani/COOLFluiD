// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_StdSetup_hh
#define COOLFluiD_FluxReconstructionMethod_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to setup the FluxReconstruction method
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class StdSetup : public FluxReconstructionSolverCom {

public: // functions

  /// Constructor
  explicit StdSetup(const std::string& name);

  /// Destructor
  ~StdSetup();

  /// Execute processing actions
  void execute();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();
  
protected: // data

  /// socket with proxy to be able to use the data handle of nodal states
  /// uniformly independently from the actual storage type being RealVector or
  /// State*
  Framework::DataSocketSource< Framework::ProxyDofIterator< RealVector >* >
    socket_nstatesProxy;

  /// socket for state's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// socket for gradients
  Framework::DataSocketSource< std::vector< RealVector > > socket_gradients;
  
  /// socket for normals
  Framework::DataSocketSource< CFreal > socket_normals;

};  // class StdSetup

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_StdSetup_hh

