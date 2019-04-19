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

/**
 * This is a standard command to deallocate data specific to a FluxReconstruction method
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 */
class StdUnSetup : public FluxReconstructionSolverCom {

public: // functions

  /// Constructor
  explicit StdUnSetup(const std::string& name);

  /// Destructor
  ~StdUnSetup() {}

  /// Execute processing actions
  void execute();
  
  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();
protected:
    
    /// unsetup nodes in states
    void unsetStateNodes();
    
protected: // protected data
  
  /// socket for unit normals in face flux points
  /// socket for size of projection vector in face flux points
  Framework::DataSocketSink<  std::vector< CFreal > > socket_faceJacobVecSizeFaceFlxPnts;
    
  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;
  
  /// socket for gradientsAV
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradientsAV;
  
  /// socket for positivity preservation values
  Framework::DataSocketSink< CFreal > socket_posPrev;

}; // class StdUnSetup

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_StdUnSetup_hh

