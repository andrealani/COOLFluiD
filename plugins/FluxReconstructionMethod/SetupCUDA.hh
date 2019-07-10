// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_SetupCUDA_hh
#define COOLFluiD_FluxReconstructionMethod_SetupCUDA_hh

//////////////////////////////////////////////////////////////////////////////


#include "FluxReconstructionMethod/StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a command to setup the FluxReconstruction method if porting to GPU is needed
 * 
 * @author Ray Vandenhoeck
 */
class SetupCUDA : public StdSetup {

public: // functions

  /// Constructor
  explicit SetupCUDA(const std::string& name);

  /// Destructor
  ~SetupCUDA();

  /// Execute processing actions
  void execute();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  
private: // private functions    

  /**
   * Creates sockets that store the needed normals
   */
  void createNormalSockets();
  
protected: // data
  
  /// storage for the normals in the solution points
  Framework::DataSocketSource< CFreal > socket_solPntNormals;
  
};  // class SetupCUDA

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_SetupCUDA_hh


