// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_SetupCUDA_hh
#define COOLFluiD_FluxReconstructionMethod_SetupCUDA_hh

//////////////////////////////////////////////////////////////////////////////


#include "FluxReconstructionMethod/SetupExtra.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a command to setup the FluxReconstruction method if porting to GPU is needed
 * 
 * @author Ray Vandenhoeck
 */
class SetupCUDA : public SetupExtra {

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
  void createCUDASockets();
  
protected: // data

  /// storage for the gradients
  Framework::DataSocketSource< CFreal > socket_gradientsCUDA;
  
  /// storage for the gradients
  Framework::DataSocketSource< CFreal > socket_gradientsAVCUDA;

  /// storage for the face directions
  Framework::DataSocketSource< CFint > socket_faceDir;
  
};  // class SetupCUDA

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_SetupCUDA_hh


