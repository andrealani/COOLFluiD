// Copyright (C) 2019 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_NSBndGradientComputerCUDA_hh
#define COOLFluiD_FluxReconstructionMethod_NSBndGradientComputerCUDA_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstruction.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Daughterclass of ConvBndCorrectionsRHSFluxReconstruction, needed to 
 * calculate the bnd gradients for NS for CUDA
 * 
 * @author Ray Vandenhoeck
 */
class NSBndGradientComputerCUDA : public ConvBndCorrectionsRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit NSBndGradientComputerCUDA(const std::string& name);

  /// Destructor
  virtual ~NSBndGradientComputerCUDA() {}

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected: //functions
  
  /**
   * Compute the correction part of the corrected gradient for the bnd face
   */
  virtual void computeGradientBndFaceCorrections();
  
protected: //data

  /// socket for gradients cuda
  Framework::DataSocketSink< CFreal > socket_gradientsCUDA;
  
   /// socket for gradients cuda
  Framework::DataSocketSink< CFreal > socket_gradientsAVCUDA;
  
  /// diffusive variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_diffusiveVarSet;
  
  /// Vector transformer from update to solution variables
  Common::SafePtr<Framework::VarSetTransformer> m_updateToSolutionVecTrans;
  
  /// matrix to store the state terms needed for the gradients (p, u, v, T)
  RealMatrix m_tempGradTerm;
  
  /// matrix to store the ghost state terms needed for the gradients (p, u, v, T)
  RealMatrix m_tempGradTermGhost;
  
  /// element states within an element in the correct format
  std::vector< RealVector* > m_tempStates;
  
  /// element ghost states in the correct format
  std::vector< RealVector* > m_tempStatesGhost;
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_NSBndGradientComputerCUDA_hh

