// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_MHDJacobGradientComputer_hh
#define COOLFluiD_FluxReconstructionMethod_MHDJacobGradientComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSJacobFluxReconstruction.hh"
#include "MHD/MHDProjectionDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Daughterclass of ConvRHSJacobFluxReconstruction, needed to calculate the 
 * gradients for implicit schemes for MHD
 * 
 * @author Ray Vandenhoeck
 */
class MHDJacobGradientComputer : public ConvRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit MHDJacobGradientComputer(const std::string& name);

  /// Destructor
  virtual ~MHDJacobGradientComputer() {}
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected: //functions
  
  /**
   * Compute the discontinuous contribution to the corrected gradients
   */
  virtual void computeGradients();
  
  /**
   * Compute the correction part of the corrected gradient
   */
  virtual void computeGradientFaceCorrections();
  
protected: //data
  
  /// diffusive variable set
  Common::SafePtr< Physics::MHD::MHDProjectionDiffVarSet > m_diffusiveVarSet;
  
  /// Vector transformer from update to solution variables
  Common::SafePtr<Framework::VarSetTransformer> m_updateToSolutionVecTrans;
  
  /// matrix to store the state terms needed for the gradients (p, u, v, T) inside element
  RealMatrix m_tempGradTerm;
  
  /// matrix to store the state terms needed for the gradients (p, u, v, T) for left neighbor
  RealMatrix m_tempGradTermL;
  
  /// matrix to store the state terms needed for the gradients (p, u, v, T) for right neighbor
  RealMatrix m_tempGradTermR;
  
  /// element states within an element in the correct format
  std::vector< RealVector* > m_tempStates;
  
  /// element states of the left neighbor in the correct format
  std::vector< RealVector* > m_tempStatesL;
  
  /// element states of the right neighbor in the correct format
  std::vector< RealVector* > m_tempStatesR;
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_MHDJacobGradientComputer_hh

