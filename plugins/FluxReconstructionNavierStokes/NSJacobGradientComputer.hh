// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_NSJacobGradientComputer_hh
#define COOLFluiD_FluxReconstructionMethod_NSJacobGradientComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvRHSJacobFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Daughterclass of ConvRHSJacobFluxReconstruction, needed to calculate the 
 * gradients for implicit schemes for NS
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 */
class NSJacobGradientComputer : public ConvRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit NSJacobGradientComputer(const std::string& name);

  /// Destructor
  virtual ~NSJacobGradientComputer() {}

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
  Common::SafePtr< Framework::DiffusiveVarSet > m_diffusiveVarSet;
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_NSJacobGradientComputer_hh

