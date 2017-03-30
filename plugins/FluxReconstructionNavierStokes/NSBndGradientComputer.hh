// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_NSBndGradientComputer_hh
#define COOLFluiD_FluxReconstructionMethod_NSBndGradientComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// Daughterclass of ConvBndCorrectionsRHSFluxReconstruction, needed to calculate the gradients
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class NSBndGradientComputer : public ConvBndCorrectionsRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit NSBndGradientComputer(const std::string& name);

  /// Destructor
  virtual ~NSBndGradientComputer() {}

protected: //functions
  
  /**
   * Compute the correction part of the corrected gradient for the bnd face
   */
  virtual void computeGradientBndFaceCorrections();
  
protected: //data
  
  /// diffusive variable set
  Common::SafePtr< Framework::DiffusiveVarSet > m_diffusiveVarSet;
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_NSBndGradientComputer_hh

