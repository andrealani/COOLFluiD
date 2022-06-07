// Copyright (C) 2016 KU Leuven, von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvDiffCLLAVJacobFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_ConvDiffCLLAVJacobFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/ConvDiffLLAVJacobFluxReconstruction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to assemble the convective, diffusive 
/// and artificial viscosity part of the system using a FluxReconstruction solver for an implicit scheme
/// @author Ray Vandenhoeck
class ConvDiffCLLAVJacobFluxReconstruction : public ConvDiffLLAVJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvDiffCLLAVJacobFluxReconstruction(const std::string& name);

  /// Destructor
  virtual ~ConvDiffCLLAVJacobFluxReconstruction() {}

  /// Execute processing actions
  virtual void execute();
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();
    
protected: //functions
  
  /**
   * Set the data for the current face necessary to calculate FI
   */
  virtual void setFaceData(CFuint faceID);

  /**
   * compute the contribution of the diffusive face term to both Jacobians
   */
  void computeBothJacobsDiffFaceTerm();

  /**
   * compute the contribution of the diffusive face term to one Jacobians
   */
  void computeOneJacobDiffFaceTerm(const CFuint side);
  
  /// compute the AV to state Jacobian
  virtual void computeEpsToStateJacobianAna() = 0;
  
  /// compute the cell CLLAV Jacobian (to state) analytically
  virtual void computeLLAVCellFluxJacobianAna(const CFreal resFactor);
  
  /// compute the Riemann CLLAV Jacobian (to state) analytically
  virtual void computeLLAVRiemannFluxJacobianAna(const CFreal resFactor);
  
  /**
   * Compute the artificial viscosity
   */
  virtual void computeEpsilon();
  
  /**
  * Store the computed artificial viscosity
  */
  virtual void storeEpsilon();
  
protected: //data
    
  /// correction function for LLAV for current cell
  std::vector< std::vector< CFreal > > m_corrFctLLAV;
  
  /// order of LLAV VCJH corr fct
  CFuint m_LLAVCorrFctOrder;
  
  /// VCJH factor of LLAV corr fct
  CFreal m_LLAVCorrFctFactor;
  
  /// alpha values for subcell AV redistribution for each side for each sol pnt 
  std::vector< std::vector< CFreal > > m_alphaValues;
  
  private:

  /// Physical data temporary vector
  RealVector m_pData;

}; // class ConvDiffCLLAVJacobFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvDiffCLLAVJacobFluxReconstruction_hh
