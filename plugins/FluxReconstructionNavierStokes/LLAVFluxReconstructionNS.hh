// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_LLAVFluxReconstructionNS_hh
#define COOLFluiD_FluxReconstructionMethod_LLAVFluxReconstructionNS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "FluxReconstructionMethod/LLAVFluxReconstruction.hh"

#include "Framework/MultiScalarTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class EulerVarSet;
      class EulerTerm;
    }
  }
  
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Command to add Localized Laplacian Artificial Viscosity near discontinuities 
 * for NS
 * 
 * @author Ray Vandenhoeck
 */
class LLAVFluxReconstructionNS : public LLAVFluxReconstruction {

public: // functions

  /// Constructor
  explicit LLAVFluxReconstructionNS(const std::string& name);

  /// Destructor
  virtual ~LLAVFluxReconstructionNS() {}
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
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
   * Set the data for the current cell necessary to calculate the residual update
   */
  virtual void setCellData();
  
  /**
   * Compute the Peclet number based on the user input
   */
  virtual CFreal computePeclet();
  
  /**
   * command to compute the coefficient to be multiplied with epsilon for the wavespeedupdates
   */
  virtual CFreal computeViscCoef(RealVector* state);

protected: //data
  
  /// corresponding diffusive variable set
  Common::SafePtr<Framework::DiffusiveVarSet> m_diffVarSet;
  
  /// backup of the gradients in the neighbouring cell
  std::vector< std::vector< std::vector< RealVector* > > > m_gradsBackUp;
  
  /// physical model
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> m_eulerVarSet;
  
  /// physical model MS
  Common::SafePtr< Framework::MultiScalarTerm< Physics::NavierStokes::EulerTerm > > m_msEulerTerm;
  
  /// number of species
  CFuint m_nbrSpecies;
  
  private:

  /// Physical data temporary vector
  RealVector m_pData;
  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_LLAVFluxReconstructionNS_hh

