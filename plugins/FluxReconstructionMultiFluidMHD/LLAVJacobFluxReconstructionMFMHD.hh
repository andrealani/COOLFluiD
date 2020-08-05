// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_LLAVJacobFluxReconstructionMFMHD_hh
#define COOLFluiD_FluxReconstructionMethod_LLAVJacobFluxReconstructionMFMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "FluxReconstructionNavierStokes/LLAVJacobFluxReconstructionNS.hh"

#include "Framework/MultiScalarTerm.hh"

#include "NavierStokes/MultiScalarVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
    
    namespace Maxwell {
      class Maxwell2DProjectionVarSet;
    }
  }
  
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Command to add Localized Laplacian Artificial Viscosity near discontinuities 
 * for implicit schemes for MFMHD
 * 
 * @author Ray Vandenhoeck
 */
class LLAVJacobFluxReconstructionMFMHD : public LLAVJacobFluxReconstructionNS {

public: // functions

  /// Constructor
  explicit LLAVJacobFluxReconstructionMFMHD(const std::string& name);

  /// Destructor
  virtual ~LLAVJacobFluxReconstructionMFMHD() {}
  
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
   * Compute the smoothness indicator
   */
  virtual void computeSmoothness();
  
  /**
   * Compute the smoothness indicator
   */
  virtual void computeSmoothness(const CFuint side);


protected: //data
  
  /// physical model var set
  Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet<Physics::Maxwell::Maxwell2DProjectionVarSet> > m_varSet;
  
  private:

  /// Physical data temporary vector
  RealVector m_pData;
  /// Physical data temporary vector
  RealVector m_pData2;
  
  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_LLAVJacobFluxReconstructionMFMHD_hh

