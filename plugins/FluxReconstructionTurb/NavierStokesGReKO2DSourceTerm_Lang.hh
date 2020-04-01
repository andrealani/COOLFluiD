#ifndef COOLFluiD_FluxReconstructionMethod_NavierStokesGReKO2DSourceTerm_Lang_hh
#define COOLFluiD_FluxReconstructionMethod_NavierStokesGReKO2DSourceTerm_Lang_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FluxReconstructionTurb/KOmega2DSourceTerm.hh"
#include "Common/SafePtr.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "KOmega/NavierStokesKOmegaVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }
  
    namespace FluxReconstructionMethod {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the 2D Gamma-Re_theta-K-Omega Source Term for FR
 *
 * @author Ray Vandenhoeck
 *
 */
class NavierStokesGReKO2DSourceTerm_Lang : 
	public KOmega2DSourceTerm {
  
public:
  
  /**
   * Constructor
   */
  NavierStokesGReKO2DSourceTerm_Lang(const std::string& name);

  /**
   * Default destructor
   */
  ~NavierStokesGReKO2DSourceTerm_Lang();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();
  
  /**
   * add the source term
   */
  void addSourceTerm(RealVector& resUpdates);
  
private: // helper functions
  
  /// corresponding diffusive variable set
  void getVorticity(const CFuint iState);

  void getStrain(const CFreal VoverRadius, const CFuint iState);

  void getRethetac(const CFreal Retheta);

  void getFlength(const CFreal Retheta);

  void getRethetat(const CFreal Tu);  
  
  void getLambda(CFreal& lambda, const CFreal theta, const CFreal viscosity, const CFuint iState);

  void getFlambda(const CFreal lambda, const CFreal Tu, CFreal& Flambda, const CFreal theta, bool Prime );
  
  void getRethetatwithPressureGradient(const CFreal viscosity,const CFreal Tu, const CFuint iState);

private: // data
  
  ///K the Farfield
  CFreal m_kamb;

  ///Omega at the Farfield
  CFreal m_omegaamb;
  
  // Flag is true if we are using Vortivity source term
  // instead of the exact source term
  bool m_SST_V;  
  // Flag is true if we are adding a sustaining terma in order
  //  to control the non-physical decay of turbulence variables 
  //  in the freestream
  bool m_SST_sust;  
  
  CFreal   m_Rethetat; 
  CFreal   m_Rethetac; 
  CFreal   m_Flength; 
  CFreal   m_vorticity;
  CFreal  m_strain;
  bool     m_PGrad; 

}; // end of class NavierStokesGReKO2DSourceTerm_Lang

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_NavierStokesGReKO2DSourceTerm_hh
