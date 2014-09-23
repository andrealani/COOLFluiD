#ifndef COOLFluiD_Numerics_FiniteVolume_NavierStokesGReKO2DSourceTerm_Lang_hh
#define COOLFluiD_Numerics_FiniteVolume_NavierStokesGReKO2DSourceTerm_Lang_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteVolumeTurb/NavierStokesKOmega2DSourceTerm.hh"
#include "Common/SafePtr.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "KOmega/NavierStokesKOmegaVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the 2D Gamma-Re_theta-K-Omega Source Term
 *
 * @author Khalil Bensassi 
 *
 */
class NavierStokesGReKO2DSourceTerm_Lang : 
	public NavierStokesKOmega2DSourceTerm<Physics::KOmega::NavierStokes2DKOmega> {
  
public:
  
  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
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
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Compute the source term
   */
  void computeSource(Framework::GeometricEntity *const element,
  		     RealVector& source,
		     RealMatrix& jacobian);
  
private: // helper functions
  
  /// corresponding diffusive variable set
  CFreal GetVorticity();

  CFreal GetStrain(CFreal& VoverRadius);

  CFreal GetRethetac(CFreal& Retheta);

  CFreal GetFlength(CFreal& Retheta);

  CFreal GetRethetat(const CFreal& Tu);  
  
  CFreal GetLambda(CFreal& Lambda, CFreal& Theta, CFreal& Viscosity);

  CFreal GetFlambda(CFreal& Lambda,const CFreal& Tu, CFreal& Flambda, CFreal& Theta, bool Prime );
  CFreal GetRethetatwithPressureGradient(CFreal& Viscosity,const CFreal& Tu);

private: // data

  /// average State
  RealVector _avState;
  
  ///K the Farfield
  CFreal _kamb;

  ///Omega at the Farfield
  CFreal _omegaamb;
  
  // Flag is true if we are using Vortivity source term
  // instead of the exact source term
  bool _SST_V;  
  // Flag is true if we are adding a sustaining terma in order
  //  to control the non-physical decay of turbulence variables 
  //  in the freestream
  bool _SST_sust;  
  
  CFreal   _Rethetat; 
  CFreal   _Rethetac; 
  CFreal   _Flength; 
  CFreal   _Vorticity;
  CFreal  _Strain;
  bool     _PGrad; 

}; // end of class NavierStokesGReKO2DSourceTerm_Lang

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NavierStokesGReKO2DSourceTerm_hh
