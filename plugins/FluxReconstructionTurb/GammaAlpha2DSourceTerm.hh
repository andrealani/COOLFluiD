#ifndef COOLFluiD_FluxReconstructionMethod_GammaAlpha2DSourceTerm_hh
#define COOLFluiD_FluxReconstructionMethod_GammaAlpha2DSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FluxReconstructionTurb/KLogOmega2DSourceTerm.hh"
#include "Common/SafePtr.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "KOmega/NavierStokesKLogOmegaVarSetTypes.hh"

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
class GammaAlpha2DSourceTerm : 
	public KLogOmega2DSourceTerm {
  
public:
  
  /**
   * Constructor
   */
  GammaAlpha2DSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~GammaAlpha2DSourceTerm();
  
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
  
  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();
  
  virtual void getSToStateJacobian(const CFuint iState);
  
  virtual void getSToGradJacobian(const CFuint iState){};
  
  virtual bool isGradDependent(){return true;};
  
private: // helper functions
    
  /// storage for local freestream M
  Framework::DataSocketSource<CFreal> socket_MInfLocal;
  
  /// storage for local freestream u
  Framework::DataSocketSource<CFreal> socket_uInfLocal;
  
  /// storage for local freestream T
  Framework::DataSocketSource<CFreal> socket_TInfLocal;
  
  /// storage for local freestream Tu
  Framework::DataSocketSource<CFreal> socket_TuInfLocal;
  
  /// storage for alpha_c-alpha
  Framework::DataSocketSource<CFreal> socket_alphaDiff;
  
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal > socket_updateCoeff;
  
  /// corresponding diffusive variable set
  void getVorticity(const CFuint iState);

  void getStrain(const CFuint iState);

  CFreal getRethetat(const CFreal Tu, bool prime);  
  
  CFreal getRethetatwithPressureGradient(const CFreal avAlpha, const CFreal rhoInfLocal, const CFreal uInfLocal, const CFreal muInfLocal, const CFreal mu, const CFreal rho, const CFreal dpx, const CFreal dpy);

  CFreal getTuInfLocal(const CFreal Rethetat, const CFreal MInfLocal, const CFreal TuLocal);

private: // data
  
  CFreal   m_Rethetat; 
  CFreal   m_Rethetac; 
  CFreal   m_Flength; 
  CFreal   m_vorticity;
  CFreal  m_strain;
  bool     m_PGrad; 
  CFreal m_FLambda;
  bool m_decouple;
  bool m_limPRe;
  bool m_addUpdateCoeff;
  bool m_addDGDA;
  CFreal  m_lambdaLim;
  
  CFuint m_order;

}; // end of class GammaAlpha2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_NavierStokesGammaAlpha2DSourceTerm_hh
