#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VarSetTransformer.hh"
#include "FiniteVolume/FVMCC_FluxSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Lax-Friedrichs flux
 *
 * @author Andrea Lani
 * @author Khalil Bensassi
 *
 */
template <class UPDATEVAR>
class AUSMFlux : public FVMCC_FluxSplitter {
public:
  
  /**
   * Constructor
   */
  AUSMFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMFlux();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configure the data from the supplied arguments.
   * @param args configuration arguments
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    FVMCC_FluxSplitter::configure(args);
  }
  
  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
  /**
   * Compute the left flux jacobian
   */
  // virtual void computeLeftJacobian();
  
  /**
   * Compute the right flux jacobian
   */
  // virtual void computeRightJacobian();
  
protected:
    
  /**
   * Compute the flux
   */
  virtual void computeMassFluxImpl(const CFuint eulerID,
				   const CFuint nbEulerEqs,
				   const std::vector<CFuint>& eulerVarIDs,
				   RealVector& result);
  
  /**
   * Compute the interface mass flux
   */
  virtual void computeMassFlux() = 0;
  
  /**
   * Compute the interface pressure flux
   */
  virtual void computePressureFlux() = 0;
  
   /**
   * Compute correction term for Incompressible flows
   */
  virtual void computeIncompCorrectionTerm() {m_addIncompCorrection = false;}
  
  /**
   * Compute the interface sound speed
   */
  void computeInterfaceSoundSpeed();
  
  /**
   * Compute the update coefficient in the standard way
   */
  void computeUpdateCoeff();

  /**
   * Compute the update coefficient in Liou's way from lambda+ and lambda-
   */
  void computeLiouUpdateCoeff();
  
  /**
   * Compute the interface sound speed in the 1st way
   */
  void computeSoundSpeed1();
  
  /**
   * Compute the interface sound speed in the 2nd way
   */
  void computeSoundSpeed2();
  
  /**
   * Compute the interface sound speed in the 3rd way
   */
  void computeSoundSpeed3();
  
  /**
   * Compute the interface sound speed in the 4th way
   */
  void computeSoundSpeed4();
  
  /**
   * Compute the interface sound speed in the 5th way (this should be chosen for TCNEQ flows)
   */
  void computeSoundSpeed5();
  
  /**
   * Applies the function Mach1+(M)
   */
  CFreal mach1Plus(const CFreal mach) const 
  {
    return 0.5*(mach + std::abs(mach));
  }
  
  /**
   * Applies the function Mach1-(M)
   */
  CFreal mach1Min(const CFreal mach) const 
  {
    return 0.5*(mach - std::abs(mach));
  }
  
  /**
   * Applies the function Mach2+(M)
   */
  CFreal mach2Plus(const CFreal mach) const 
  {
    return 0.25*pow(mach + 1.0, 2.0);
  }
  
  /**
   * Applies the function Mach2-(M)
   */
  CFreal mach2Min(const CFreal mach) const 
  {
    return -0.25*pow(mach - 1.0, 2.0);
  }
  
  /// Compute abs of the linearized flux jacobian
  void computeLinearizedAbsJacob();
  
  /// Tell if the incompressible correction needs to be applied
  bool addIncompCorrection() const {return m_addIncompCorrection;}
  
  /// COmpute the flux in a decoupled case 
  virtual void computeDecoupled(RealVector& result);
  
protected:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> m_updateVarSet;
  
  /// flag telling if the incompressible correction must be added
  bool m_addIncompCorrection;  
  
  /// left data  
  RealVector* m_lData;
  
  /// right data
  RealVector* m_rData;
    
  /// normal velocity left
  CFreal m_unL;
  
  /// normal velocity right
  CFreal m_unR;
  
  /// interface sound speed
  CFreal m_a12;
  
  /// mach number left
  CFreal m_mL;
  
  /// mach number right
  CFreal m_mR;
    
  /// interface mass flux
  CFreal m_mflux12;
  
  /// interface pressure flux
  CFreal m_p12;
  
  /// interface pressure flux
  CFreal m_mincomp;
  
  /// temporary unit normal
  RealVector m_tempUnitNormal;
  
  /// array of physical data 
  RealVector m_pdata;
  
  /// matrix of right eigenvectors
  RealMatrix   _rightEv;

  /// matrix of left eigenvectors
  RealMatrix   _leftEv;

  /// vector of eigenvalues
  RealVector   _eValues;
  
  /// vector of eigenvalues
  RealVector   _absEvalues;

  /// abs of the jacobian matrix
  RealMatrix   _absJacob;
  
  /// right jacobian matrix
  RealMatrix   _jRight;
  
  /// left jacobian matrix
  RealMatrix   _jLeft;
  
  /// jacobian matrix
  RealMatrix   _jacob;
  
  /// jacobian matrix
  RealMatrix   _jacobLeftTransf;
  
  /// jacobian matrix
  RealMatrix   _jacobRightTransf;
  
  /// dummy jacobian matrix
  RealMatrix   _jacobDummy;
  
  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;
  
  /// user defined coefficient for the calculation method of a12
  CFuint m_choiceA12;

  /// flag telling if Liou's way of computing the update coeff 
  /// imposing positivity has to be used
  bool m_useLiouUpdateCoeff;
  
  /// flag telling to compute fluxes in a decoupled manner
  bool m_useDecoupled;
  
}; // end of class AUSMFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMFlux_hh
