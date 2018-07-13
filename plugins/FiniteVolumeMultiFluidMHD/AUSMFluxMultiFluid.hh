#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMFluxMultiFluid_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMFluxMultiFluid_hh

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
 * This class computes the AUSM flux for the fluid and 
 * Steger-Warming for Maxwell equations
 *
 * @author Alejandro Alvarez Laguna 
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>                                       
class AUSMFluxMultiFluid : public FVMCC_FluxSplitter {
public:

 
  /**
   * Constructor
   */
  AUSMFluxMultiFluid(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMFluxMultiFluid();
  
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
  virtual void computeLeftJacobian();
  
  /**
   * Compute the right flux jacobian
   */
  virtual void computeRightJacobian();
  
protected:
  
  /**
   * Compute the flux
   */
  virtual void computeMassFluxImpl(RealVector& result);
  
  /**
   * Compute the interface mass flux
   */
  virtual void computeMassFlux() = 0;
  
  /**
   * Compute the interface pressure flux
   */
  virtual void computePressureFlux() = 0;
  
//   /**
//    * Compute the interface Electromagnetic flux
//    */
//   virtual void computeMaxwellFlux() = 0;  
  
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
  
   /**
   * Compute the Aplus Matrix
   */
  void computeMatrixAplus();   
  
   /**
   * Compute the Aminus Matrix
   */
  void computeMatrixAminus();  
  
protected:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> m_updateVarSet;
  
  /// left data  
  RealVector* m_lData;
  
  /// right data
  RealVector* m_rData;
  
  /// interface sound speed
  RealVector m_a12Vec;
  
  /// mach number left
  RealVector m_mL;
  
  /// mach number right
  RealVector m_mR;
      
  /// normal velocity left
  RealVector m_unL;
  
  /// normal velocity right
  RealVector m_unR;
    
  /// interface mass flux
  RealVector m_mflux12;
  
  /// interface pressure flux
  RealVector m_p12;
  
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
  
  /// vector with the electromagnetic field variables (LEFT)
  RealVector _EMField_l;

  /// vector with the electromagnetic field variables (RIGHT)
  RealVector _EMField_r;
  
  /// A plus Matrix
  RealMatrix   _Aplus;
  
  /// A minus Matrix
  RealMatrix   _Aminus; 
  
  /// A minus Matrix
  RealVector   _resultEM;   
  
  /// Array with the convected quatities of the convective flux of the species from the left
  RealVector _psi_l;
  
  /// Array with the convected quatities of the convective flux of the species from the right 
  RealVector _psi_r;  
  
  /// Array with the molar mass of the species
  RealVector _m_i;
    
  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;
  
   /// user defined coefficient for the calculation method of a12
  CFuint m_choiceA12;

  /// Coeff for numerical viscosity
  CFreal m_coeff;

  /// Coeff for numerical viscosity in the magnetic field
  CFreal m_Bdiss;

  /// Coeff for numerical viscosity in the magnetic field
  CFreal m_Ediss; 
 
  /// flag telling if Liou's way of computing the update coeff 
  /// imposing positivity has to be used
  bool m_useLiouUpdateCoeff; 
  
  /// flag for using the MacCormack scaling in Maxwell equations 
  /// the default is true
  bool m_useMacCormackScaling;   

  virtual CFreal getMacCormackCoeff(){return m_coeff;}

  virtual CFreal getMagneticDissCoeff(){return m_Bdiss;}

  virtual CFreal getElectricDissCoeff(){return m_Ediss;}
  
  virtual bool getUseMacCormackScaling(){return m_useMacCormackScaling;}

}; // end of class AUSMFluxMultiFluid


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMFluxMultiFluid.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMFluxMultiFluid_hh
