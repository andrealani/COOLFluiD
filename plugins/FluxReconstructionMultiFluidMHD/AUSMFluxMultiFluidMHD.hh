#ifndef COOLFluiD_FluxReconstructionMethod_AUSMFluxMultiFluidMHD_hh
#define COOLFluiD_FluxReconstructionMethod_AUSMFluxMultiFluidMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement the Flux Reconstruction space discretization.
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent an AUSM flux
 *
 * This class defines the Module FluxReconstructionMultiFluidMHD
 */

template <class UPDATEVAR>
class AUSMFluxMultiFluid : public RiemannFlux {
public:  // methods

  //Constructor
  
  AUSMFluxMultiFluid(const std::string& name);

  //Destructor
  
  ~AUSMFluxMultiFluid();

  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState,
                                  Framework::State& rState,
                                  const RealVector& normal);

  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState, RealVector& lExtraVars,
                                  Framework::State& rState, RealVector& rExtraVars,
                                  const RealVector& normal);
                                  

  /// Gets the Class name
  static std::string getClassName()
  {
    return "AUSMFluxMultiFluid";
  }
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Set up private data and data
  virtual void setup();
  
  /// Unset up private data and data
  virtual void unsetup();

protected: // helper function
  
  /**
   * Compute the flux
   */
  virtual void computeMassFluxImpl(RealVector& m_rFlux, const RealVector normal);
  
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

/**
   * Compute the Aplus Matrix
   */
  void computeMatrixAplus(const RealVector normal);   
  
   /**
   * Compute the Aminus Matrix
   */
  void computeMatrixAminus(const RealVector normal);  

protected: // data
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> m_updateVarSet;
  
  /// left data  
  RealVector* m_lData;
  
  /// right data
  RealVector* m_rData;

  /// user defined coefficient for the calculation method of a12
  CFuint m_choiceA12;
  
  /// normal velocity left
  RealVector m_unL;
  
  /// normal velocity right
  RealVector m_unR;
  
  /// interface sound speed
  RealVector m_a12;
  
  /// mach number left
  RealVector m_mL;
  
  /// mach number right
  RealVector m_mR;
    
  /// interface mass flux
  RealVector m_mflux12;
  
  /// interface pressure flux
  RealVector m_p12;

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

  /// Coeff for numerical viscosity
  CFreal m_coeff;

  /// Coeff for numerical viscosity in the magnetic field
  CFreal m_Bdiss;

  /// Coeff for numerical viscosity in the electric field
  CFreal m_Ediss; 

  /// flag for using the MacCormack scaling in Maxwell equations 
  /// the default is true
  bool m_useMacCormackScaling;   

  virtual CFreal getMacCormackCoeff(){return m_coeff;}

  virtual CFreal getMagneticDissCoeff(){return m_Bdiss;}
  
  virtual bool getUseMacCormackScaling(){return m_useMacCormackScaling;}

}; // end of class AUSMFluxMultiFluid

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMultiFluidMHD/AUSMFluxMultiFluidMHD.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMultiFluidMHD_hh
