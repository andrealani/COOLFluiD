#ifndef COOLFluiD_FluxReconstructionMethod_AUSMFlux_hh
#define COOLFluiD_FluxReconstructionMethod_AUSMFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent an AUSM flux
 *
 * @author Ray Vandenhoeck
 * @author Alexander Papen
 */

template <class UPDATEVAR>
class AUSMFlux : public RiemannFlux {
  
  typedef UPDATEVAR UPVAR;

public:  // methods

  /// Constructor
  AUSMFlux(const std::string& name);

  /// Destructor
  ~AUSMFlux();

  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState, RealVector& lExtraVars,
                                  Framework::State& rState, RealVector& rExtraVars,
                                  const RealVector& normal);
                                  
  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState,
                                  Framework::State& rState,
                                  const RealVector& normal);


  /// Gets the Class name
  static std::string getClassName()
  {
    return "AUSMFlux";
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
  virtual void computeMassFluxImpl(RealVector& result, const RealVector normal);
  
  /**
   * Compute the interface mass flux
   */
  virtual void computeMassFlux() = 0;
  
  /**
   * Compute the interface pressure flux
   */
  virtual void computePressureFlux() = 0;

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

protected: // data
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> m_updateVarSet;
  
  /// user defined coefficient for the calculation method of a12
  CFuint m_choiceA12;
  
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

}; // class AUSMFlux

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/AUSMFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_AUSMFlux_hh
