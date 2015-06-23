#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/AUSMFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the AUSM flux
 *
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>
class AUSMPlusUpFlux : public AUSMFlux<UPDATEVAR> {
public:
  
  /**
   * Constructor
   */
  AUSMPlusUpFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMPlusUpFlux();
  
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
    AUSMFlux<UPDATEVAR>::configure(args);
  }
  
protected:

  /**
   * Compute the interface mass flux
   */
  virtual void computeMassFlux();
  
  /**
   * Compute the interface pressure flux
   */
  virtual void computePressureFlux();
  
  /// Correct the given Mach number
  virtual CFreal correctMachInf(CFreal oldMach) const
  {
    return oldMach;
  }
  
  /// compute flux derivative with respect to pressure
  virtual void dFdP(CFuint side, CFuint iVar, CFreal* row);
  
  /// compute flux derivative with respect to density
  virtual void dFdRho(CFuint side, CFuint iVar, CFreal* row);
  
  /// compute flux derivative with respect to temperature
  virtual void dFdT(CFuint side, CFuint iVar, CFreal* row);
  
  /// compute flux derivative with respect to velocity u
  virtual void dFdU(CFuint side, CFuint iVar, CFreal* row);
  
  /// compute flux derivative with respect to velocity v
  virtual void dFdV(CFuint side, CFuint iVar, CFreal* row);
  
  /// compute flux derivative with respect to velocity w
  virtual void dFdW(CFuint side, CFuint iVar, CFreal* row);
  
  /// compute flux derivative with respect to k
  virtual void dFdK(CFuint side, CFuint iVar, CFreal* row);
  
private:
  
  /// preconditioning coefficient 
  CFreal m_fa;
  
  /// P5 plus coefficient 
  CFreal m_P5Plus;
  
  /// P5 minus coefficient 
  CFreal m_P5Minus;
    
  /// user defined coefficient for Ku
  CFreal m_coeffKu;
  
  /// user defined coefficient for Kp
  CFreal m_coeffKp;
  
  /// user defined coefficient for sigma
  CFreal m_coeffSigma;

  /// mach infinity
  CFreal m_machInf;
  
  /// beta  coefficient
  CFreal m_beta;
  
}; // end of class AUSMPlusUpFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMPlusUpFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFlux_hh
