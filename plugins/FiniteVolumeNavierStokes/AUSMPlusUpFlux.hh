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
  
private:
  
  /// preconditioning coefficient 
  CFreal m_fa;
  
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
