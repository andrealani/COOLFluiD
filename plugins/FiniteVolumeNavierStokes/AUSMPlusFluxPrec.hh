#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMPlusFluxPrec_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMPlusFluxPrec_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/AUSMFluxPrec.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the AUSM+ flux with preconditioning
 *
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>
class AUSMPlusFluxPrec : public AUSMFluxPrec<UPDATEVAR> {
public:
  
  /**
   * Constructor
   */
  AUSMPlusFluxPrec(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMPlusFluxPrec();
  
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
    AUSMFluxPrec<UPDATEVAR>::configure(args);
  }

protected:

  /**
   * Compute the interface mass flux & the pressure diffusion interface mass flux of the preconditionned AUSM scheme
   */
  virtual void computeMassFlux();
  
  /**
   * Compute the interface pressure flux
   */
  virtual void computePressureFlux();
  
private:
  
  /// beta  coefficient
  CFreal m_beta;
  
  /// alpha  coefficient
  CFreal m_alpha;

  //// Kinematic Viscosity used to preconditioned
  CFreal m_nu;

}; // end of class AUSMPlusFluxPrec

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMPlusFluxPrec.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMPlusFluxPrec_hh
