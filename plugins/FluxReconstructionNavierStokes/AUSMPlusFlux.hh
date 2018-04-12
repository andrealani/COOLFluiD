#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_AUSMPlusFlux_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_AUSMPlusFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/AUSMFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the AUSM+ flux
 *
 * @author Ray Vandenhoeck
 * @author Alexander Papen
 *
 */
template <class UPDATEVAR>
class AUSMPlusFlux : public AUSMFlux<UPDATEVAR> {
public:
  
  /**
   * Constructor
   */
  AUSMPlusFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMPlusFlux();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  

protected:

  /**
   * Compute the interface mass flux
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

}; // end of class AUSMPlusFlux

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/AUSMPlusFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_AUSMPlusFlux_hh
