#ifndef COOLFluiD_FluxReconstructionMethod_AUSMPlusUpFlux_hh
#define COOLFluiD_FluxReconstructionMethod_AUSMPlusUpFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/AUSMFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the AUSM plus up flux
 *
 * @author Ray Vandenhoeck
 * @author Alexander Papen
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

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/AUSMPlusUpFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFlux_hh
