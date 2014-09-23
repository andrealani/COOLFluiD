#ifndef COOLFluiD_Numerics_FiniteVolume_RoeFluxTurb_hh
#define COOLFluiD_Numerics_FiniteVolume_RoeFluxTurb_hh

//////////////////////////////////////////////////////////////////////////////

#include "RoeFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Roe flux corresponding to the Euler
 * physical model 2D (in conservative variables)
 * The lax-friedrichs flux is used for the turbulent variables
 *
 * @author Thomas Wuilbaut
 *
 */
class RoeFluxTurb : public RoeFlux {
public:
  
  /**
   * Constructor
   */
  RoeFluxTurb(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RoeFluxTurb();

  /**
   * Set up private data
   */
  virtual void setup()
  {
    RoeFlux::setup();
  }

  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
private: // helper function

  /**
   * Compute the flux for jacobian evaluation
   */
  void computePerturbCase(RealVector& result);
  
  /**
   * Compute the flux for RHS evaluation
   */
  void computeUnperturbCase(RealVector& result);
  
protected:
  
  /// diffusion reduction coefficient defined interactively
  CFreal _currentDiffRedCoeff;
  
}; // end of class RoeFluxTurb

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RoeFluxTurb_hh
