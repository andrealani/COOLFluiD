#ifndef COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFRNEQ_hh
#define COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFRNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "AeroCoef/NavierStokesSkinFrictionHeatFluxFR.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace AeroCoef {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the skin friction and the heat flux for NavierStokes
 * simulations with
 * @see CellCenterFVM
 *
 * @author Andrea Lani
 * @author Rayan Dhib
 *
 */

class NavierStokesSkinFrictionHeatFRNEQ : public NavierStokesSkinFrictionHeatFluxFR {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokesSkinFrictionHeatFRNEQ(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesSkinFrictionHeatFRNEQ();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected:

  /**
   * Computes the Stanton number of the heat flux into the wall:
   *   StantonNumberID 0: St = q/(rho_inf u_inf^3)
   *   StantonNumberID 1: St = q/(rho_inf u_inf (H_inf - h_w)),
   * with H_inf = TotalEnthalpyInf and h_w the static enthalpy at the wall
   * @param heatFlux     heat flux into the wall q
   * @param temperature  wall temperature
   * @param flxIdx       index of the flux point
   */
  CFreal computeStantonNumber(CFreal heatFlux, CFreal temperature, CFuint flxIdx);

  /// Tells whether the states hold transition model variables: never for NEQ
  bool hasTransitionLayout() const
  {
    return false;
  }

  /**
   * Compute dimensional pressure, density and temperature
   */
  virtual void computeDimensionalPressDensTemp(CFreal& pDim, CFreal& rhoDim, CFreal& TDim, CFuint flxIdx);
  
protected:
  
  // temporary vibrational temperature
  RealVector _tempVib;

  /// freestream total enthalpy [J/kg], with chemistry and internal modes
  CFreal m_totalEnthalpyInf;
  
}; // end of class NavierStokesSkinFrictionHeatFRNEQ

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFRNEQ_hh
