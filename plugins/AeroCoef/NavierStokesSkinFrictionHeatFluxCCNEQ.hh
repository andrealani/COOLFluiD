#ifndef COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxCCNEQ_hh
#define COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxCCNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "AeroCoef/NavierStokesSkinFrictionHeatFluxCC.hh"

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
 *
 */

class NavierStokesSkinFrictionHeatFluxCCNEQ : public NavierStokesSkinFrictionHeatFluxCC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokesSkinFrictionHeatFluxCCNEQ(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesSkinFrictionHeatFluxCCNEQ();
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected:
  
  /**
   * Compute dimensional pressure, density and temperature
   */
  virtual void computeDimensionalPressDensTemp(CFreal& pDim, CFreal& rhoDim, CFreal& TDim);
  
protected:
  
  // temporary vibrational temperature
  RealVector _tempVib;
  
}; // end of class NavierStokesSkinFrictionHeatFluxCCNEQ

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxCCNEQ_hh
