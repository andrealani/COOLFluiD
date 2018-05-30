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
   * Compute dimensional pressure, density and temperature
   */
  virtual void computeDimensionalPressDensTemp(CFreal& pDim, CFreal& rhoDim, CFreal& TDim, CFuint flxIdx);
  
protected:
  
  // temporary vibrational temperature
  RealVector _tempVib;
  
}; // end of class NavierStokesSkinFrictionHeatFRNEQ

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFRNEQ_hh
