#ifndef COOLFluiD_Numerics_FiniteVolume_HUSFlux2D_hh
#define COOLFluiD_Numerics_FiniteVolume_HUSFlux2D_hh

//////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/OsherFlux2D.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents the HUS flux splitting corresponding
 * to the Euler physical model 2D (in conservative variables)
 *
 * @author Janos Molnar
 * @author Andrea Lani
 & @author Justin Elszasz
 *
 */
template <class UPDATEVAR>   
class HUSFlux2D : public OsherFlux2D<UPDATEVAR> {
public: // classes
  
  /**
   * Constructor
   */
  HUSFlux2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~HUSFlux2D();
  
  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();
  
private: // helper functions
  
  /**
   * Compute the flux assuming full euler equations + species equations
   * in chemical non equilibrium
   */
  virtual void CNEQ(RealVector& result);
  
  /**
   * Compute the flux assuming thermo-chemical non equilibrium
   */
  virtual void TCNEQ(RealVector& result);
  
}; // end of class HUSFlux2D

//////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#include "HUSFlux2D.ci"

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_HUSFlux2D_hh
