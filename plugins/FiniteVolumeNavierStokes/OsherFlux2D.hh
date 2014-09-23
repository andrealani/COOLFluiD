#ifndef COOLFluiD_Numerics_FiniteVolume_OsherFlux2D_hh
#define COOLFluiD_Numerics_FiniteVolume_OsherFlux2D_hh

//////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_FluxSplitter.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents the Osher flux splitting corresponding
 * to the Euler physical model 2D (in conservative variables)
 
 * @author Justin Elszasz
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>   
class OsherFlux2D : public FVMCC_FluxSplitter {
public: // classes

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  OsherFlux2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~OsherFlux2D();
  
  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();
  
  /**
   * Compute the flux in the current face
   */
  virtual void compute(RealVector& result);
  
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
  
protected:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> _updateVarSet;
  
  /// user defined flar for computation order
  bool _isNatural;
  
  /// coefficient for computation order
  CFreal _beta; 
  
  // temporary unit normal
  RealVector _tempUnitNormal;
  
  // vectors to calculate the rest of the fluxes according to
    // Van Leer (v is for vector, L: plus, R: minus)
  RealVector _vfluxL;
  RealVector _vfluxR;
  
  // vectors for the species and vibrational energy state data
  RealVector _YsLS;
  RealVector _YsRS;
  RealVector _evLS;
  RealVector _evRS;
  
}; // end of class OsherFlux2D

//////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#include "OsherFlux2D.ci"

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_OsherFlux2D_hh
