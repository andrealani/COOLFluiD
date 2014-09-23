#ifndef COOLFluiD_Numerics_FiniteVolume_VanLeer2D_hh
#define COOLFluiD_Numerics_FiniteVolume_VanLeer2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_FluxSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Van Leer flux splitting corresponding
 * to the Euler physical model 2D (in conservative variables)
 *
 * @author Janos Molnar
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>
class VanLeer2D : public FVMCC_FluxSplitter {
public: // classes
  
    /**
   * Constructor
   */
  VanLeer2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~VanLeer2D();

  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();
  
  /**
   * Compute the flux in the current face
   */
  virtual void compute(RealVector& result);
  
protected:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> _updateVarSet;
  
  // temporary unit normal
  RealVector _tempUnitNormal;
  
  // vectors to calculate the rest of the fluxes according to
  // Van Leer (v is for vector, L: plus, R: minus)
  RealVector _vfluxL;
  RealVector _vfluxR;
  
}; // end of class VanLeer2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "VanLeer2D.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_VanLeer2D_hh
