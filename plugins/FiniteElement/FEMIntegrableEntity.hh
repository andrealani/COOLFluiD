#ifndef COOLFluiD_Numerics_FiniteElement_FEMIntegrableEntity_hh
#define COOLFluiD_Numerics_FiniteElement_FEMIntegrableEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElementMethodData.hh"
#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the basic interface for a physics-dependent
 * integrable entity
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
class FEMIntegrableEntity : public FiniteElementMethodStrategy {

public:

  /// Type for the provider of this abstract class
  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,FEMIntegrableEntity> PROVIDER;

public:

  /**
   * Constructor
   */
  FEMIntegrableEntity(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~FEMIntegrableEntity();

  /**
   * Set up private data and data
   */
  virtual void setup();

  /**
   * Compute a contour integration of a physical quantity
   * in the given cell
   */
  virtual void integrateContourInGeo(Framework::GeometricEntity* const geo,
                                     RealVector& result);

  /**
   * Compute a volume integration of a physical quantity
   * in the given cell
   */
  virtual void integrateVolumeInGeo(Framework::GeometricEntity* const geo,
                                    RealVector& result);

  /**
   * Compute a contour integration of a physical quantity
   * in the given cell
   */
  virtual void integrateContourInGeo(Framework::GeometricEntity* const geo,
                                     RealMatrix& result);

  /**
   * Compute a volume integration of a physical quantity
   * in the given cell
   */
  virtual void integrateVolumeInGeo(Framework::GeometricEntity* const geo,
                                    RealMatrix& result);

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "FEMIntegrableEntity";
  }

protected:

  Common::SafePtr<LocalElementData> _localElemData;

}; // end of class FEMIntegrableEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_FEMIntegrableEntity_hh
