#ifndef COOLFluiD_Numerics_FiniteElement_FEM_VolumeIntegrator_hh
#define COOLFluiD_Numerics_FiniteElement_FEM_VolumeIntegrator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VolumeIntegrator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class provides the functionality of a
 * GeometricEntity volume integrator.
 *
 * @author Thomas Wuilbaut
 *
 */
class FEM_VolumeIntegrator : public Framework::VolumeIntegrator {
public:

  /**
   * Constructor
   */
  FEM_VolumeIntegrator();

  /**
   * Destructor
   */
  ~FEM_VolumeIntegrator();

  /**
   * Set up the integrator for the simulation
   */
  void setup()
  {
    Framework::VolumeIntegrator::setup();
  }

  /**
   * PreCompute all the data linked to the cell.
   */
  void precomputeCellData(Framework::GeometricEntity* geo);

  /**
   * Compute a volume integration for on a GeometricEntity
   * where the functor is a FEMEntity and filling the LocalElementData.
   */
  template <class FUNCTOR, class RETURN>
  void integrateGeneralFEMEntityOnGeoEnt(FUNCTOR& functor,
                                         RETURN& result);

  /**
   * Compute a volume integration for on a GeometricEntity
   * where the functor is a FEMEntity and filling the LocalElementData.
   */
  template <class FUNCTOR, class RETURN>
  void integrateFaceFEMEntityOnGeoEnt(FUNCTOR& functor,
                                         RETURN& result);

  /**
   * Compute a volume integration for on a GeometricEntity
   * where the functor is a FEMEntity and filling the LocalElementData.
   */
  template <class FUNCTOR, class RETURN>
  void integrateFastGeneralFEMEntityOnGeoEnt(FUNCTOR& functor,
                                             RETURN& result);

private:

  /**
   * Copy constructor
   */
  FEM_VolumeIntegrator(const FEM_VolumeIntegrator&);

  /**
   * assignment operator
   */
  FEM_VolumeIntegrator& operator= (const FEM_VolumeIntegrator&);

private: //data

  CFuint _lastPrecomputedCellID;

}; // end class FEM_VolumeIntegrator

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FEM_VolumeIntegrator.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_FEM_VolumeIntegrator_hh
