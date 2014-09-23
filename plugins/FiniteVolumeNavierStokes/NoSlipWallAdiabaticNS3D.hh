#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNS3D_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNS3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes { class Euler3DVarSet; }

  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the adiabatic no-slip will bc
   *
   * @author Andrea Lani
   * @author Thomas Wuilbaut
   *
   */
class NoSlipWallAdiabaticNS3D : public FVMCC_BC {

public:

  /**
   * Constructor
   */
  NoSlipWallAdiabaticNS3D(const std::string& name);

  /**
   * Default destructor
   */
  ~NoSlipWallAdiabaticNS3D();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

 private:

  /// physical model var set
  Common::SafePtr<Physics::NavierStokes::Euler3DVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

}; // end of class NoSlipWallAdiabaticNS3D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNS3D_hh
