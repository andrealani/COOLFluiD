#ifndef COOLFluiD_Numerics_FiniteVolume_SubOutletNSTurb3D_hh
#define COOLFluiD_Numerics_FiniteVolume_SubOutletNSTurb3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "NavierStokes/Euler3DPvt.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a subsonic outlet command
 *
 * @author Thomas Wuilbaut
 * @author Milan Zaloudek
 *
 */
class SubOutletNSTurb3D : public FVMCC_BC {
public:

  typedef Physics::NavierStokes::MultiScalarVarSet<
  Physics::NavierStokes::Euler3DPvt<Physics::NavierStokes::Euler3DVarSet> > ConvTurb3DVarSet;
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubOutletNSTurb3D(const std::string& name);

  /**
   * Default destructor
   */
  ~SubOutletNSTurb3D();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

private: // data

  /// physical model (in conservative variables)
  Common::SafePtr<ConvTurb3DVarSet> _varSetTurb;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// static pressure
  CFreal _pressure;

}; // end of class SubOutletEuler3D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubOutletNSTurb3D_hh
