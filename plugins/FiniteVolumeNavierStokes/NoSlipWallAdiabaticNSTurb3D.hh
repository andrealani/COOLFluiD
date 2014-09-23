#ifndef COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNSTurb3D_hh
#define COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNSTurb3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "NavierStokes/Euler3DPvt.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class NavierStokes3DVarSet;
    }
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the no slip wall bc
   * for K-Omega and SA turbulence models
   *
   * @author Thomas Wuilbaut
   * @author Milan Zaloudek
   *
   */
class NoSlipWallAdiabaticNSTurb3D : public FVMCC_BC {

public:

  typedef Physics::NavierStokes::MultiScalarVarSet<
  Physics::NavierStokes::Euler3DPvt<Physics::NavierStokes::Euler3DVarSet> > ConvTurb3DVarSet;
    
  typedef Physics::NavierStokes::NavierStokesTurbVarSet<
    Physics::NavierStokes::NavierStokes3DVarSet, 0> DiffTurb3DVarSet;
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  NoSlipWallAdiabaticNSTurb3D(const std::string& name);

  /**
   * Default destructor
   */
  ~NoSlipWallAdiabaticNSTurb3D();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

 private:

  /// physical model convective variable set
  Common::SafePtr<ConvTurb3DVarSet> _varSetTurb;

  /// physical model diffusive variable set
  Common::SafePtr<DiffTurb3DVarSet> _diffVarTurb;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// X-component of a velocity vector of the wall
  CFreal _xWallVelocity;

  /// Y-component of a velocity vector of the wall
  CFreal _yWallVelocity;

  /// Z-component of a velocity vector of the wall
  CFreal _zWallVelocity;

  /// handle to the distances from the states to the walls
  Framework::DataSocketSink< CFreal> socket_wallDistance;


}; // end of class NoSlipWallAdiabaticNSTurb3D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_NoSlipWallAdiabaticNSTurb3D_hh
