#ifndef COOLFluiD_Numerics_FiniteVolume_UnsteadyNoSlipWallAdiabaticNSTurb2D_hh
#define COOLFluiD_Numerics_FiniteVolume_UnsteadyNoSlipWallAdiabaticNSTurb2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/Node.hh"
#include "NavierStokes/Euler2DPuvt.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
      class NavierStokes2DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the UnsteadySlipWall bc
   *
   * @author Thomas Wuilbaut
   *
   */
class UnsteadyNoSlipWallAdiabaticNSTurb2D : public FVMCC_BC {

public:

  typedef Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DPuvt> ConvTurb2DVarSet;
 
  typedef Physics::NavierStokes::NavierStokesTurbVarSet<
    Physics::NavierStokes::NavierStokes2DVarSet, 0> DiffTurb2DVarSet;

  /**
   * Constructor
   */
  UnsteadyNoSlipWallAdiabaticNSTurb2D(const std::string& name);

  /**
   * Default destructor
   */
  ~UnsteadyNoSlipWallAdiabaticNSTurb2D();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private://function

  /**
   * Compute the mesh speed at the ghost state
   */
  void computeGhostStateSpeed(Framework::GeometricEntity *const face);

 private://data

  /// physical model convective variable set
  Common::SafePtr<ConvTurb2DVarSet> _varSetTurb;

  /// physical model diffusive variable set
  Common::SafePtr<DiffTurb2DVarSet> _diffVarTurb;
  
  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// handle to the nodes at their past position
  Framework::DataSocketSink< Framework::Node*> socket_pastNodes;

  /// handle to the nodes at their future position
  Framework::DataSocketSink< Framework::Node*> socket_futureNodes;

  /// Temporary Storage of coordinates
  RealVector _coord;

  /// Speed of the mesh at the face
  RealVector _speed;

  /// handle to the distances from the states to the walls
  Framework::DataSocketSink< CFreal> socket_wallDistance;

}; // end of class UnsteadyNoSlipWallAdiabaticNSTurb2D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_UnsteadyNoSlipWallAdiabaticNSTurb2D_hh
