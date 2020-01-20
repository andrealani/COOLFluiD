#ifndef COOLFluiD_FluxReconstructionMethod_BCNoSlipWallTurb2D_hh
#define COOLFluiD_FluxReconstructionMethod_BCNoSlipWallTurb2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "NavierStokes/Euler2DPuvt.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a no-slip wall boundary condition for the 2D RANS equations
 *
 * @author Ray Vandenhoeck
 */
class BCNoSlipWallTurb2D : public BCStateComputer {

public:  // methods
    
  typedef Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DPuvt> ConvTurb2DVarSet;
  typedef Physics::NavierStokes::NavierStokesTurbVarSet<Physics::NavierStokes::NavierStokes2DVarSet, 0> DiffTurb2DVarSet;


  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCNoSlipWallTurb2D(const std::string& name);

  /// Destructor
  ~BCNoSlipWallTurb2D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCNoSlipWallTurb2D";
  }

  /// Set up private data and data
  void setup();

  /**
   * Sets the ghost states in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostStates(const std::vector< Framework::State* >& intStates,
                          std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& normals,
                          const std::vector< RealVector >& coords);

  /**
   * Sets the ghost gradients in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                             std::vector< std::vector< RealVector* > >& ghostGrads,
                             const std::vector< RealVector >& normals,
                             const std::vector< RealVector >& coords);

protected: // data
    
  /// conv physical model var set
  Common::SafePtr<ConvTurb2DVarSet> m_varSetTurb;
  
  /// diff physical model var set
  Common::SafePtr<DiffTurb2DVarSet> m_diffVarTurb;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;
  
  /// wall static temperature
  CFreal m_wallT;
  
  /// wall heat flux
  CFreal m_wallQ;
  
  /// boolean telling if the wall has constant heat flux
  bool m_heatFlux;
  
  /// iteration after which is changed to an isothermal wall BC
  CFuint m_changeToIsoT;

  /// X-component of a velocity vector of the wall
  CFreal m_xWallVelocity;

  /// Y-component of a velocity vector of the wall
  CFreal m_yWallVelocity;
  
  /// turb intensity at the wall
  CFreal m_wallK;

}; // class BCNoSlipWallTurb2D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallTurb2D_hh

