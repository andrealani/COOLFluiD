#ifndef COOLFluiD_FluxReconstructionMethod_BCNoSlipWallTurb3D_hh
#define COOLFluiD_FluxReconstructionMethod_BCNoSlipWallTurb3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/Euler3DPvt.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler3DVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a no-slip wall boundary condition for the 3D RANS equations
 *
 * Ghost state of the mean flow, with u_i the interior state at the wall flux point (the trace of the
 * cell polynomial): interior density, velocity reflected about the wall velocity, temperature
 * reflected about the wall temperature on an isothermal wall (T_g = 2 T_wall - T_i, at least
 * 0.01 T_wall) or copied on a heat-flux wall, pressure from density and temperature. The wall model of
 * the turbulence variables uses the wall values (interior pressure, wall or interior temperature).
 * LegacyGhost = true restores the previous ghost (the wall state).
 *
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class BCNoSlipWallTurb3D : public BCStateComputer {

public:  // methods
    
  typedef Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler3DPvt<Physics::NavierStokes::Euler3DVarSet> > ConvTurb3DVarSet;
  typedef Physics::NavierStokes::NavierStokesTurbVarSet<Physics::NavierStokes::NavierStokes3DVarSet, 0> DiffTurb3DVarSet;


  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCNoSlipWallTurb3D(const std::string& name);

  /// Destructor
  ~BCNoSlipWallTurb3D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCNoSlipWallTurb3D";
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

  /**
   * Sets the boundary values g_b of the gradient variables at the flux points:
   * the gradient variables extrapolated to the flux points with the wall
   * velocity, the wall temperature of an isothermal wall, and KWall and the
   * wall value of log-omega the ghost state carries. The gradient variables are
   * the update variables, Puvt plus the turbulence variables.
   */
  void computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                          const std::vector< Framework::State* >& intStates,
                          const std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& unitNormals,
                          const std::vector< RealVector >& flxPntCoords,
                          std::vector< RealVector* >& bndGradVars);

  /**
   * Sets the boundary states U_b the diffusive boundary flux and the transport
   * properties are evaluated at: the wall state, that is the interior pressure,
   * the wall velocity, the wall temperature of an isothermal wall or the
   * interior one of a heat-flux wall, and the wall values of the turbulence
   * variables the ghost state carries.
   */
  void computeBndStates(const std::vector< Framework::State* >& intStates,
                        const std::vector< Framework::State* >& ghostStates,
                        const std::vector< RealVector >& unitNormals,
                        const std::vector< RealVector >& flxPntCoords,
                        std::vector< RealVector* >& bndStates);

  /**
   * Sets the boundary gradients q_b of the diffusive boundary flux: the compact face gradients, q_b = q, without their normal component for
   * gamma and Re_theta of the four-equation models (no diffusive flux of the
   * transition variables through the wall).
   */
  void computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                       std::vector< std::vector< RealVector* > >& bndGrads,
                       const std::vector< RealVector* >& bndStates,
                       const std::vector< RealVector >& unitNormals,
                       const std::vector< RealVector >& flxPntCoords);

  /**
   * Sets the normal temperature gradient from the prescribed wall heat flux q,
   * see prescribeNSWallHeatFlux; q = 0 gives an adiabatic wall. Only when the
   * heat flux is prescribed (HeatFlux true).
   */
  void constrainBndGrads(const RealVector& bndState,
                         std::vector< RealVector* >& bndGrads,
                         const RealVector& unitNormal,
                         const RealVector& flxPntCoord);

protected: // data
    
  /// conv physical model var set
  Common::SafePtr<ConvTurb3DVarSet> m_varSetTurb;
  
  /// diff physical model var set
  Common::SafePtr<DiffTurb3DVarSet> m_diffVarTurb;

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

  /// use the previous ghost state instead of the reflected one
  bool m_legacyGhost;

  /// X-component of a velocity vector of the wall
  CFreal m_xWallVelocity;

  /// Y-component of a velocity vector of the wall
  CFreal m_yWallVelocity;
  
  /// Z-component of a velocity vector of the wall
  CFreal m_zWallVelocity;
  
  /// turb intensity at the wall
  CFreal m_wallK;
  
  /// distance of the first sol pnt to the wall
  CFreal m_wallDist;
  
  /// omega wall multiplication factor
  CFreal m_omegaWallFactor;
  
  /// iteration at which to impose theoretical omegaWall
  CFuint m_imposeOmegaWallIter;
  
  /// previous iteration
  CFuint m_prevIter;

}; // class BCNoSlipWallTurb3D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallTurb3D_hh

