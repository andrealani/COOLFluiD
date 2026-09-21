#ifndef COOLFluiD_FluxReconstructionMethod_BCSubInletTurbUVT2D_hh
#define COOLFluiD_FluxReconstructionMethod_BCSubInletTurbUVT2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "NavierStokes/Euler2DPuvt.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

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
 * This class represents a subsonic inlet boundary condition for the 2D RANS equations
 *
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class BCSubInletTurbUVT2D : public BCStateComputer {

public:  // methods
    
  typedef Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DPuvt> ConvTurb2DVarSet;

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCSubInletTurbUVT2D(const std::string& name);

  /// Destructor
  ~BCSubInletTurbUVT2D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCSubInletTurbUVT2D";
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
   * the gradient variables extrapolated to the flux points with, from the ghost
   * state, the prescribed velocity, temperature and turbulence variables. The
   * gradient variables are the update variables, Puvt plus the turbulence
   * variables.
   */
  void computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                          const std::vector< Framework::State* >& intStates,
                          const std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& unitNormals,
                          const std::vector< RealVector >& flxPntCoords,
                          std::vector< RealVector* >& bndGradVars);

  /**
   * Sets the boundary states U_b the diffusive boundary flux and the transport
   * properties are evaluated at: the ghost state: interior pressure, prescribed velocity, temperature and turbulence variables.
   */
  void computeBndStates(const std::vector< Framework::State* >& intStates,
                        const std::vector< Framework::State* >& ghostStates,
                        const std::vector< RealVector >& unitNormals,
                        const std::vector< RealVector >& flxPntCoords,
                        std::vector< RealVector* >& bndStates);

  /**
   * Sets the boundary gradients q_b of the diffusive boundary flux: the compact face gradients, q_b = q.
   */
  void computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                       std::vector< std::vector< RealVector* > >& bndGrads,
                       const std::vector< RealVector* >& bndStates,
                       const std::vector< RealVector >& unitNormals,
                       const std::vector< RealVector >& flxPntCoords);

protected: // data
    
  /// physical model var set
  Common::SafePtr<ConvTurb2DVarSet> m_varSetTurb;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// x velocity
  CFreal m_uinf;

  /// y velocity
  CFreal m_vinf;

  /// static temperature
  CFreal m_temperature;

  /// Turbulent Variables K, Omega...
  std::vector<CFreal> m_turbVars;

}; // class BCSubInletTurbUVT2D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCSubInletTurbUVT2D_hh

