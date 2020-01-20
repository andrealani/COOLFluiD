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

