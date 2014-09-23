#ifndef COOLFluiD_Numerics_SpectralFD_BCSubInletEulerRhoV3D_hh
#define COOLFluiD_Numerics_SpectralFD_BCSubInletEulerRhoV3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler3DVarSet;
    }
  }

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a subsonic inlet boundary condition for the 3D Euler/Navier-Stokes equations
 * Mass density and velocity components are imposed
 *
 * @author Kris Van den Abeele
 */
class BCSubInletEulerRhoV3D : public BCStateComputer {

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCSubInletEulerRhoV3D(const std::string& name);

  /// Destructor
  ~BCSubInletEulerRhoV3D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCSubInletEulerRhoV3D";
  }

  /// Set up private data and data
  void setup();

  /**
   * Configures the command.
   */
  void configure ( Config::ConfigArgs& args );

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

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler3DVarSet> m_eulerVarSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of for the variable names
  std::vector<std::string> m_vars;

  /// the VectorialFunction to use to parse the user conditions
  Framework::VectorialFunction m_vFunction;

  /// input variables (rho, vx, vy, vz)
  RealVector m_inputVars;

  /// reference value for mass density for non-dimensionalization
  CFreal m_rhoRef;

  /// reference value for velocity for non-dimensionalization
  CFreal m_velRef;

}; // class BCSubInletEulerRhoV3D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_BCSubInletEulerRhoV3D_hh

