#ifndef COOLFluiD_Numerics_SpectralFD_BCSubInletLinEulerRhoV3D_hh
#define COOLFluiD_Numerics_SpectralFD_BCSubInletLinEulerRhoV3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace LinearizedEuler {
      class LinEuler3DVarSet;
    }
  }

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a subsonic inlet boundary condition for the 3D Linearized Euler equations
 * Mass density and velocity components are imposed
 *
 * @author Ghader Ghorbaniasl
 */
class BCSubInletLinEulerRhoV3D : public BCStateComputer {

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCSubInletLinEulerRhoV3D(const std::string& name);

  /// Destructor
  ~BCSubInletLinEulerRhoV3D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCSubInletLinEulerRhoV3D";
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

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();


protected: // data

  /// the socket stores the data of the mean flow
  Framework::DataSocketSink<RealVector> socket_meanflow;

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::LinearizedEuler::LinEuler3DVarSet> m_linEulerVarSet;

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


}; // class BCSubInletLinEulerRhoV3D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_BCFarField1DCharLinEuler3D_hh

