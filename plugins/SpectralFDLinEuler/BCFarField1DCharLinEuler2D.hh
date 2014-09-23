#ifndef COOLFluiD_Numerics_SpectralFD_BCFarField1DCharLinEuler2D_hh
#define COOLFluiD_Numerics_SpectralFD_BCFarField1DCharLinEuler2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace LinearizedEuler {
      class LinEuler2DVarSet;
    }
  }

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a simple far field boundary condition for the 2D linearized Euler equations
 * It is based on the characteristics of the 1D linearized Euler equations. They are not sufficient
 * for real linearized Euler computations, because they still cause quite a lot of reflections. Should
 * be combined with a sponge zone.
 *
 * @author Kris Van den Abeele
 */
class BCFarField1DCharLinEuler2D : public BCStateComputer {

public:  // methods

  /// Constructor
  BCFarField1DCharLinEuler2D(const std::string& name);

  /// Destructor
  ~BCFarField1DCharLinEuler2D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCFarField1DCharLinEuler2D";
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
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();


protected: // data

  /// the socket stores the data of the mean flow
  Framework::DataSocketSink<RealVector> socket_meanflow;

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::LinearizedEuler::LinEuler2DVarSet> m_linEulerVarSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// tangential velocity component
//   RealVector m_tanVel;

}; // class BCFarField1DCharLinEuler2D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_BCFarField1DCharLinEuler2D_hh

