#ifndef COOLFluiD_Numerics_SpectralFV_BCMirrorEuler2D_hh
#define COOLFluiD_Numerics_SpectralFV_BCMirrorEuler2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFV/BCStateComputer.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a mirror boundary condition for the 2D Euler/Navier-Stokes equations
 *
 * @author Kris Van den Abeele
 */
class BCMirrorEuler2D : public BCStateComputer {

public:  // methods

  /// Constructor
  BCMirrorEuler2D(const std::string& name);

  /// Destructor
  ~BCMirrorEuler2D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCMirrorEuler2D";
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

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> m_eulerVarSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

}; // class BCMirrorEuler2D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_BCMirrorEuler2D_hh

