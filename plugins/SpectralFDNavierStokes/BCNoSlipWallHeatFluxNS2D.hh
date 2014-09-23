#ifndef COOLFluiD_Numerics_SpectralFD_BCNoSlipWallHeatFluxNS2D_hh
#define COOLFluiD_Numerics_SpectralFD_BCNoSlipWallHeatFluxNS2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a no-slip-wall boundary condition with imposed heat flux
 * for the 2D Euler/Navier-Stokes equations.
 *
 * @author Kris Van den Abeele
 */
class BCNoSlipWallHeatFluxNS2D : public BCStateComputer {

public:  // methods

  /// Constructor
  BCNoSlipWallHeatFluxNS2D(const std::string& name);

  /// Destructor
  ~BCNoSlipWallHeatFluxNS2D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCNoSlipWallHeatFluxNS2D";
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

  /// imposed heat flux
  /// @todo add option for this variable, with imposed function for heat flux as a function of coordinates

}; // class BCNoSlipWallHeatFluxNS2D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_BCNoSlipWallHeatFluxNS2D_hh

