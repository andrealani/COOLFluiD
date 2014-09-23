#ifndef COOLFluiD_Numerics_SpectralFV_BCNoSlipWallHeatFluxNS3D_hh
#define COOLFluiD_Numerics_SpectralFV_BCNoSlipWallHeatFluxNS3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFV/BCStateComputer.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler3DVarSet;
    }
  }

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a no-slip-wall boundary condition with imposed heat flux
 * for the 3D Euler/Navier-Stokes equations.
 *
 * @author Kris Van den Abeele
 */
class BCNoSlipWallHeatFluxNS3D : public BCStateComputer {

public:  // methods

  /// Constructor
  BCNoSlipWallHeatFluxNS3D(const std::string& name);

  /// Destructor
  ~BCNoSlipWallHeatFluxNS3D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCNoSlipWallHeatFluxNS3D";
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
  Common::SafePtr<Physics::NavierStokes::Euler3DVarSet> m_eulerVarSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// imposed heat flux
  /// @todo add option for this variable, with imposed function for heat flux as a function of coordinates

}; // class BCNoSlipWallHeatFluxNS3D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_BCNoSlipWallHeatFluxNS3D_hh

