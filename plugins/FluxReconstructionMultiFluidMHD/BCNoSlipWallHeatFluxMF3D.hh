#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallHeatFluxMF3D_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallHeatFluxMF3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

class BCNoSlipWallHeatFluxMF3D : public BCStateComputer {

public:  // methods

  /// Constructor
  BCNoSlipWallHeatFluxMF3D(const std::string& name);

  /// Destructor
  ~BCNoSlipWallHeatFluxMF3D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCNoSlipWallHeatFluxMF3D";
  }
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

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
  Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet<Physics::Maxwell::Maxwell2DProjectionVarSet> > m_varSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// Electrical conductivity
  CFreal _Ez0;

  /// boolean telling if the wall has constant heat flux
  bool m_heatFlux;
  
  bool m_needsSpatCoord;

  /// temperatures at the wall
  std::vector<CFreal> _wallTemp;

  /// wall static temperature
  /*std::vector<CFreal>*/ CFreal m_wallT;
  
  /// wall heat flux
  /*std::vector<CFreal>*/ CFreal m_wallQ;
  
  /// iteration after which is changed to an isothermal wall BC
  CFuint m_changeToIsoT;

}; // class BCNoSlipWallHeatFluxMF3D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallHeatFluxMF3D_hh


