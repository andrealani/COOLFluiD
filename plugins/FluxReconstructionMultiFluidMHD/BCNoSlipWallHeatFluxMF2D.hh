#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallHeatFluxMF2D_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallHeatFluxMF2D_hh

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

class BCNoSlipWallHeatFluxMF2D : public BCStateComputer {

public:  // methods

  /// Constructor
  BCNoSlipWallHeatFluxMF2D(const std::string& name);

  /// Destructor
  ~BCNoSlipWallHeatFluxMF2D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCNoSlipWallHeatFluxMF2D";
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

  /// temperatures at the wall
  std::vector<CFreal> _wallTemp;
  
  /// Electrical conductivity
  CFreal _Ez0;

  /// boolean telling if the wall has constant heat flux
  bool m_heatFlux;
  
  bool m_needsSpatCoord;

  /// wall static temperature
  CFreal m_wallT;  //std::vector<CFreal>
  
  /// wall heat flux
  CFreal m_wallQ;  //std::vector<CFreal>
  
  /// iteration after which is changed to an isothermal wall BC
  CFuint m_changeToIsoT;

}; // class BCNoSlipWallHeatFluxMF2D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallHeatFluxMF2D_hh

