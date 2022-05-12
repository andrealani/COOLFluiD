#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCSuperOutletProj3DMHD_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCSuperOutletProj3DMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an outlet boundary condition
 * for the 3D MHD equations.
 *
 * @author Ray Vandenhoeck
 */
class BCSuperOutletProj3DMHD : public BCStateComputer {

public:  // methods

  /// Constructor
  BCSuperOutletProj3DMHD(const std::string& name);

  /// Destructor
  ~BCSuperOutletProj3DMHD();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCSuperOutletProj3DMHD";
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

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> m_varSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;
  
  CFreal m_refPhi;
  
  bool m_imposeBDecay;
  
  /// inner states
  std::vector< RealVector > m_tempStates;

}; // class BCSuperOutletProj3DMHD

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCSuperOutletProj3DMHD_hh

