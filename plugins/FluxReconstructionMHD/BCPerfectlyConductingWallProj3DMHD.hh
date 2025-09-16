#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCPerfectlyConductingWallProj3DMHD_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCPerfectlyConductingWallProj3DMHD_hh

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
 * This class represents a perfectly conducting wall boundary condition
 * for the 3D MHD equations.
 *
 * @author Rayan Dhib
 */
class BCPerfectlyConductingWallProj3DMHD : public BCStateComputer {

public:  // methods

  /// Constructor
  BCPerfectlyConductingWallProj3DMHD(const std::string& name);

  /// Destructor
  ~BCPerfectlyConductingWallProj3DMHD();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCPerfectlyConductingWallProj3DMHD";
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

  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// number of flux pnts on a face
  CFuint m_nbrFaceFlxPnts;
  
  /// inner states
  std::vector< RealVector > m_tempStates;

}; // class BCPerfectlyConductingWallProj3DMHD

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_BCPerfectlyConductingWallProj3DMHD_hh
