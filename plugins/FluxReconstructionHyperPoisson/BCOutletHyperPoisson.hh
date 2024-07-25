#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCOutletHyperPoisson_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCOutletHyperPoisson_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace HyperPoisson {
      class HyperPoisson3DVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an outlet boundary condition
 * for the 3D Hyperbolized Poisson equations.
 *
 * @author Rayan Dhib
 */
class BCOutletHyperPoisson : public BCStateComputer {

public:  // methods

  /// Constructor
  BCOutletHyperPoisson(const std::string& name);

  /// Destructor
  ~BCOutletHyperPoisson();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCOutletHyperPoisson";
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
  Common::SafePtr<Physics::HyperPoisson::HyperPoisson3DVarSet> m_varSet;

  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// number of flux pnts on a face
  CFuint m_nbrFaceFlxPnts;
  
  CFreal m_refPhi;
  
  /// inner states
  std::vector< RealVector > m_tempStates;

}; // class BCOutletHyperPoisson

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCOutletHyperPoisson_hh

