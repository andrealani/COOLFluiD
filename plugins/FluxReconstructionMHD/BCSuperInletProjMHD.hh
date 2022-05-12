#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCSuperInletProjMHD_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCSuperInletProjMHD_hh

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
 * This class represents an inlet boundary condition
 * for the 3D MHD equations.
 *
 * @author Ray Vandenhoeck
 */
class BCSuperInletProjMHD : public BCStateComputer {

public:  // methods

  /// Constructor
  BCSuperInletProjMHD(const std::string& name);

  /// Destructor
  ~BCSuperInletProjMHD();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCSuperInletProjMHD";
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
  
  /// map TRS name -> initial solution array that will be used as BC value
  Common::CFMap<std::string, RealVector*> m_initialSolutionMap;
  
  ///bnd density value
  CFreal m_rhoBC;
  
  ///bnd p value
  CFreal m_pBC;

}; // class BCSuperInletProjMHD

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCSuperInletProjMHD_hh

