#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCSubInletTtPtAlphaEIWRhoiViTi3D_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCSubInletTtPtAlphaEIWRhoiViTi3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "Framework/VectorialFunction.hh"
#include "Common/BadValueException.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MultiFluidMHD {      
      class DiffMFMHD2DVarSet;
      class DiffMFMHD3DVarSet;
    }
    
    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
    namespace Maxwell {
      class Maxwell2DProjectionVarSet;
      class Maxwell3DVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Subsonic Inlet imposing the Velocity and Temperature
   * Maxwell Equations: Perfectly Conducting Wall Condition
   */
class BCSubInletTtPtAlphaEIWRhoiViTi3D : public BCStateComputer {

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCSubInletTtPtAlphaEIWRhoiViTi3D(const std::string& name);

  /// Destructor
  ~BCSubInletTtPtAlphaEIWRhoiViTi3D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCSubInletTtPtAlphaEIWRhoiViTi3D";
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
  Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet<Physics::Maxwell::Maxwell3DVarSet> > m_varSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// total temperature
  std::vector<CFreal>     m_tTotal;

  /// total pressure
  std::vector<CFreal>     m_pTotal;

  /// alpha (for v/u (Y/X))
  std::vector<CFreal>     m_alphaXY;

  /// alpha (for w/u (Z/X))
  std::vector<CFreal>     m_alphaXZ;

}; // class BCSubInletTtPtAlphaEIWRhoiViTi3D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCSubInletTtPtAlphaEIWRhoiViTi3D_hh

