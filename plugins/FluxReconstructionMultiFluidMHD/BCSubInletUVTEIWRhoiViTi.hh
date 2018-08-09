#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCSubInletUVTEIWRhoiViTi_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCSubInletUVTEIWRhoiViTi_hh

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
    }
    
    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
    namespace Maxwell {
      class Maxwell2DProjectionVarSet;
    }
  }
    
  namespace FluxReconstructionMethod {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Subsonic Inlet imposing the Velocity and Temperature
   * Maxwell Equations: Perfectly Conducting Wall Condition
   * 
   * @author Alejandro Alvarez
   *
   */

class BCSubInletUVTEIWRhoiViTi : public BCStateComputer { 

public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  BCSubInletUVTEIWRhoiViTi(const std::string& name);
  
  /**
   * Default destructor
   */
  ~BCSubInletUVTEIWRhoiViTi();
  
  /**
   * Configures the command.
   */
  void configure ( Config::ConfigArgs& args );
  
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
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

 protected:
  
  /// array for temporary u,v,T
  RealVector m_uvT;
  
  /// checks if an function is used in the inlet
  bool m_useFunction;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;
    
  /// physical model var set
  Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet<Physics::Maxwell::Maxwell2DProjectionVarSet> > m_varSet;
  
//   ///x velocity of the Species
//   std::vector<CFreal> _ui;  
//   
//   ///y velocity of the Species
//   std::vector<CFreal> _vi;
//   
//   /// temperatures at Inlet
//   std::vector<CFreal> _Ti;
  
  /// storage for the temporary boundary point coordinates
  RealVector m_bCoord; 
  
  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of string to hold the functions
  std::vector<std::string> m_vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction m_vFunction;  
  
    
}; // end of class BCSubInletUVTEIWRhoiViTi

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_BCSubInletUVTEIWRhoiViTi_hh
