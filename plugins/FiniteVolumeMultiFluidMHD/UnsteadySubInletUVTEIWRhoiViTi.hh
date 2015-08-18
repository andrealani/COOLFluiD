#ifndef COOLFluiD_Numerics_FiniteVolume_UnsteadySubInletUVTEIWRhoiViTi_hh
#define COOLFluiD_Numerics_FiniteVolume_UnsteadySubInletUVTEIWRhoiViTi_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Common/BadValueException.hh"


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
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Subsonic Inlet imposing the Velocity and Temperature
   * Maxwell Equations: Perfectly Conducting Wall Condition
   * 
   * @author Alejandro Alvarez
   *
   */

class UnsteadySubInletUVTEIWRhoiViTi : public FVMCC_BC {

public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  UnsteadySubInletUVTEIWRhoiViTi(const std::string& name);
  
  /**
   * Default destructor
   */
  ~UnsteadySubInletUVTEIWRhoiViTi();
  
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
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

 protected:
  
  /// array for temporary u,v,T
  RealVector _uvT;
  
  /// checks if an function is used in the inlet
  bool _useFunction;
    
  /// physical model var set
  Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet<Physics::Maxwell::Maxwell2DProjectionVarSet> > _updateVarSet;
  
//   ///x velocity of the Species
//   std::vector<CFreal> _ui;  
//   
//   ///y velocity of the Species
//   std::vector<CFreal> _vi;
//   
//   /// temperatures at Inlet
//   std::vector<CFreal> _Ti;

  /// Vector for coordinates + time
  RealVector _variables;
  
  /// storage for the temporary boundary point coordinates
  RealVector _bCoord; 
  
  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;    
    
}; // end of class UnsteadySubInletUVTEIWRhoiViTi

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_UnsteadySubInletUVTEIWRhoiViTi_hh
