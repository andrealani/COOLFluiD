#ifndef COOLFluiD_Numerics_FiniteVolume_SuperInletPhi3D_hh
#define COOLFluiD_Numerics_FiniteVolume_SuperInletPhi3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Common/BadValueException.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace MultiFluidMHD {      
      class DiffMFMHD3DVarSet;
    }
    
    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
    namespace Maxwell {
      class Maxwell3DProjectionVarSet;
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

class SuperInletPhi3D : public FVMCC_BC { 

public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  SuperInletPhi3D(const std::string& name);
  
  /**
   * Default destructor
   */
  ~SuperInletPhi3D();
  
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
  
  /// array for temporary variables
  RealVector _variables;
  
  /// checks if an function is used in the inlet
  bool _useFunction;
    
  /// physical model var set
  Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet<Physics::Maxwell::Maxwell3DProjectionVarSet> > _updateVarSet;
  
  /// storage for the temporary boundary point coordinates
  RealVector _bCoord; 
  
  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;    
    
  /// option to change to Silver Muller in EM field
  bool _silverMuller;

  /// The electrons are subsonic
  bool _isSubsonic;

}; // end of class SuperInletPhi3D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SuperInletPhi3D_hh
