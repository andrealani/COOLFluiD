#ifndef COOLFluiD_Numerics_FiniteVolume_MagnetosphericWall_hh
#define COOLFluiD_Numerics_FiniteVolume_MagnetosphericWall_hh

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

class MagnetosphericWall : public FVMCC_BC { 

public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  MagnetosphericWall(const std::string& name);
  
  /**
   * Default destructor
   */
  ~MagnetosphericWall();
  
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
  RealVector _rhoT;
  
  /// checks if an function is used in the inlet
  bool _useFunction;
    
  /// physical model var set
  Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet<Physics::Maxwell::Maxwell3DProjectionVarSet> > _updateVarSet;
  
//   ///x velocity of the Species
//   std::vector<CFreal> _ui;  
//   
//   ///y velocity of the Species
//   std::vector<CFreal> _vi;
//   
//   /// temperatures at Inlet
//   std::vector<CFreal> _Ti;
  
  /// storage for the temporary boundary point coordinates
  RealVector _bCoord; 
  
  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;    
    
}; // end of class MagnetosphericWall

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MagnetosphericWall_hh
