#ifndef COOLFluiD_Numerics_FiniteVolume_SubOutletEuler3D_hh
#define COOLFluiD_Numerics_FiniteVolume_SubOutletEuler3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class EulerVarSet;
    }
  }
  
  namespace Numerics {

    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a subsonic outlet command 
 * 
 * @author Mehmet Sarp Yalim
 * @author Andrea Lani
 *
 */
class SubOutletEuler3D : public FVMCC_BC { 
public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  SubOutletEuler3D(const std::string& name);
  
  /**
   * Default destructor
   */
  ~SubOutletEuler3D();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  void setup();
  
  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);
      
private: // data
  
  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> _varSet;
  
  /// physical model data
  RealVector _dataInnerState;
  
  /// physical model data
  RealVector _dataGhostState;
  
  /// static pressure
  CFreal _pressure;
  
}; // end of class SubOutletEuler3D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubOutletEuler3D_hh 
