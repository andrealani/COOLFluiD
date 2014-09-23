#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DUVT_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DUVT_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }
  
  namespace Numerics {

    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a subsonic inlet command with the initial conditions given for tTotal, pTotal and alpha 
   * 
   * @author Mehmet Sarp Yalim
   * @author Andrea Lani
   *
   */
class SubInletEuler2DUVT : public FVMCC_BC { 
public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletEuler2DUVT(const std::string& name);
  
  /**
   * Default destructor
   */
  ~SubInletEuler2DUVT();
  
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
  
  /// physical model var set
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;
    
  /// physical model data
  RealVector _dataInnerState;
  
  /// physical model data
  RealVector _dataGhostState;
  
  /// x velocity
  CFreal                                   _uinf;
  
  /// y velocity
  CFreal                                   _vinf;
  
  /// static temperature
  CFreal                                   _temperature;
  
}; // end of class SubInletEuler2DUVT

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DUVT_hh 
