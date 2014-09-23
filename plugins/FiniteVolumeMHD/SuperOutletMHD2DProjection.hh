#ifndef COOLFluiD_Numerics_FiniteVolume_SuperOutletMHD2DProjection_hh
#define COOLFluiD_Numerics_FiniteVolume_SuperOutletMHD2DProjection_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD2DProjectionVarSet;
    }
  }
  
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a supersonic outlet command in 2D
   * for projection scheme
   * 
   * @author Andrea Lani
   * @author Mehmet Sarp Yalim
   *
   */
class SuperOutletMHD2DProjection : public FVMCC_BC {
public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  SuperOutletMHD2DProjection(const std::string& name);
  
  /**
   * Default destructor
   */
  ~SuperOutletMHD2DProjection();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Configures this object by complementing the 
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

private: //data

  /// phi value that is to be fixed
  CFreal _refPhi;

  /// physical model var set
  Common::SafePtr<Physics::MHD::MHD2DProjectionVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

}; // end of class SuperOutletMHD2DProjection

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SuperOutletMHD2DProjection_hh
