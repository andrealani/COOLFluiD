#ifndef COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DTanakaPFix_hh
#define COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DTanakaPFix_hh

//////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the mirror bc in 3D
   * as in Powell et. al. JCP Vol.154 1999 pp.284-309
   * for Powell's source term scheme according to Tanaka's approach
   *
   * @author Mehmet Sarp Yalim
   * @author Andrea Lani
   *
   */
class MirrorMHD3DTanakaPFix : public FVMCC_BC {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  MirrorMHD3DTanakaPFix(const std::string& name);

  /**
   * Default destructor
   */
  ~MirrorMHD3DTanakaPFix();

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
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
  
 private:

  /// physical model var set
  Common::SafePtr<Physics::MHD::MHD3DVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// density value that is to be fixed
  CFreal _rhoFixed;  

  /// pressure value that is to be fixed
  CFreal _pFixed;

}; // end of class MirrorMHD3DTanakaPFix

//////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DTanakaPFix_hh
