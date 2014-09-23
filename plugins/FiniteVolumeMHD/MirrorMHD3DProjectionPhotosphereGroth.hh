#ifndef COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DProjectionPhotosphereGroth_hh
#define COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DProjectionPhotosphereGroth_hh

//////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the bc in 3D
   * for projection scheme in the vicinity of the solar surface 
   * (i.e. photosphere) according to Groth et al., JGR Space Physics, 
   * Vol. 105, No. A11, pp. 25053-25078, 2000
   *
   * @author Mehmet Sarp Yalim
   *
   */
class MirrorMHD3DProjectionPhotosphereGroth : public FVMCC_BC {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  MirrorMHD3DProjectionPhotosphereGroth(const std::string& name);

  /**
   * Default destructor
   */
  ~MirrorMHD3DProjectionPhotosphereGroth();

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
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// density value that is to be fixed
  CFreal _rhoFixed;  

  /// transformation matrix from Cartesian to spherical coordinate system for inner state
  RealMatrix _cartesianSphericalTMInnerState;

  /// transformation matrix from spherical to Cartesian coordinate system for inner state
  RealMatrix _sphericalCartesianTMInnerState; 

  /// transformation matrix from Cartesian to spherical coordinate system for ghost state
  RealMatrix _cartesianSphericalTMGhostState;

  /// transformation matrix from spherical to Cartesian coordinate system for ghost state
  RealMatrix _sphericalCartesianTMGhostState;

}; // end of class MirrorMHD3DProjectionPhotosphereGroth

//////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DProjectionPhotosphereGroth_hh
