#ifndef COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DProjectionPhotosphere_hh
#define COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DProjectionPhotosphere_hh

//////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/DataSocketSink.hh"

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
   * for projection scheme on volumetric heating modelling of the solar wind
   * in the vicinity of the solar surface (i.e. photosphere)
   *
   * @author Mehmet Sarp Yalim
   *
   */
class MirrorMHD3DProjectionPhotosphere : public FVMCC_BC {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  MirrorMHD3DProjectionPhotosphere(const std::string& name);

  /**
   * Default destructor
   */
  ~MirrorMHD3DProjectionPhotosphere();

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

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

 private:

  /// physical model var set
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> _varSet;

  /// socket for the PFSS magnetic field components in Cartesian coordinates to be computed once in the setup phase
  Framework::DataSocketSink<std::vector<CFreal> > socket_BPFSS;

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

}; // end of class MirrorMHD3DProjectionPhotosphere

//////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DProjectionPhotosphere_hh
