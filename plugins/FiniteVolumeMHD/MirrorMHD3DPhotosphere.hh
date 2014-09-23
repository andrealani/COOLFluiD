#ifndef COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DPhotosphere_hh
#define COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DPhotosphere_hh

//////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/DataSocketSink.hh"

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
   * This class represents a command that applies the bc in 3D
   * for Powell's source term in the vicinity of the solar surface (i.e. photosphere)
   *
   * @author Mehmet Sarp Yalim
   *
   */
class MirrorMHD3DPhotosphere : public FVMCC_BC {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  MirrorMHD3DPhotosphere(const std::string& name);

  /**
   * Default destructor
   */
  ~MirrorMHD3DPhotosphere();

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
   * Compute the coronal magnetic field in ghost cells using PFSS technique 
   */
  void computePFSSMagneticField(const RealVector& stateCoordsSpherical,
                                RealVector& BPFSSSpherical);

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
  Common::SafePtr<Physics::MHD::MHD3DVarSet> _varSet;

  /// socket for the real part of the PFSS spherical harmonics coefficients, Alm, to be used in the boundary conditions
  Framework::DataSocketSink<std::vector<CFreal> > socket_Almreal;

  /// socket for the imaginary part of the PFSS spherical harmonics coefficients, Alm, to be used in the boundary conditions
  Framework::DataSocketSink<std::vector<CFreal> > socket_Almimg;

  /// socket for the real part of the PFSS spherical harmonics coefficients, Blm, to be used in the boundary conditions
  Framework::DataSocketSink<std::vector<CFreal> > socket_Blmreal;

  /// socket for the imaginary part of the PFSS spherical harmonics coefficients, Blm, to be used in the boundary conditions
  Framework::DataSocketSink<std::vector<CFreal> > socket_Blmimg;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// density value that is to be fixed
  CFreal _rhoFixed;  

  /// temperature value that is to be fixed
  CFreal _TFixed;

  /// transformation matrix from Cartesian to spherical coordinate system for inner state
  RealMatrix _cartesianSphericalTMInnerState;

  /// transformation matrix from spherical to Cartesian coordinate system for inner state
  RealMatrix _sphericalCartesianTMInnerState; 

  /// transformation matrix from Cartesian to spherical coordinate system for ghost state
  RealMatrix _cartesianSphericalTMGhostState;

  /// transformation matrix from spherical to Cartesian coordinate system for ghost state
  RealMatrix _sphericalCartesianTMGhostState;

}; // end of class MirrorMHD3DPhotosphere

//////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DPhotosphere_hh
