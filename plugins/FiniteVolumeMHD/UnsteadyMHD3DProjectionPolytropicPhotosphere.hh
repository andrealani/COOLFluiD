#ifndef COOLFluiD_Numerics_FiniteVolume_UnsteadyMHD3DProjectionPolytropicPhotosphere_hh
#define COOLFluiD_Numerics_FiniteVolume_UnsteadyMHD3DProjectionPolytropicPhotosphere_hh

//////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionPolytropicVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies a time-dependent bc in 3D
   * for projection scheme in the vicinity of the solar surface (i.e. photosphere)
   * on the polytropic modelling of the solar wind
   *
   * @author Mehmet Sarp Yalim
   *
   */
class UnsteadyMHD3DProjectionPolytropicPhotosphere : public FVMCC_BC {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  UnsteadyMHD3DProjectionPolytropicPhotosphere(const std::string& name);

  /**
   * Default destructor
   */
  ~UnsteadyMHD3DProjectionPolytropicPhotosphere();

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
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

 private:

  /// physical model var set
  Common::SafePtr<Physics::MHD::MHD3DProjectionPolytropicVarSet> _varSet;

  /// socket for the PFSS magnetic field components in Cartesian coordinates to be computed once in the setup phase
  Framework::DataSocketSink<std::vector<CFreal> > socket_BPFSS;

  /// socket for the PFSS magnetic field components in Cartesian coordinates to be computed once in the setup phase
  Framework::DataSocketSource<std::vector<CFreal> > socket_BPFSSGhostBegin;

  /// socket for the PFSS magnetic field components in Cartesian coordinates to be computed once in the setup phase
  Framework::DataSocketSource<std::vector<CFreal> > socket_BPFSSGhostEnd;

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

  /// begin time of the overall unsteady simulation in case of restart (in seconds)
  CFreal _restartTime;

  /// end time of the overall unsteady simulation (in seconds)
  CFreal _endTime;

  /// initial time of the overall unsteady simulation (in seconds)
  CFreal _initTime;

}; // end of class UnsteadyMHD3DProjectionPolytropicPhotosphere

//////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_UnsteadyMHD3DProjectionPolytropicPhotosphere_hh
