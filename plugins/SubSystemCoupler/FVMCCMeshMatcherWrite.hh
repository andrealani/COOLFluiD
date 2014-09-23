#ifndef COOLFluiD_Numerics_SubSystemCoupler_FVMCCMeshMatcherWrite_hh
#define COOLFluiD_Numerics_SubSystemCoupler_FVMCCMeshMatcherWrite_hh

//////////////////////////////////////////////////////////////////////////////

#include "SubSysCouplerData.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to set the
   * match between meshes
   *
   * @author Thomas Wuilbaut
   *
   */

class FVMCCMeshMatcherWrite : public CouplerCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit FVMCCMeshMatcherWrite(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCCMeshMatcherWrite();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Executes the command.
   */
  virtual void execute();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // data

  /**
   * Concrete execution of the command.
   */
  virtual void executeWrite(const CFuint iProc);

  /**
   * Projects the point on the closest face
   * Pairs the projected point with the closest face
   */
  virtual void nodeToElementPairing(const RealVector& coord, CFint& nodeID, RealVector& coordProj);

  /**
   * Writing to a file the acceptance status of the points
   */
  virtual void writeIsAcceptedFile(const std::string dataHandleName);

protected: // data

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;

  /// Name of the interface
  std::string _interfaceName;

  /// Distance between point to project and the closest face
  /// (when the projection is not on any face)
  CFreal _minimumDistanceOffFace  ;

  /// Distance between point to project and the closest face
  /// (when the projection is not on a face)
  CFreal _minimumDistanceOnFace  ;

  /// Temporary variable to count the number of accepted states
  CFuint _nbAcceptedStates;

  SubSysCouplerData::GeoEntityIdx _matchingFace;

  RealVector _shapeFunctionAtCoord;

  bool _shiftStates;

}; // class FVMCCMeshMatcherWrite

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_FVMCCMeshMatcherWrite_hh

