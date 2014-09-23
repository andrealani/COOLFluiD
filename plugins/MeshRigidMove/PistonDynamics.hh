#ifndef COOLFluiD_Numerics_MeshRigidMove_PistonDynamics_hh
#define COOLFluiD_Numerics_MeshRigidMove_PistonDynamics_hh

//////////////////////////////////////////////////////////////////////////////

#include "RigidMoveData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshRigidMove {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to setup the MeshData.
   *
   * @author Khalil Bensassi 
   * @author Andrea Lani 
   *
   */
class PistonDynamics : public RigidMoveCom {
public:

  /**
   * Constructor.
   */
   PistonDynamics(const std::string& name); 

  /**
   * Destructor.
   */
  ~PistonDynamics();

  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Set up private data
   */
  void setup();

  /**
   * Execute Processing actions
   */
  void execute();

private:
  
  /// compute the volumes
  void computeVolumes();
  
private: // data
  
  CFreal  getdeltaP();
  
  CFreal getAcceleration();
  
  void WritePistonData();
 

  void WriteSpaceTimeData(const CFreal& Xpostion,
			  const CFreal& LogP,
			  CFuint&  Imax, 
			  CFuint& nbface, 
			  std::string filename);  
  
  // the sink socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// storage of ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;
  
  /// storage of volumes
  Framework::DataSocketSink<CFreal> socket_volumes;
  
  /// builder for faces
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceTrsGeoBuilder;

  // builder of GeometricEntity's with Node's
  Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> _geoWithNodesBuilder;

  //Acceleration of the piston at past state
  CFreal _Acceleration;
  CFreal _Pfront; 
  CFreal _Pback;
  CFreal D_piston;
  CFreal L_piston;
  CFreal M_piston;
  CFreal Min_piston;
  CFreal Max_piston;

 //CFreal _oldAcceleration;

  //Speed of the piston at past state
  CFreal _Vpiston;

  //Location of the piston
  CFreal _totalDisplacement;

   CFreal _LogP;

}; // class PistonDynamics

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshRigidMove_PistonDynamics_hh

