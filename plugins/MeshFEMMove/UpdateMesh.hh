#ifndef COOLFluiD_Numerics_MeshFEMMove_UpdateMesh_hh
#define COOLFluiD_Numerics_MeshFEMMove_UpdateMesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "FEMMoveData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshFEMMove {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to setup the MeshData.
   *
   * @author Thomas Wuilbaut
   *
   */
class UpdateMesh : public FEMMoveCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit UpdateMesh(const std::string& name);

  /**
   * Destructor.
   */
  ~UpdateMesh();

  /**
   * Configure
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

private:

  /**
   * Move the mesh of the other namespace using the FEM solution
   */
  void moveMeshNodes();

  /**
   * Write on screen for each show-rate
   */
  void writeOnScreen();

private: // data

  // the sink socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  // the sink socket to the data handle of the node's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for the nodes used in the subsystem
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_otherNodes;

  /// Chronometer to time the simulation
  Common::Stopwatch<Common::CPUTime> _stopwatch;

  /// Precision for outputting to screen
  CFint _precision;

}; // class UpdateMesh

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshFEMMove_UpdateMesh_hh

