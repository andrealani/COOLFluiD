#ifndef COOLFluiD_Numerics_MeshRigidMove_UpdateMesh_hh
#define COOLFluiD_Numerics_MeshRigidMove_UpdateMesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "RigidMoveData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshRigidMove {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to setup the MeshData.
   *
   * @author Thomas Wuilbaut
   *
   */
class UpdateMesh : public RigidMoveCom {
public:

  /**
   * Constructor.
   */
  explicit UpdateMesh(const std::string& name) :
    RigidMoveCom(name),
    socket_nodes("nodes")
  {
  }

  /**
   * Destructor.
   */
  ~UpdateMesh()
  {
  }

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

private: // data

  // the sink socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;


}; // class UpdateMesh

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshRigidMove_UpdateMesh_hh

