#ifndef COOLFluiD_Numerics_MeshRigidMove_Wedge_hh
#define COOLFluiD_Numerics_MeshRigidMove_Wedge_hh

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
class Wedge : public RigidMoveCom {
public:

  /**
   * Constructor.
   */
  explicit Wedge(const std::string& name) :
    RigidMoveCom(name),
    socket_nodes("nodes")
  {
  }

  /**
   * Destructor.
   */
  ~Wedge()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // data

  // the sink socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

}; // class Wedge

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshRigidMove_Wedge_hh

