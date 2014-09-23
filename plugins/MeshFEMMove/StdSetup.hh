#ifndef COOLFluiD_Numerics_MeshFEMMove_StdSetup_hh
#define COOLFluiD_Numerics_MeshFEMMove_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "MeshFEMMove/FEMMoveData.hh"

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
class StdSetup : public FEMMoveCom {
public:

  /**
   * Constructor.
   */
  explicit StdSetup(const std::string& name);

  /**
   * Destructor.
   */
  ~StdSetup()
  {
  }

  /**
   * Configure
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();


private: //helper functions

  /**
   * Create the FEM connectivity if it does not exist yet
   */
  void createExtraConnectivity();

protected: //member data

  /// socket for the states used for the FEM
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_femStates;

  /// socket for nodes used for the FEM
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_femNodes;

  /// socket for the states used in the subsystem
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_otherStates;

  /// socket for the nodes used in the subsystem
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_otherNodes;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshFEMMove_StdSetup_hh

