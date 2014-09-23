#ifndef COOLFluiD_DiscontGalerkin_StdSetup_hh
#define COOLFluiD_DiscontGalerkin_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "DiscontGalerkin/DiscontGalerkinSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a standard command to setup the DiscontGalerkin method
 * @author Martin Holik
 * @author Vit Dolejsi
 */
class StdSetup : public DiscontGalerkinSolverCom {

public: // functions

  /// Constructor
  explicit StdSetup(const std::string& name);

  /// Destructor
  ~StdSetup();

  /// Execute processing actions
  void execute();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

  /// Node to state ID map
  std::vector< CFuint > m_nodeIDToStateID;

protected: // data

  /**
   * socket with proxy to be able to use the data handle of nodal states
   * uniformly independently from the actual storage type being RealVector or
   * State*
   */
  Framework::DataSocketSource< Framework::ProxyDofIterator< RealVector >* >
    socket_nstatesProxy;

  /// socket for state's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// socket for node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;

  /// socket for integration indexes
  Framework::DataSocketSource< std::vector<CFuint> >
    socket_integrationIndex;

  /// socket for normals on face
  Framework::DataSocketSource< RealVector >
    socket_normals;

  /// socket for visualization - state to node connection
  Framework::DataSocketSource< CFuint >
    socket_techStatesToNodes;

  /// socket for visualization - connection techNodes to statesID
  Framework::DataSocketSource< std::vector < CFuint > >
    socket_techNodesToStates;

  /// socket for visualization - coordinates of nodes for visualization
  Framework::DataSocketSource< RealVector >
    socket_techNodesCoordinates;

};  // class StdSetup

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_StdSetup_hh

