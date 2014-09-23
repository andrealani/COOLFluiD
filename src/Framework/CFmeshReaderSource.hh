// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CFmeshReaderSource_hh
#define COOLFluiD_Framework_CFmeshReaderSource_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ConnectivityTable.hh"
#include "Common/CFMap.hh"

#include "Framework/BaseCFMeshFileSource.hh"
#include "Framework/Storage.hh"
#include "Framework/State.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This base class provides an interface to access the data needed for the
/// mesh creation and manipulation, common to all the I/O operation
/// @author Andrea Lani
/// @author Tiago Quintino
/// @author Dries Kimpe
class Framework_API CFmeshReaderSource : public BaseCFMeshFileSource {
public: // functions

  /// Constructor
  CFmeshReaderSource();

  /// Destructor
  ~CFmeshReaderSource();

  /// Releases all temporary memory created
  void releaseMemory();

  /// Set the mesh related data
  void setDataSockets(DataSocketSink<State*,GLOBAL> statesSocket,
                      DataSocketSink<Node*,GLOBAL> nodesSocket,
                      Common::SafePtr<DynamicDataSocketSet<> > dynamicSocket);

  /// Set the mesh related data
  void setDataSockets(DataSocketSink<State*,GLOBAL> statesSocket,
                      DataSocketSink<Node*,GLOBAL> nodesSocket);

  /// Resize the nodes
  void resizeNodes(const CFuint nbNodes);

  /// Resize the states
  void resizeStates(const CFuint nbStates);

  /// resizeExtraVars
  void resizeExtraVars();
  
  /// Set the node corresponding to the given ID
  void setNode(const CFuint nodeID, const RealVector& value);

  /// Set the past node corresponding to the given ID
  void setPastNode(const CFuint nodeID, const RealVector& value);

  /// Set the intermediate node corresponding to the given ID
  void setInterNode(const CFuint nodeID, const RealVector& value);

  /// Set the extra variables corresponding to the given ID
  void setNodalExtraVar(const CFuint nodeID, const RealVector& value);

  /// Set the state corresponding to the given ID
  void setState(const CFuint stateID, const RealVector& value);

  /// Set the past state corresponding to the given ID
  void setPastState(const CFuint stateID, const RealVector& value);
  /// Set the intermediate state corresponding to the given ID
  void setInterState(const CFuint stateID, const RealVector& value);

  /// Set the state corresponding to the given ID
  void setStateExtraVar(const CFuint stateID, const RealVector& value);

  /// Set the extra variables corresponding to the given ID
  void setExtraVar(const RealVector& value);

  /// Prepare the storage for extra nodal variables
  void prepareNodalExtraVars();

  /// Prepare the storage for extra state variables
  void prepareStateExtraVars();
  
  /// Prepare the storage for extra state variables
  void prepareExtraVars();

  /// Set the node corresponding to the given ID
  /// (preallocated memory) -- DO NOT USE
  void setNode(const CFuint nodeID, CFreal * Mem, const RealVector & D, bool IsUpdatable);

  /// Set the state corresponding to the given ID
  /// (preallocated memory) -- DO NOT USE
  void setState(const CFuint stateID, CFreal * Mem, const RealVector & D, bool IsUpdatable);

  /// Set the node corresponding to the given ID
  Node* createNode(const CFuint nodeID, const RealVector& value);

  /// Set the state corresponding to the given ID
  State* createState(const CFuint stateID, const RealVector& value);

  /// Set the node corresponding to the given ID
  /// (preallocated memory) -- DO NOT USE
  Node* createNode(const CFuint nodeID, CFreal * Mem, const RealVector & D, bool IsUpdatable);

  /// Set the state corresponding to the given ID
  /// (preallocated memory) -- DO NOT USE
  State* createState(const CFuint stateID, CFreal * Mem, const RealVector & D, bool IsUpdatable);

  /// Get the node corresponding to the given ID
  /// @post the template parameter allows to return
  ///         const RealVector& or const Node*&
  template <class DOFTYPE>
  const DOFTYPE* getNode(const CFuint nodeID) const
  {
    return (socket_nodes.getDataHandle())[nodeID];
  }

  /// Get the state corresponding to the given ID
  /// @post the template parameter allows to return
  ///         const RealVector& or const State*&
  template <class DOFTYPE>
  const DOFTYPE* getState(const CFuint stateID) const
  {
    return socket_states.getDataHandle()[stateID];
  }

  /// for the MeshDataSourceInterface
  template <class OUTPUT>
  void getStateCopy (const CFuint StateID, OUTPUT Out) const
  {
      *Out = *getState<RealVector>(StateID);
  }

  /// from the MeshDataSourceInterface
  template <class OUTPUT>
  void getNodeCopy (const CFuint NodeID, OUTPUT Out) const
  {
      *Out = *getNode<RealVector>(NodeID);
  }

  /// Resize the element-node connectivity
  void resizeElementNode(const std::valarray<CFuint>& nbCols)
  {
    _elementNode->resize(nbCols);
  }

  /// Resize the element-state connectivity
  void resizeElementState(const std::valarray<CFuint>& nbCols)
  {
    _elementState->resize(nbCols);
  }

  /// Get the number of nodes in the given element
  CFuint getNbNodesInElement(const CFuint iElem) const
  {
    return _elementNode->nbCols(iElem);
  }

  /// Get the number of states in the given element
  CFuint getNbStatesInElement(const CFuint iElem) const
  {
    return _elementState->nbCols(iElem);
  }

  /// Set the element-node
  void setElementNode(const CFuint iElem,
  	      const CFuint iNode,
  	      const CFuint value)
  {
    (*_elementNode)(iElem, iNode) = value;
  }

  /// Set the element-state
  void setElementState(const CFuint iElem,
  	       const CFuint iState,
  	       const CFuint value)
  {
    (*_elementState)(iElem, iState) = value;
  }

  /// Get the element-node
  CFuint getElementNode(const CFuint iElem,
  const CFuint iNode) const
  {
    return (*_elementNode)(iElem, iNode);
  }

  /// Gets the table of connectivity of element to state
  Common::SafePtr<Common::ConnectivityTable<CFuint> >  getElementStateTable()
  {
    return _elementState;
  }
  
  /// Gets the table of connectivity of element to node
  Common::SafePtr<Common::ConnectivityTable<CFuint> >  getElementNodeTable()
  {
    return _elementNode;
  }

  /// Get the element-state
  CFuint getElementState(const CFuint iElem,
  		 const CFuint iState) const
  {
    return (*_elementState)(iElem, iState);
  }

  /// return true if the state is owned by this cpu
  bool isLocalState (const CFuint StateNum) const
  {
      cf_assert (StateNum < _StateOwnership.size());
      return _StateOwnership[StateNum];
  }

  /// set the state ownership
  void setLocalState (const CFuint StateNum, bool IsLocal=true)
  {
      cf_assert (StateNum < _StateOwnership.size());
      _StateOwnership[StateNum]=IsLocal;
  }

  /// return true if the node is owned by this cpu
  bool isLocalNode (const CFuint NodeNum) const
  {
      cf_assert (NodeNum < _NodeOwnership.size());
      return _NodeOwnership[NodeNum];
  }

  /// set the node ownership
  void setLocalNode (const CFuint NodeNum, bool IsLocal=true)
  {
      cf_assert (NodeNum < _NodeOwnership.size());
      _NodeOwnership[NodeNum]=IsLocal;
  }


  /// Set Local To global node mapping
  void setNodeLocalToGlobal (const CFuint NodeNum, const CFuint Global)
  {
      cf_assert (NodeNum < _NodeLocalToGlobal.size());
      _NodeLocalToGlobal[NodeNum]=Global;
  }

  /// get the Local To Global node mapping
  CFuint getNodeLocalToGlobal (const CFuint NodeNUm) const
  {
      cf_assert (NodeNUm < _NodeLocalToGlobal.size());
      return _NodeLocalToGlobal[NodeNUm];
  }

  /// get the Local To Global state mapping
  CFuint getStateLocalToGlobal (const CFuint StateNum) const
  {
      cf_assert (StateNum < _StateLocalToGlobal.size());
      return _StateLocalToGlobal[StateNum];
  }

  /// set the Local To Global State mapping
  void setStateLocalToGlobal (const CFuint StateNum, const CFuint Global)
  {
      cf_assert (StateNum < _StateLocalToGlobal.size());
      _StateLocalToGlobal[StateNum]=Global;
  }

  /// return the nodes and the states of the element in a PairNodeState
  /// structure pointed to by IteratorOut
  void getElementData (const CFuint EleNum, PairNodeState & Out) const;

  /// Accessor to the DataHandle of the Nodes
  /// MeshDataBuilder's need it
  /// @return DataHandle to the "nodes"
  DataHandle<Node*,GLOBAL> getNodesHandle()   { return socket_nodes.getDataHandle(); }

  /// Accessor to the DataHandle of the States
  /// MeshDataBuilder's need it
  /// @return DataHandle to the "states"
  DataHandle<State*,GLOBAL> getStatesHandle() { return socket_states.getDataHandle(); }

private: // data

  /// handle for the node list
  DataSocketSink<Node*,GLOBAL>  socket_nodes;

  /// handle for the state list
  DataSocketSink<State*,GLOBAL> socket_states;

  /// table for the element-node connectivity
  Common::SafePtr<Common::ConnectivityTable<CFuint> > _elementNode;

  /// table for the element-state connectivity
  Common::SafePtr<Common::ConnectivityTable<CFuint> > _elementState;

  /// array of data handles for nodal extra data
  std::vector<Common::SafePtr<DataSocketSource<CFreal> > > _nodalExtraDataSourceSockets;

  /// array of data handles for state extra data
  std::vector<Common::SafePtr<DataSocketSource<CFreal> > > _stateExtraDataSourceSockets;
  
  /// array of data handles for extra data
  std::vector<Common::SafePtr<DataSocketSource<CFreal> > > _extraDataSourceSockets;

  /// mapping between local and global states
  std::vector<CFuint> _StateLocalToGlobal;

  /// mapping between local and global nodes
  std::vector<CFuint> _NodeLocalToGlobal;

  /// state ownership flags
  std::vector<bool> _StateOwnership;

  /// node ownership flags
  std::vector<bool> _NodeOwnership;

  /// total number of state related extra variables
  CFuint _totalNbStateExtraVars;

  /// total number of node related extra variables
  CFuint _totalNbNodalExtraVars;
  
  /// total number of extra variables
  CFuint _totalNbExtraVars;

}; // end of class CFmeshReaderSource

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshReaderSource_hh
