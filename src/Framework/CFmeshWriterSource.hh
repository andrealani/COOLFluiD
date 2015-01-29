// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CFmeshWriterSource_hh
#define COOLFluiD_Framework_CFmeshWriterSource_hh

//////////////////////////////////////////////////////////////////////////////

#include "BaseCFMeshFileSource.hh"
#include "Storage.hh"
#include "ElementTypeData.hh"
#include "TopologicalRegionSet.hh"
#include "Framework/State.hh"
#include "DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This base class provides an interface to access the data needed for the
/// mesh creation and manipulation, common to all the I/O operation
/// @author Andrea Lani
class Framework_API CFmeshWriterSource : public BaseCFMeshFileSource {
public: // functions

  /// Constructor
  CFmeshWriterSource();

  /// Destructor
  ~CFmeshWriterSource();

  /// Releases all temporary memory created
  void releaseMemory();

  /// Set the mesh related data
  void setMeshData();

  /// Resize the nodes
  void resizeNodes(const CFuint nbNodes)
  {
    _nodes.resize(nbNodes);
    IndexList<Node>::getList().reset();
    cf_assert(_nodes.size() == nbNodes);
  }

  /// Resize the states
  void resizeStates(const CFuint nbStates)
  {
    _states.resize(nbStates);
    IndexList<State>::getList().reset();
    cf_assert(_states.size() == nbStates);
  }

  /// Get the pastNodes corresponding to the given ID
  const RealVector* getPastNode(const CFuint nodeID) const;

  /// Get the pastStates corresponding to the given ID
  const RealVector* getPastState(const CFuint stateID) const;

  /// Get a vector of extra nodal vars corresponding to the given ID
  RealVector& getExtraNodalValues(const CFuint nodeID);
  
  /// Get the interNodes corresponding to the given ID
  const RealVector* getInterNode(const CFuint nodeID) const;

  /// Get the interStates corresponding to the given ID
  const RealVector* getInterState(const CFuint stateID) const;

  /// Get a vector of extra state vars corresponding to the given ID
  RealVector& getExtraStateValues(const CFuint stateID);
  
   /// Get the vector of extra vars
  RealVector& getExtraValues();

  /// Get the node corresponding to the given ID
  /// @post the template parameter allows to return
  ///         const RealVector& or const Node*&
  CFreal* getNode(const CFuint nodeID)
  {
    return &(*_nodes[nodeID])[0];
  }
  
  /// Get the state corresponding to the given ID
  /// @post the template parameter allows to return
  ///         const RealVector& or const State*&
  CFreal* getState(const CFuint stateID)
  {
    return &(*_states[stateID])[0];
  }
  
  /// Get the state corresponding to the given ID
  /// @post the template parameter allows to return
  ///         const RealVector& or const State*&
  const State* getCoordState(const CFuint stateID) const
  {
    return _states[stateID];
  }

  /// Get the number of nodes in the given element
  CFuint getNbNodesInElement(const CFuint iElem) const
  {
    return _cellNodes->nbCols(iElem);
  }

  /// Get the number of states in the given element
  CFuint getNbStatesInElement(const CFuint iElem) const
  {
    return _cellStates->nbCols(iElem);
  }

  /// Get the element-node
  CFuint getElementNode(const CFuint iElem,
  const CFuint iNode) const
  {
    return (*_cellNodes)(iElem,iNode);
  }

  /// Get the element-state
  CFuint getElementState(const CFuint iElem,
  		 const CFuint iState) const
  {
    return (*_cellStates)(iElem,iState);
  }

  /// Prepare the storage for extra nodal variables
  void prepareNodalExtraVars();

  /// Prepare the storage for extra variables
  void prepareStateExtraVars();
  
  /// Prepare the storage for extra state variables
  void prepareExtraVars();

  /// Sets the extra sockets needed for outputting the extra data
  void setExtraDataSockets(Common::SafePtr<DynamicDataSocketSet<> > dynamicSocket);

private: //helper functions

  /// Set additional mesh data to fully set the
  /// data needed to write CFmesh file
  void copyMeshData();

  /// Set additional TRS data to fully set the
  /// data needed to write CFmesh filet
  void copyTrsData();

  /// Set additional TRS data to fully set the
  /// data needed to write CFmesh filet
  void copyExtraData();

private: //data

  /// handle for the node list
  DataHandle<Node*,GLOBAL>  _nodes;

  /// handle for the state list
  DataHandle<State*,GLOBAL> _states;

  /// cell-node connectivity
  Common::SafePtr<Common::ConnectivityTable<CFuint> > _cellNodes;

  /// cell-state connectivity
  Common::SafePtr<Common::ConnectivityTable<CFuint> > _cellStates;

  /// array of data handles for nodal extra data
  std::vector<Common::SafePtr<DataSocketSink<CFreal> > > _nodalExtraDataSinkSockets;

  /// array of data handles for state extra data
  std::vector<Common::SafePtr<DataSocketSink<CFreal> > > _stateExtraDataSinkSockets;
  
  /// array of data handles for extra data
  std::vector<Common::SafePtr<DataSocketSink<CFreal> > > _extraDataSinkSockets;


}; // end of class CFmeshWriterSource

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshWriterSource_hh
