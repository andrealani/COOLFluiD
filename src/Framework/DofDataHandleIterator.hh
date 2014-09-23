// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DofDataHandleIterator_hh
#define COOLFluiD_Framework_DofDataHandleIterator_hh

//////////////////////////////////////////////////////////////////////////////

#include "ProxyDofIterator.hh"
#include "Storage.hh"
#include "Framework/State.hh"
#include "Common/CFMultiMap.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This function returns a pointer type
/// @author Andrea Lani
template <class PTR, class TYPE>
PTR* getPtr(TYPE* data)
{
  return static_cast<PTR*>(data);
}

template <class TYPE>
TYPE* getPtr(TYPE& data)
{
  return &data;
}

//////////////////////////////////////////////////////////////////////////////

/// This class represents a proxy for a degree of freedom iterator based
/// on DataHandle of an ARRAY type convertable to RETURNTYPE. It
/// offers a uniform interface
/// to use underlying different
/// storages.
/// @author Andrea Lani
template <class RETURNTYPE, class ARRAY, class TAG = LOCAL>
class DofDataHandleIterator : public ProxyDofIterator<RETURNTYPE> {
public:
  
  /// Default constructor
  DofDataHandleIterator(DataHandle<ARRAY*, TAG>& handle) :
    ProxyDofIterator<RETURNTYPE>(),
    _handle(handle)
  {
  }
  
  /// Default destructor
  ~DofDataHandleIterator()
  {
  }
  
  /// Gets the dof corresponding to the given ID
  RETURNTYPE* getState(const CFuint nodeID)
  {
    return getPtr<RETURNTYPE, ARRAY>(_handle[nodeID]);
  }
  
  /// Gets the dof corresponding to the given ID
  RETURNTYPE* getNode(const CFuint nodeID)
  {
    return CFNULL;
  }
  
  /// Gets the ID of the node in this dof
  CFuint getNodeLocalID(const CFuint nodeID)
  {
    return nodeID;
  }
  
  /// Gets the ID of the state in this dof
  CFuint getStateLocalID(const CFuint nodeID)
  {
    return nodeID;
  }

  /// Gets the size
  CFuint getSize() const
  {
    return _handle.size();
  }

private:
  
  /// handle to the states
  DataHandle< ARRAY*, TAG> _handle;
  
}; // end of class DofDataHandleIterator

//////////////////////////////////////////////////////////////////////////////

/// This class represents a proxy for a degree of freedom iterator based
/// on DataHandle of an State pointer.
/// storages.
/// @author Andrea Lani
template <class RETURNTYPE, class TAG>
class DofDataHandleIterator<RETURNTYPE, State, TAG> : 
  public ProxyDofIterator<RETURNTYPE> {
public:

  /// Default constructor
  DofDataHandleIterator(DataHandle<State*, TAG>& handle,
  std::vector<CFuint> *const nodeIdToStateId) :
    ProxyDofIterator<RETURNTYPE>(),
    _handle(handle),
    _nodeIdToStateId(nodeIdToStateId)
  {
  }
  
  /// Default destructor
  ~DofDataHandleIterator()
  {
  }
  
  /// Gets the dof corresponding to the given ID
  RETURNTYPE* getState(const CFuint nodeID)
  {
    return getPtr<RETURNTYPE, State>
      (_handle[(*_nodeIdToStateId)[nodeID]]);
  }
  
  /// Gets the dof corresponding to the given ID
  RETURNTYPE* getNode(const CFuint nodeID)
  {
    return getPtr<RETURNTYPE, Node>
      (&_handle[(*_nodeIdToStateId)[nodeID]]->getCoordinates());
  }
  
  /// Gets the ID of the node in this dof
  CFuint getNodeLocalID(const CFuint nodeID)
  {
    return  _handle[(*_nodeIdToStateId)[nodeID]]->getCoordinates().getLocalID();
  }

  /// Gets the ID of the state in this dof
  CFuint getStateLocalID(const CFuint nodeID)
  {
    return _handle[(*_nodeIdToStateId)[nodeID]]->getLocalID();
  }
  
  /// Gets the size
  CFuint getSize() const
  {
    return _handle.size();
  }
  
private:
  
  /// handle to the states
  DataHandle<State*, TAG> _handle;
  
  /// mapping between local nodeID and local stateID
  std::vector<CFuint> *const _nodeIdToStateId;
  
}; // end of class DofDataHandleIterator

//////////////////////////////////////////////////////////////////////////////

/// This class represents a proxy for a degree of freedom iterator based
/// on DataHandle of an State pointer.
/// storages.
/// @author Andrea Lani
template <class TAG>
class DofDataHandleIterator<CFreal, State, TAG> : 
  public ProxyDofIterator<CFreal> {
public:

  /// Default constructor
  DofDataHandleIterator(DataHandle<State*, TAG>& handle,
			std::vector<CFuint> *const nodeIdToStateId) :
    ProxyDofIterator<CFreal>(),
    _handle(handle),
    _nodeIdToStateId(nodeIdToStateId)
  {
  }
  
  /// Default destructor
  ~DofDataHandleIterator()
  {
  }
  
  /// Gets the dof corresponding to the given ID
  CFreal* getState(const CFuint nodeID)
  {
    return (&(*_handle[(*_nodeIdToStateId)[nodeID]])[0]);
  }
  
  /// Gets the dof corresponding to the given ID
  CFreal* getNode(const CFuint nodeID)
  {
    return (&(_handle[(*_nodeIdToStateId)[nodeID]]->getCoordinates())[0]);
  }
  
  /// Gets the ID of the node in this dof
  CFuint getNodeLocalID(const CFuint nodeID)
  {
    return _handle[(*_nodeIdToStateId)[nodeID]]->getCoordinates().getLocalID();
  }
  
  /// Gets the ID of the state in this dof
  CFuint getStateLocalID(const CFuint nodeID)
  {
    return _handle[(*_nodeIdToStateId)[nodeID]]->getLocalID();
  }
  
  /// Gets the size
  CFuint getSize() const
  {
    return _handle.size();
  }
  
private:
  
  /// handle to the states
  DataHandle<State*, TAG> _handle;
  
  /// mapping between local nodeID and local stateID
  std::vector<CFuint> *const _nodeIdToStateId;
  
}; // end of class DofDataHandleIterator

//////////////////////////////////////////////////////////////////////////////

/// This class represents a proxy for a degree of freedom iterator based
/// on DataHandle of an State pointer.
/// storages.
/// @author Andrea Lani
template <class TAG>
class DofDataHandleIterator<CFreal, CFreal, TAG> : 
  public ProxyDofIterator<CFreal> {
public:

  /// Default constructor
  DofDataHandleIterator(std::vector<CFreal>* statesNodes,
			Common::CFMultiMap<CFuint,CFuint>* mapIDs,
			CFuint stateStride, 
			CFuint nodeStride) :
    ProxyDofIterator<CFreal>(),
    _statesNodes(statesNodes),
    _mapIDs(mapIDs),
    _ss(stateStride),
    _stride(_ss + nodeStride)
  {
  }
  
  /// Default destructor
  ~DofDataHandleIterator()
  {
  }
  
  /// Gets the state vector corresponding to the given ID
  CFreal* getState(const CFuint nodeID)
  {
    return (&(*_statesNodes)[getStateLocalID(nodeID)*_stride]);
  }
  
  /// Gets the node coordinates corresponding to the given ID
  CFreal* getNode(const CFuint nodeID)
  {
    return (&(*_statesNodes)[getNodeLocalID(nodeID)*_stride + _ss]);
  }
  
  /// Gets the ID of the node in this dof
  CFuint getNodeLocalID(const CFuint nodeID)
  {
    bool found = false;
    return _mapIDs->find(nodeID, found).first->second;
    cf_assert(found);
  }
  
  /// Gets the ID of the state in this dof
  CFuint getStateLocalID(const CFuint nodeID)
  { 
    bool found = false;
    return _mapIDs->find(nodeID, found).first->second;
    cf_assert(found);
  }
  
  /// Gets the size
  CFuint getSize() const
  {
    return _statesNodes->size()/_stride;
  }
  
private:
  
  /// array of states
  std::vector<CFreal>* _statesNodes;

  // map the input ID with a corresponding (local) state ID
  Common::CFMultiMap<CFuint, CFuint>* _mapIDs;
  
  /// states stride
  const CFuint _ss;
  
  /// states+nodes stride
  const CFuint _stride;
  
}; // end of class DofDataHandleIterator

//////////////////////////////////////////////////////////////////////////////

/// This class represents a proxy for a degree of freedom iterator based
/// on DataHandle of an ARRAY type convertable to RealVector. It
/// offers a uniform interface
/// to use underlying different
/// storages.
/// @author Andrea Lani
template <class ARRAY, class TAG>
class DofDataHandleIterator<ARRAY, ARRAY, TAG> :
  public ProxyDofIterator<ARRAY> {
public:
  
  /// Default constructor
  DofDataHandleIterator(DataHandle<ARRAY, TAG>& handle) :
    ProxyDofIterator<ARRAY>(),
    _handle(handle)
  {
  }

  /// Default destructor
  ~DofDataHandleIterator()
  {
  }

  /// Gets the dof corresponding to the given ID
  ARRAY* getState(const CFuint nodeID)
  {
    return getPtr<ARRAY>(_handle[nodeID]);
  }
  
  /// Gets the dof corresponding to the given ID
  ARRAY* getNode(const CFuint nodeID)
  {
    return CFNULL;
  }

  /// Gets the ID of the node in this dof
  CFuint getNodeLocalID(const CFuint nodeID)
  {
    return nodeID;
  }

  /// Gets the ID of the state in this dof
  CFuint getStateLocalID(const CFuint nodeID)
  {
    return nodeID;
  }

  /// Gets the size
  CFuint getSize() const
  {
    return _handle.size();
  }

private:
  
  /// handle to the states
  DataHandle< ARRAY, TAG> _handle;
  
}; // end of class DofDataHandleIterator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DofDataHandleIterator_hh
