// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CFmeshReaderWriterSource_hh
#define COOLFluiD_Framework_CFmeshReaderWriterSource_hh

//////////////////////////////////////////////////////////////////////////////

#include "BaseCFMeshFileSource.hh"

#include "MathTools/RealVector.hh"
#include "Common/Table.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class Node;
    class State;

//////////////////////////////////////////////////////////////////////////////

/// This base class provides an interface to access the data needed for the
/// mesh creation and manipulation, common to all the I/O operation
/// @author Andrea Lani
class Framework_API CFmeshReaderWriterSource : public BaseCFMeshFileSource {
public: // functions

  /// Constructor
  CFmeshReaderWriterSource();

  /// Destructor
  ~CFmeshReaderWriterSource();

  /// Releases all temporary memory created
  void releaseMemory();

  /// Prepare the storage for extra nodal variables
  void prepareNodalExtraVars()
  {
    // throw Common::NotImplementedException (FromHere(),"CFmeshReaderWriterSource::prepareNodalExtraVars()");
  }

  /// Prepare the storage for extra state variables
  void prepareStateExtraVars()
  {
    // throw Common::NotImplementedException (FromHere(),"CFmeshReaderWriterSource::prepareStateExtraVars()");
  }
  
  /// Prepare the storage for extra state variables
  void prepareExtraVars()
  {
    // throw Common::NotImplementedException (FromHere(),"CFmeshReaderWriterSource::prepareExtraVars()");
  }

  /// Resize the nodes
  void resizeNodes(const CFuint nbNodes)
  {
    // to safely resize FREE the memory first
    //(this is NOT done by std::resize())
    if (_nodes.size() > 0) {
      std::vector<CFreal>().swap(_nodes);
    }
    _nodes.resize(nbNodes*_dimension);
  }

  /// Resize the states
  void resizeStates(const CFuint nbStates)
  {
    // to safely resize FREE the memory first
    //(this is NOT done by std::resize())
    if (_states.size() > 0) {
      std::vector<CFreal>().swap(_states);
    }
    _states.resize(nbStates*_nbEquations);
  }
  
  /// Does not actually create a Node
  /// Set the node corresponding to the given ID
  /// @pre value.size() ==  (*_nodes)[nodeID].size()
  void setNode(const CFuint nodeID, const RealVector& value)
  {
    cf_assert(value.size() == _dimension);
    const CFuint start = nodeID*_dimension;
    for (CFuint i=0; i < _dimension; ++i) {
      _nodes[start+i] = value[i];
    }
  }
  
  /// Does not actually create a State
  /// Set the state corresponding to the given ID
  void setState(const CFuint stateID, const RealVector& value)
  {
    cf_assert(value.size() == _nbEquations);
    const CFuint start = stateID*_nbEquations;
    for (CFuint i=0; i < _nbEquations; ++i) {
      _states[start+i] = value[i];
    }
  }
  
  /// Set the extra variables corresponding to the given ID
  void setNodalExtraVar(const CFuint nodeID, const RealVector& value)
  {
    //nothing to be done here
  }


  /// Set the state corresponding to the given ID
  void setStateExtraVar(const CFuint stateID, const RealVector& value)
  {
    //nothing to be done here
  }

  /// Get the pastNodes corresponding to the given ID
  const RealVector* getPastNode(const CFuint nodeID) const
  {
    cf_assert(false);
    return (CFNULL);
  }

  /// Get the pastStates corresponding to the given ID
  const RealVector* getPastState(const CFuint stateID) const
  {
    cf_assert(false);
    return (CFNULL);
  }

  /// Get the interNodes corresponding to the given ID
  const RealVector* getInterNode(const CFuint nodeID) const
  {
    cf_assert(false);
    return (CFNULL);
  }

  /// Get the interStates corresponding to the given ID
  const RealVector* getInterState(const CFuint stateID) const
  {
    cf_assert(false);
    return (CFNULL);
  }
  /// Set the extra vars corresponding to the given ID
  RealVector& getExtraNodalValues(const CFuint nodeID)
  {
    //nothing to be done here
    return _dummyVector;
  }

  /// Set the extra vars corresponding to the given ID
  RealVector& getExtraStateValues(const CFuint stateID)
  {
    //nothing to be done here
    return _dummyVector;
  }

  /// Get the node corresponding to the given ID
  /// @post the template parameter allows to return
  ///         const RealVector& or const Node*&
  CFreal* getNode(const CFuint nodeID)
  {
    // return pointer and not reference to RealVector
    // because otherwise code compiled with icc breaks
    return &_nodes[nodeID*_dimension];
  }

  /// Get the state corresponding to the given ID
  /// @post the template parameter allows to return
  ///         const RealVector& or const State*&
  CFreal* getState(const CFuint stateID)
  {
    // return pointer and not reference to RealVector
    // because otherwise code compiled with icc breaks
    return &_states[stateID*_nbEquations];
  }
  
  /// Resize the element-node connectivity
  void resizeElementNode(const std::valarray<CFuint>& nbCols)
  {
    _elementNode.resize(nbCols);
  }

  /// Resize the element-state connectivity
  void resizeElementState(const std::valarray<CFuint>& nbCols)
  {
    _elementState.resize(nbCols);
  }

  /// Get the number of nodes in the given element
  CFuint getNbNodesInElement(const CFuint iElem) const
  {
    return _elementNode.nbCols(iElem);
  }

  /// Get the number of states in the given element
  CFuint getNbStatesInElement(const CFuint iElem) const
  {
    return _elementState.nbCols(iElem);
  }

  /// Set the element-node
  void setElementNode(const CFuint iElem,
          const CFuint iNode,
          const CFuint value)
  {
    _elementNode(iElem, iNode) = value;
  }

  /// Set the element-state
  void setElementState(const CFuint iElem,
           const CFuint iState,
           const CFuint value)
  {
    _elementState(iElem, iState) = value;
  }

  /// Get the element-node
  CFuint getElementNode(const CFuint iElem,
      const CFuint iNode) const
  {
    return _elementNode(iElem, iNode);
  }

  /// Get the element-state
  CFuint getElementState(const CFuint iElem,
        const CFuint iState) const
  {
    return _elementState(iElem, iState);
  }

  /// Copy element node connectivity to the given Table
  void copyElementNodeTo(Common::Table<CFuint>& table)
  {
    table = _elementNode;
  }

  /// Copy element state connectivity to the given Table
  void copyElementStateTo(Common::Table<CFuint>& table)
  {
    table = _elementState;
  }

  /// Get the node list
  Common::SafePtr< std::vector<CFreal> > getNodeList()
  {
    return &_nodes;
  }

  /// Get the state list
  Common::SafePtr< std::vector<CFreal> > getStateList()
  {
    return &_states;
  }

  /// Get the element node connectivity
  Common::SafePtr<Common::Table<CFuint> > getElementNode()
  {
    return &_elementNode;
  }

  /// Get the element state connectivity
  Common::SafePtr<Common::Table<CFuint> > getElementState()
  {
    return &_elementState;
  }

private: //data

  /// storage for the node list
  std::vector<CFreal>  _nodes;
  
  /// storage for the state list
  std::vector<CFreal>  _states;
  
  /// table for the element-node connectivity
  Common::Table<CFuint>    _elementNode;

  /// table for the element-state connectivity
  Common::Table<CFuint>     _elementState;

  ///dummy vector
  RealVector _dummyVector;

}; // end of class CFmeshReaderWriterSource

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshReaderWriterSource_hh
