// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_Node_hh
#define COOLFluiD_Framework_Node_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"

#include "Framework/IndexedObject.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a geometrical Node
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API Node :
    public RealVector,
    public IndexedObject<Node> {

public:

  typedef CFreal GTYPE;

  /// Constructor. Constructs new Node from input coordinates.
  /// @param isOnMesh is a boolean indicating if this node is on the mesh or not.
  explicit Node (bool isOnMesh = false);

  /// Constructor taking preallocated memory
  /// does not initialize the memory!!
  Node (CFreal * data, bool isOnMesh);

  /// Constructor. Constructs new Node from input coordinates.
  /// @param inCoord Input coordinates
  /// @param isOnMesh is a boolean indicating if this node is on the mesh or not.
  Node (const RealVector& inCoord, bool isOnMesh);

  /// Default destructor
  ~Node();

  /// Copy Constructor: Constructs new Node with the inner state equal to inNode.
  /// @param inNode Node to be copied
  Node(const Node& inNode);

  /// Set isOnMesh to true if the node is a meshpoint
  inline void setIsOnMesh(bool isOnMesh);

  /// Set isOwnedByState to true if the node is a meshpoint
  inline void setIsOwnedByState(bool isOnMesh);

  /// Tell if the node is a meshpoint
  inline bool isOnMesh();

  /// Returns true if the node is a ghost node (for the parallellisation)
  bool isParUpdatable() const
  {
      return _flags.parUpdatable;
  };

  /// Set the updatable attribute of the node (for the parallellisation)
  void setParUpdatable(bool up)
  {
    _flags.parUpdatable = up;
  }

  /// Tell if the node is owned by a state
  /// @pre only one owner per state
  inline bool isOwnedByState();

  /// Sets internal Data in the Node to the input data by copying it.
  /// @param inData Input numerical data
  inline void setData(const RealVector& inData);

  /// Gets Data in a RealVector pointer
  /// @return this State as a RealVector
  inline RealVector* getData();

  /// Overloading of the assignment operator "=".
  inline const Node& operator= (const CFreal& value);

  /// Overloading of the assignment operator "=".
  inline const Node& operator= (const RealVector& other);

  /// Overloading of the assignment operator "=".
  inline const Node& operator= (const Node& inNode);

private:

  /// Status bits of the node
  struct
  {
    /// flag indicates if node is updatable in a parallel simulation
    unsigned parUpdatable : 1;

    /// flag indicates if the node belong to a mesh or is a simple point
    unsigned isOnMesh     : 1;

  } _flags;

}; //end of class Node

//////////////////////////////////////////////////////////////////////////////

inline void Node::setIsOnMesh(bool isOnMesh)
{
  _flags.isOnMesh = isOnMesh;
}

//////////////////////////////////////////////////////////////////////////////

inline bool Node::isOnMesh()
{
  return _flags.isOnMesh;
}

//////////////////////////////////////////////////////////////////////////////

inline void Node::setIsOwnedByState(bool isOwnedByState)
{
  _flags.isOnMesh = !isOwnedByState;
}

//////////////////////////////////////////////////////////////////////////////

inline bool Node::isOwnedByState()
{
  return !_flags.isOnMesh;
}

//////////////////////////////////////////////////////////////////////////////

inline void Node::setData(const RealVector& inData)
{
  operator=(inData);
}

//////////////////////////////////////////////////////////////////////////////

inline RealVector * Node::getData()
{
  return this;
}

//////////////////////////////////////////////////////////////////////////////

inline const Node& Node::operator= (const CFreal& value)
{
  RealVector::operator= (value);
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

inline const Node& Node::operator= (const RealVector& other)
{
  cf_assert(&other != this);
  cf_assert(size() == other.size());
  RealVector::operator= (other);
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

inline const Node& Node::operator= (const Node& inNode)
{
  //  cf_assert(&inNode != this);
  cf_assert(size() == inNode.size());
  RealVector::operator= (inNode);
  _flags.isOnMesh = inNode._flags.isOnMesh;
  _flags.parUpdatable = inNode._flags.parUpdatable;
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_Node_hh
