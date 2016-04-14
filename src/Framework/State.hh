// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_State_hh
#define COOLFluiD_Framework_State_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a solution State.
/// @author Andrea Lani
/// @author Tiago Quintino
/// @author Dries Kimpe
class Framework_API State :
    public RealVector,
    public IndexedObject<State>
{
public:

  typedef CFreal GTYPE;

  /// Default constructor
  State();

  /// Constructor taking preallocated memory
  /// (does not initialise the memory!)

  State (CFreal* data);

  /// Constructor
  /// @param data initialising state
  /// @param  flag that tells if the current state is a ghost one
  State(const RealVector& data, const bool isGhost = false);

  /// Destructor
  ~State();

  /// Copy Constructor: Constructs new State with the inner state equal to inState.
  /// @param inNode Node to be copied
  State(const State& inState);

  /// Sets internal Data in the State to the input data by copying it.
  /// Only copies the data.
  /// @param inData Input numerical data
  inline void copyData(const RealVector& inData);

  /// Clones all the data, flags, coordinates and indexes into this object
  /// @param other State to be fully cloned
  inline void clone(const State& other);

  /// Gets Data in a RealVector pointer
  /// @return this State as a RealVector
  inline RealVector* getData();

  /// Overloading of the assignment operator "=".
  /// @param value constant value
  inline const State& operator= (const CFreal& value);

  /// Overloading of the assignment operator "=".
  /// Only copies the data.
  /// @param other RealVector with the data to be copied
  inline const State& operator= (const RealVector& other);

  /// Overloading of the assignment operator "=".
  /// Copies the data, the coordinates into a new Node, but does not copy the flags neither the indexes.
  /// @param other State with the data, and coordinates to be copied
  inline const State& operator= (const State& other);

  /// Set the space coordinates of the state
  /// @pre coordinate != CFNULL
  /// @param coordinate space coordinates of the state
  void setSpaceCoordinates(Node* const coordinate)
  {
    cf_assert(coordinate != CFNULL);
    _coordinate = coordinate;
  }

  /// Reset the space coordinates of the state
  /// @post coordinate == CFNULL
  void resetSpaceCoordinates()
  {
    _coordinate = CFNULL;
  }

  /// Get the coordinates
  Node& getCoordinates() const
  {
    cf_assert(_coordinate != CFNULL);
    return *_coordinate;
  }

  /// Get the pointer to the node
  Node* getNodePtr() const
  {
    cf_assert(_coordinate != CFNULL);
    return _coordinate;
  }


  /// Tells if the state is a ghost one
  bool isGhost() const
  {
    return _flags.ghost;
  }

  /// Set the ghost flag on the State
  void setGhost (bool ghost)
  {
    _flags.ghost = ghost;
  }

  /// Check the validity of the current data held by this state.
  /// @return true if the data is invalid.
  bool isValid() const;

  /// Prints the coordinates in this state to the supplied ostream
  void printCoord(std::ostream& out) const;

  /// Reads the coordinates in this state from the supplied istream
  void readCoord(std::istream& in);

  /// Dumps all the data in the state
  void dump(std::ostream& out) const;

  /// Overloading of the stream operator "<<" for the output
  Framework_API friend std::ostream& operator<<(std::ostream& out, const State& state);

  /// Overloading of the stream operator ">>" for the input
  Framework_API friend std::istream& operator>>(std::istream& in, State& state);

  /// Returns true if the state is a ghost state (for the parallellisation)
  bool isParUpdatable() const
  {
      return _flags.parUpdatable;
  };

  /// Set the updatable attribute of the state (for the parallellisation)
  void setParUpdatable(bool up)
  {
    _flags.parUpdatable = up;
  }

  /// Set generic Flag on
  void setGenericFlagOn()
  {
    _flags.generic = true;
  }

  /// Set generic Flag off
  void setGenericFlagOff()
  {
    _flags.generic = false;
  }

  /// Tell if the State is flagged on
  bool isGenericFlagOn() const
  {
    return _flags.generic;
  }

private:

  // the position of this State
  Node*      _coordinate;

  /// Status bits of the state
  struct
  {
    /// flag indicates if state is updatable in a parallel simulation
    unsigned parUpdatable : 1;
    /// flag indicates if state is a ghost state, used in some SpaceMethod's
    unsigned ghost        : 1;
    /// flag used in generic algorithms
    unsigned generic      : 1;

  } _flags;

}; // end of class State

//////////////////////////////////////////////////////////////////////////////


inline void State::copyData(const RealVector& inData)
{
  operator=(inData);
}

//////////////////////////////////////////////////////////////////////////////

inline RealVector* State::getData()
{
  return this;
}

//////////////////////////////////////////////////////////////////////////////

inline const State& State::operator= (const CFreal& value)
{
  RealVector::operator= (value);
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

inline const State& State::operator= (const RealVector& other)
{
  cf_assert(&other != this);
  cf_assert(size() == other.size());
  RealVector::operator= (other);
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

inline const State& State::operator= (const State& inState)
{
  // cf_assert(&inState != this);
  cf_assert(size() == inState.size());

  // copys all data and coordinates using assignment operator
  RealVector::operator= (inState);
  if (inState._coordinate != CFNULL) {
    if (_coordinate == CFNULL) {
      _coordinate = new Node(false);
    }
    cf_assert(_coordinate->size() == inState._coordinate->size());
    *(_coordinate) = *(inState._coordinate);
  }

  return *this;
}

//////////////////////////////////////////////////////////////////////////////

inline void State::clone(const State& inState)
{
  // copys all data and coordinates using assignment operator
  *this = inState;
  // copy the flags
  _flags = inState._flags;
  // copy the indexes
  IndexedObject<State>::operator= (inState);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_State_hh
