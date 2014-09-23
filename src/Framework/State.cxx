// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

State::State() :
  RealVector(PhysicalModelStack::getActive()->getNbEq()),
  IndexedObject<State>(),
  _coordinate(CFNULL)
{
  _flags.ghost = false;
  _flags.parUpdatable = false;
  _flags.generic = false;
}

//////////////////////////////////////////////////////////////////////////////

State::State(CFreal * ptr) :
  RealVector(PhysicalModelStack::getActive()->getNbEq(), ptr),
  IndexedObject<State>(),
  _coordinate(CFNULL)
{
  cf_assert ( ptr != CFNULL );
  _flags.ghost = false;
  _flags.parUpdatable = false;
  _flags.generic = false;
}

//////////////////////////////////////////////////////////////////////////////

State::State(const RealVector& data,  const bool isGhost) :
  RealVector(data),
  IndexedObject<State>(),
  _coordinate(CFNULL)
{
  _flags.ghost = isGhost;
  _flags.parUpdatable = false;
  _flags.generic = false;
}

//////////////////////////////////////////////////////////////////////////////

State::State(const State& inState) :
  RealVector(inState),                  // copy the data
  IndexedObject<State>(inState)         // copy the indexes
{
  // copy the coordinate into a new Node
  const bool ownedByMesh = false;
  if (inState._coordinate != CFNULL) {
    _coordinate = new Node(*inState._coordinate,ownedByMesh);
  }
  else {
    _coordinate = CFNULL;
  }

  // copy the flags
  _flags = inState._flags;
}

//////////////////////////////////////////////////////////////////////////////

State::~State()
{
  if (_coordinate != CFNULL){
    if(_coordinate->isOwnedByState())
    {
      deletePtr(_coordinate);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

bool State::isValid() const
{
  return PhysicalModelStack::getActive()->validate(*this);
}

//////////////////////////////////////////////////////////////////////////////

void State::printCoord(std::ostream& out) const
{
  cf_assert(_coordinate != CFNULL);
  cf_assert(_coordinate->size() == PhysicalModelStack::getActive()->getDim());
  for (CFuint j = 0; j < PhysicalModelStack::getActive()->getDim(); ++j) {
    out.precision(12);
    out << (*_coordinate)[j] << " " ;
  }
}

//////////////////////////////////////////////////////////////////////////////

void State::dump(std::ostream& out) const
{
  out << "State Dump: " << this << " -----------------------\n";
  out << "Value: " << *this << "\n";
  if(_coordinate != CFNULL) {
    out << "Coordinates: "; printCoord(out); out << "\n";
  } else { CFout << "No Coordinates\n";}
  out << "ParUpdatable: " << _flags.parUpdatable << "\n";
  out << "Ghost: " << _flags.ghost << "\n";
  out << "Generic: " << _flags.generic << "\n";
  if(hasLocalID())  out << "LocalID: " << getLocalID() << "\n";
  else out << "No LocalID\n";
  if(hasGlobalID()) out << "GlobalID: " << getGlobalID() << "\n";
  else out << "No GlobalID\n";
  out << "----------------------------------------------" << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void State::readCoord(std::istream& in)
{
  cf_assert(_coordinate != CFNULL);
  cf_assert(_coordinate->size() == PhysicalModelStack::getActive()->getDim());
  for (CFuint j = 0; j < PhysicalModelStack::getActive()->getDim(); ++j) {
    in >> (*_coordinate)[j];
  }
}

//////////////////////////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream& out, const State& state)
{
  for (CFuint j = 0; j < state.size(); ++j) {
    out.precision(12);
    out << state[j] << " " ;
  }
  return out;
}

//////////////////////////////////////////////////////////////////////////////

std::istream& operator>>(std::istream& in, State& state)
{
  for (CFuint j = 0; j < state.size(); ++j) {
    in >> state[j];
  }
  return in;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
