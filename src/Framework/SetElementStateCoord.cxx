// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/SetElementStateCoord.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

SetElementStateCoord::SetElementStateCoord() :
  Common::OwnedObject(),
  _tempCoord(PhysicalModelStack::getActive()->getDim())
{
}

//////////////////////////////////////////////////////////////////////////////

void SetElementStateCoord::setIsoParamStateCoord(const vector<Node*>& nodes,
                                                 vector<State*>& states)
{
  cf_assert(states.size() == nodes.size());
  for (CFuint iState = 0; iState < states.size(); ++iState) {
    states[iState]->setSpaceCoordinates(nodes[iState]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void SetElementStateCoord::setLagrangeP1LagrangeP0StateCoord
(const vector<Node*>& nodes, vector<State*>& states)
{
  cf_assert(states.size() == 1);
  cf_assert(nodes.size() > 1);

  _tempCoord = 0.0;

  vector<Node*>::const_iterator it;
  for (it = nodes.begin(); it != nodes.end(); ++it) {
    _tempCoord += *(*it);
  }
  _tempCoord /= nodes.size();

  Node* node = new Node(_tempCoord,false);
  states[0]->setSpaceCoordinates(node);
}

//////////////////////////////////////////////////////////////////////////////

void SetElementStateCoord::updateLagrangeP1LagrangeP0StateCoord
(const vector<Node*>& nodes, vector<State*>& states)
{
  cf_assert(states.size() == 1);
  cf_assert(nodes.size() > 1);

  _tempCoord = 0.0;

  vector<Node*>::const_iterator it;
  for (it = nodes.begin(); it != nodes.end(); ++it) {
    _tempCoord += *(*it);
  }
  _tempCoord /= nodes.size();

  Node* node = states[0]->getNodePtr();
  (*node) = _tempCoord;
  states[0]->setSpaceCoordinates(node);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
