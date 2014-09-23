// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/Node.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Node::Node(bool isOnMesh) :
  RealVector(PhysicalModelStack::getActive()->getDim()),
  IndexedObject<Node>()
{
  _flags.isOnMesh = isOnMesh;
  _flags.parUpdatable = true;

}

//////////////////////////////////////////////////////////////////////////////

Node::Node(const RealVector& inCoord, bool isOnMesh) :
  RealVector(inCoord),
  IndexedObject<Node>()
{
  _flags.isOnMesh = isOnMesh;
  _flags.parUpdatable = true;
}

//////////////////////////////////////////////////////////////////////////////

Node::Node(CFreal * ptr, bool isOnMesh) :
  RealVector(PhysicalModelStack::getActive()->getDim(), ptr),
  IndexedObject<Node>()
{
  cf_assert (ptr != CFNULL);

  _flags.isOnMesh = isOnMesh;
  _flags.parUpdatable = true;
}

//////////////////////////////////////////////////////////////////////////////

Node::Node(const Node& inNode) :
  RealVector(inNode),
  IndexedObject<Node>()
{
  _flags.isOnMesh = inNode._flags.isOnMesh;
  _flags.parUpdatable = inNode._flags.parUpdatable;

}

//////////////////////////////////////////////////////////////////////////////

Node::~Node()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
