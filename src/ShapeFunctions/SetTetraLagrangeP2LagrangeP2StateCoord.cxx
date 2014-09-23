// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/SetTetraLagrangeP2LagrangeP2StateCoord.hh"
#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"
#include "ShapeFunctions/ShapeFunctions.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP2.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTetraP1.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTetraLagrangeP2LagrangeP2StateCoord,
               Framework::SetElementStateCoord,  ShapeFunctionsLib>
SetTetraLagrangeP2LagrangeP2StateCoord("TetraLagrangeP2LagrangeP2");

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP2LagrangeP2StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{

  cf_assert(states.size() == 10);
  cf_assert(nodes.size() == 10);

  const CFuint nbNodes = nodes.size();
  const CFuint nbDim = (*(nodes[0])).size();
  cf_assert(nbDim == DIM_3D);

  LagrangeShapeFunctionTetraP2::getStatesMappedCoordinates(_tempMappedCoords);

  //Set the existing nodes to the states at the mesh vertices
  for (CFuint iState = 0; iState < 4; ++iState) {
    states[iState]->setSpaceCoordinates(nodes[iState]);
// CFout << "State Coordinate[" <<iState<<"]: " << states[iState]->getCoordinates() << "\n";
//   CFout << "State LocalID[" <<iState<<"]: " << states[iState]->getLocalID() << "\n";
  }
//   CFout << "Node Coord[0]: " << (*(nodes[0])) << "\n";
//   CFout << "Node Coord[1]: " << (*(nodes[1])) << "\n";
//   CFout << "Node Coord[2]: " << (*(nodes[2])) << "\n";
//   CFout << "Node Coord[3]: " << (*(nodes[3])) << "\n";

  for(CFuint iState = 4; iState < states.size() ; ++iState)
  {
// CFout << "State Mapped Coordinate[" <<iState<<"]: " << _tempMappedCoords[iState] << "\n";
    //Compute Coordinates of the state
    LagrangeShapeFunctionTetraP1::computeShapeFunction(_tempMappedCoords[iState], _tempShapeFunctions);

    for(CFuint iDim = 0; iDim < nbDim; ++iDim)
    {
      _tempCoord[iDim] = _tempShapeFunctions[0]*(*(nodes[0]))[iDim];

      for(CFuint iNode = 1; iNode < nbNodes; ++iNode)
      {
        _tempCoord[iDim] += _tempShapeFunctions[iNode]*(*(nodes[iNode]))[iDim];
      }
    }

    ///@todo here we are duplicating some of the nodes, try to improve!!
    //Create node and assign set it to the state
    Framework::Node* node = new Framework::Node(_tempCoord,false);
    states[iState]->setSpaceCoordinates(node);
// CFout << "State Coordinate[" <<iState<<"]: " << states[iState]->getCoordinates() << "\n";
// CFout << "State LocalID[" <<iState<<"]: " << states[iState]->getLocalID() << "\n";

  }

}

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP2LagrangeP2StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
  cf_assert(states.size() == 10);
  cf_assert(nodes.size() == 10);

  const CFuint nbNodes = nodes.size();
  const CFuint nbDim = (*(nodes[0])).size();
  cf_assert(nbDim == DIM_3D);

  LagrangeShapeFunctionTetraP2::getStatesMappedCoordinates(_tempMappedCoords);

  //Set the existing nodes to the states at the mesh vertices
  for (CFuint iState = 0; iState < 4; ++iState) {
    states[iState]->setSpaceCoordinates(nodes[iState]);
//  CFout << "State Coordinate[" <<iState<<"]: " << states[iState]->getCoordinates() << "\n";
//   CFout << "State LocalID[" <<iState<<"]: " << states[iState]->getLocalID() << "\n";
  }
//   CFout << "Node Coord[0]: " << (*(nodes[0])) << "\n";
//   CFout << "Node Coord[1]: " << (*(nodes[1])) << "\n";
//   CFout << "Node Coord[2]: " << (*(nodes[2])) << "\n";
//   CFout << "Node Coord[3]: " << (*(nodes[3])) << "\n";

  for(CFuint iState = 4; iState < states.size() ; ++iState)
  {
//  CFout << "State Mapped Coordinate[" <<iState<<"]: " << _tempMappedCoords[iState] << "\n";
    //Compute Coordinates of the state
    LagrangeShapeFunctionTetraP1::computeShapeFunction(_tempMappedCoords[iState], _tempShapeFunctions);

    for(CFuint iDim = 0; iDim < nbDim; ++iDim)
    {
      _tempCoord[iDim] = _tempShapeFunctions[0]*(*(nodes[0]))[iDim];

      for(CFuint iNode = 1; iNode < nbNodes; ++iNode)
      {
        _tempCoord[iDim] += _tempShapeFunctions[iNode]*(*(nodes[iNode]))[iDim];
      }
    }

    //Assign the node value to the node associated to the state
    Framework::Node* node = states[iState]->getNodePtr();
    (*node) = _tempCoord;
    states[iState]->setSpaceCoordinates(node);
//  CFout << "State Coordinate[" <<iState<<"]: " << states[iState]->getCoordinates() << "\n";
//   CFout << "State LocalID[" <<iState<<"]: " << states[iState]->getLocalID() << "\n";

  }

}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
