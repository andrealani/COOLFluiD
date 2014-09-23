// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/SetTriagLagrangeP1LagrangeP3StateCoord.hh"
#include "Framework/State.hh"
#include "Framework/Node.hh"
#include "Environment/ObjectProvider.hh"
#include "ShapeFunctions/ShapeFunctions.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP3.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTriagLagrangeP1LagrangeP3StateCoord,
           Framework::SetElementStateCoord,  ShapeFunctionsLib>
SetTriagLagrangeP1LagrangeP3StateCoord("TriagLagrangeP1LagrangeP3");

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP1LagrangeP3StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{

  cf_assert(states.size() == 10);
  cf_assert(nodes.size() == 3);

  const CFuint nbNodes = nodes.size();
  const CFuint nbDim = (*(nodes[0])).size();
  cf_assert(nbDim == DIM_2D);

  LagrangeShapeFunctionTriagP3::getStatesMappedCoordinates(_tempMappedCoords);

  //Set the existing nodes to the states at the mesh vertices
  for (CFuint iState = 0; iState < 3; ++iState) {
    states[iState]->setSpaceCoordinates(nodes[iState]);
  }


  for(CFuint iState = 3; iState < states.size() ; ++iState)
  {
    //Compute Coordinates of the state
    LagrangeShapeFunctionTriagP1::computeShapeFunction(_tempMappedCoords[iState], _tempShapeFunctions);

    for(CFuint iDim = 0; iDim < nbDim; ++iDim)
    {
      _tempCoord[iDim] = _tempShapeFunctions[0]*(*(nodes[0]))[iDim];

      for(CFuint iNode = 1; iNode < nbNodes; ++iNode)
      {
        _tempCoord[iDim] += _tempShapeFunctions[iNode]*(*(nodes[iNode]))[iDim];
      }
    }

    //Create node and assign set it to the state
    Framework::Node* node = new Framework::Node(_tempCoord,false);
    states[iState]->setSpaceCoordinates(node);
  }

}

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP1LagrangeP3StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{

  cf_assert(states.size() == 10);
  cf_assert(nodes.size() == 3);

  const CFuint nbNodes = nodes.size();
  const CFuint nbDim = (*(nodes[0])).size();
  cf_assert(nbDim == DIM_2D);

  LagrangeShapeFunctionTriagP3::getStatesMappedCoordinates(_tempMappedCoords);

  //Set the existing nodes to the states at the mesh vertices
  for (CFuint iState = 0; iState < 3; ++iState) {
    states[iState]->setSpaceCoordinates(nodes[iState]);
  }


  for(CFuint iState = 3; iState < states.size() ; ++iState)
  {
    //Compute Coordinates of the state
    LagrangeShapeFunctionTriagP1::computeShapeFunction(_tempMappedCoords[iState], _tempShapeFunctions);

    for(CFuint iDim = 0; iDim < nbDim; ++iDim)
    {
      _tempCoord[iDim] = _tempShapeFunctions[0]*(*(nodes[0]))[iDim];

      for(CFuint iNode = 1; iNode < nbNodes; ++iNode)
      {
        _tempCoord[iDim] += _tempShapeFunctions[iNode]*(*(nodes[iNode]))[iDim];
      }
    }

    //Create node and assign set it to the state
    Framework::Node* node = new Framework::Node(_tempCoord,false);
    states[iState]->setSpaceCoordinates(node);
  }

}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
