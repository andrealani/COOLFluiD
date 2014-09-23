// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/SetTriagLagrangeP2LagrangeP1StateCoord.hh"
#include "Framework/State.hh"
#include "Framework/Node.hh"
#include "Environment/ObjectProvider.hh"
#include "ShapeFunctions/ShapeFunctions.hh"
#include "ShapeFunctions/LagrangeShapeFunctionTriagP2.hh"
//#include "ShapeFunctions/LagrangeShapeFunctionTriagP1.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTriagLagrangeP2LagrangeP1StateCoord,
           Framework::SetElementStateCoord,  ShapeFunctionsLib>
SetTriagLagrangeP2LagrangeP1StateCoord("TriagLagrangeP2LagrangeP1");

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP2LagrangeP1StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{

  cf_assert(states.size() == 3);
  cf_assert(nodes.size() == 6);

  cf_assert((*(nodes[0])).size() == DIM_2D);

  //LagrangeShapeFunctionTriagP2::getStatesMappedCoordinates(_tempMappedCoords);

  //Set the existing nodes to the states at the mesh vertices
  for (CFuint iState = 0; iState < 3; ++iState) {
    states[iState]->setSpaceCoordinates(nodes[iState]);
  }


}

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP2LagrangeP1StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{

  cf_assert(states.size() == 3);
  cf_assert(nodes.size() == 6);

  cf_assert((*(nodes[0])).size() == DIM_2D);


  //Set the existing nodes to the states at the mesh vertices
  for (CFuint iState = 0; iState < 3; ++iState) {
    states[iState]->setSpaceCoordinates(nodes[iState]);
  }



}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
