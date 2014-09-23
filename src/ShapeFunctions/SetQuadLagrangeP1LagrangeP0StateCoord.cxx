// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/SetQuadLagrangeP1LagrangeP0StateCoord.hh"
#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"
#include "ShapeFunctions/ShapeFunctions.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetQuadLagrangeP1LagrangeP0StateCoord,
               Framework::SetElementStateCoord,
               ShapeFunctionsLib>
setQuadLagrangeP1LagrangeP0StateCoord("QuadLagrangeP1LagrangeP0");

//////////////////////////////////////////////////////////////////////////////

void SetQuadLagrangeP1LagrangeP0StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  setLagrangeP1LagrangeP0StateCoord(nodes, states);
}

//////////////////////////////////////////////////////////////////////////////

void SetQuadLagrangeP1LagrangeP0StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
  updateLagrangeP1LagrangeP0StateCoord(nodes, states);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
