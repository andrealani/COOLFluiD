#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SetQuadLagrangeP1SpectralFDP1StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetQuadLagrangeP1SpectralFDP1StateCoord,
               Framework::SetElementStateCoord,
               SpectralFDModule>
SetQuadLagrangeP1SpectralFDP1StateCoord("QuadLagrangeP1SpectralFDP1");

//////////////////////////////////////////////////////////////////////////////

void SetQuadLagrangeP1SpectralFDP1StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 4);
  cf_assert(nodes.size() == 4);

  // assign nodes to the states
  states[0]->setSpaceCoordinates(nodes[0]);
  states[1]->setSpaceCoordinates(nodes[3]);
  states[2]->setSpaceCoordinates(nodes[1]);
  states[3]->setSpaceCoordinates(nodes[2]);
}

//////////////////////////////////////////////////////////////////////////////

void SetQuadLagrangeP1SpectralFDP1StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
