#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SetHexaLagrangeP1SpectralFDP0StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetHexaLagrangeP1SpectralFDP0StateCoord,
               Framework::SetElementStateCoord,
               SpectralFDModule>
SetHexaLagrangeP1SpectralFDP0StateCoord("HexaLagrangeP1SpectralFDP0");

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1SpectralFDP0StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                                          vector<Framework::State*>& states)
{
  cf_assert(states.size() == 1);
  cf_assert(nodes.size() == 8);

  _tempCoord = ((*nodes[0]) + (*nodes[1]) + (*nodes[2]) + (*nodes[3]) +
                (*nodes[4]) + (*nodes[5]) + (*nodes[6]) + (*nodes[7]))/8.0;

  // create node and assign to the state
  Framework::Node* node = new Framework::Node(_tempCoord,false);
  states[0]->setSpaceCoordinates(node);
}

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1SpectralFDP0StateCoord::update(const vector<Framework::Node*>& nodes,
                                                     vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
