#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SetHexaLagrangeP2SpectralFDP2StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetHexaLagrangeP2SpectralFDP2StateCoord,
               Framework::SetElementStateCoord,
               SpectralFDModule>
SetHexaLagrangeP2SpectralFDP2StateCoord("HexaLagrangeP2SpectralFDP2");

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP2SpectralFDP2StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 27);
  cf_assert(nodes.size() == 27);

  states[ 0]->setSpaceCoordinates(nodes[ 0]);
  states[18]->setSpaceCoordinates(nodes[ 1]);
  states[24]->setSpaceCoordinates(nodes[ 2]);
  states[ 6]->setSpaceCoordinates(nodes[ 3]);
  states[ 2]->setSpaceCoordinates(nodes[ 4]);
  states[20]->setSpaceCoordinates(nodes[ 5]);
  states[26]->setSpaceCoordinates(nodes[ 6]);
  states[ 8]->setSpaceCoordinates(nodes[ 7]);

  // edges
  states[ 1]->setSpaceCoordinates(nodes[13]);
  states[ 3]->setSpaceCoordinates(nodes[11]);
  states[ 5]->setSpaceCoordinates(nodes[25]);
  states[ 7]->setSpaceCoordinates(nodes[19]);

  states[19]->setSpaceCoordinates(nodes[15]);
  states[21]->setSpaceCoordinates(nodes[ 9]);
  states[23]->setSpaceCoordinates(nodes[23]);
  states[25]->setSpaceCoordinates(nodes[17]);

  states[ 9]->setSpaceCoordinates(nodes[ 8]);
  states[11]->setSpaceCoordinates(nodes[22]);
  states[15]->setSpaceCoordinates(nodes[10]);
  states[17]->setSpaceCoordinates(nodes[24]);

  // faces
  states[ 4]->setSpaceCoordinates(nodes[20]);
  states[22]->setSpaceCoordinates(nodes[16]);

  states[12]->setSpaceCoordinates(nodes[12]);
  states[14]->setSpaceCoordinates(nodes[26]);

  states[10]->setSpaceCoordinates(nodes[14]);
  states[16]->setSpaceCoordinates(nodes[18]);

  // center
  states[13]->setSpaceCoordinates(nodes[21]);
}

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP2SpectralFDP2StateCoord::update(const vector<Framework::Node*>& nodes,
                                                     vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
