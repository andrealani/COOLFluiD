#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SetHexaLagrangeP1SpectralFDP2StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetHexaLagrangeP1SpectralFDP2StateCoord,
               Framework::SetElementStateCoord,
               SpectralFDModule>
SetHexaLagrangeP1SpectralFDP2StateCoord("HexaLagrangeP1SpectralFDP2");

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1SpectralFDP2StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 27);
  cf_assert(nodes.size() == 8);

  states[0]->setSpaceCoordinates(nodes[0]);
  states[18]->setSpaceCoordinates(nodes[1]);
  states[24]->setSpaceCoordinates(nodes[2]);
  states[6]->setSpaceCoordinates(nodes[3]);
  states[2]->setSpaceCoordinates(nodes[4]);
  states[20]->setSpaceCoordinates(nodes[5]);
  states[26]->setSpaceCoordinates(nodes[6]);
  states[8]->setSpaceCoordinates(nodes[7]);

  // edges
  _tempCoord = 0.5*(*nodes[0] + *nodes[4]);
  states[1]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  _tempCoord = 0.5*(*nodes[0] + *nodes[3]);
  states[3]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  _tempCoord = 0.5*(*nodes[4] + *nodes[7]);
  states[5]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  _tempCoord = 0.5*(*nodes[3] + *nodes[7]);
  states[7]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));

  _tempCoord = 0.5*(*nodes[1] + *nodes[5]);
  states[19]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  _tempCoord = 0.5*(*nodes[1] + *nodes[2]);
  states[21]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  _tempCoord = 0.5*(*nodes[5] + *nodes[6]);
  states[23]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  _tempCoord = 0.5*(*nodes[2] + *nodes[6]);
  states[25]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));

  _tempCoord = 0.5*(*nodes[0] + *nodes[1]);
  states[9]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  _tempCoord = 0.5*(*nodes[4] + *nodes[5]);
  states[11]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  _tempCoord = 0.5*(*nodes[2] + *nodes[3]);
  states[15]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  _tempCoord = 0.5*(*nodes[6] + *nodes[7]);
  states[17]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));

  // faces
  _tempCoord = 0.25*(*nodes[0] + *nodes[3] + *nodes[4] + *nodes[7]);
  states[4]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  _tempCoord = 0.25*(*nodes[1] + *nodes[2] + *nodes[5] + *nodes[6]);
  states[22]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));

  _tempCoord = 0.25*(*nodes[0] + *nodes[1] + *nodes[2] + *nodes[3]);
  states[12]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  _tempCoord = 0.25*(*nodes[4] + *nodes[5] + *nodes[6] + *nodes[7]);
  states[14]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));

  _tempCoord = 0.25*(*nodes[0] + *nodes[1] + *nodes[4] + *nodes[5]);
  states[10]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  _tempCoord = 0.25*(*nodes[2] + *nodes[3] + *nodes[6] + *nodes[7]);
  states[16]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));

  // center
  _tempCoord = ((*nodes[0]) + (*nodes[1]) + (*nodes[2]) + (*nodes[3]) +
                (*nodes[4]) + (*nodes[5]) + (*nodes[6]) + (*nodes[7]))/8.0;
  states[13]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
  /// @todo some nodes on the faces are duplicated this way, try to avoid
}

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1SpectralFDP2StateCoord::update(const vector<Framework::Node*>& nodes,
                                                     vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
