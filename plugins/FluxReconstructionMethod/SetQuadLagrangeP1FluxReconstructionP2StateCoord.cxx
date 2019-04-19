#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetQuadLagrangeP1FluxReconstructionP2StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetQuadLagrangeP1FluxReconstructionP2StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetQuadLagrangeP1FluxReconstructionP2StateCoord("QuadLagrangeP1FluxReconstructionP2");

//////////////////////////////////////////////////////////////////////////////

void SetQuadLagrangeP1FluxReconstructionP2StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 9);
  cf_assert(nodes.size() == 4);

//  states[0]->setSpaceCoordinates(nodes[0]);
//  states[6]->setSpaceCoordinates(nodes[1]);
//  states[2]->setSpaceCoordinates(nodes[3]);
//  states[8]->setSpaceCoordinates(nodes[2]);
//
//  _tempCoord = 0.5*(*nodes[0] + *nodes[3]);
//  states[1]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
//  _tempCoord = 0.5*(*nodes[0] + *nodes[1]);
//  states[3]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
//  _tempCoord = 0.25*(*nodes[0] + *nodes[1] + *nodes[2] + *nodes[3]);
//  states[4]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
//  _tempCoord = 0.5*(*nodes[2] + *nodes[3]);
//  states[5]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
//  _tempCoord = 0.5*(*nodes[1] + *nodes[2]);
//  states[7]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
//  /// @todo some nodes on the faces are duplicated this way, try to avoid
}

//////////////////////////////////////////////////////////////////////////////

void SetQuadLagrangeP1FluxReconstructionP2StateCoord::update(const vector<Framework::Node*>& nodes,
                                                     vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
