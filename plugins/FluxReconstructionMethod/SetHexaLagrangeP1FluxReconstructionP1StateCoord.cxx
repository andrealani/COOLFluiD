#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetHexaLagrangeP1FluxReconstructionP1StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetHexaLagrangeP1FluxReconstructionP1StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetHexaLagrangeP1FluxReconstructionP1StateCoord("HexaLagrangeP1FluxReconstructionP1");

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1FluxReconstructionP1StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 8);
  cf_assert(nodes.size() == 8);

  // assign nodes to the states
//  states[0]->setSpaceCoordinates(nodes[0]);
//  states[1]->setSpaceCoordinates(nodes[4]);
//  states[2]->setSpaceCoordinates(nodes[3]);
//  states[3]->setSpaceCoordinates(nodes[7]);
//  states[4]->setSpaceCoordinates(nodes[1]);
//  states[5]->setSpaceCoordinates(nodes[5]);
//  states[6]->setSpaceCoordinates(nodes[2]);
//  states[7]->setSpaceCoordinates(nodes[6]);
}

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1FluxReconstructionP1StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
