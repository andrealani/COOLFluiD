#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetPrismLagrangeP2FluxReconstructionP1StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetPrismLagrangeP2FluxReconstructionP1StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetPrismLagrangeP2FluxReconstructionP1StateCoord("PrismLagrangeP2FluxReconstructionP1");

//////////////////////////////////////////////////////////////////////////////

void SetPrismLagrangeP2FluxReconstructionP1StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 6);
  cf_assert(nodes.size() == 18);

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

void SetPrismLagrangeP2FluxReconstructionP1StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
