#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetPrismLagrangeP1FluxReconstructionP0StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetPrismLagrangeP1FluxReconstructionP0StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetPrismLagrangeP1FluxReconstructionP0StateCoord("PrismLagrangeP1FluxReconstructionP0");

//////////////////////////////////////////////////////////////////////////////

void SetPrismLagrangeP1FluxReconstructionP0StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                                          vector<Framework::State*>& states)
{
  cf_assert(states.size() == 1);
  cf_assert(nodes.size() == 6);

  _tempCoord = ((*nodes[0]) + (*nodes[1]) + (*nodes[2]) + (*nodes[3]) +
                (*nodes[4]) + (*nodes[5]))/6.0;

  // create node and assign to the state
  Framework::Node* node = new Framework::Node(_tempCoord,false);
  states[0]->setSpaceCoordinates(node);
}

//////////////////////////////////////////////////////////////////////////////

void SetPrismLagrangeP1FluxReconstructionP0StateCoord::update(const vector<Framework::Node*>& nodes,
                                                     vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
