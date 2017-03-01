#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetQuadLagrangeP1FluxReconstructionP1StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetQuadLagrangeP1FluxReconstructionP1StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetQuadLagrangeP1FluxReconstructionP1StateCoord("QuadLagrangeP1FluxReconstructionP1");

//////////////////////////////////////////////////////////////////////////////

void SetQuadLagrangeP1FluxReconstructionP1StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 4);
  cf_assert(nodes.size() == 4);

  // assign nodes to the states
  states[0]->setSpaceCoordinates(nodes[0]);
  states[1]->setSpaceCoordinates(nodes[3]);
  states[2]->setSpaceCoordinates(nodes[1]);
  states[3]->setSpaceCoordinates(nodes[2]);
  
//   for (CFuint i = 0; i < 4; ++i)
//   {
//     states[i]->setParUpdatable(nodes[i]->isParUpdatable());
//   }
  
//   CFLog(VERBOSE, "Owned by state: " << nodes[1]->isOwnedByState() << ", parUpdatable: " << nodes[1]->isParUpdatable() << ", on mesh: " << nodes[1]->isOnMesh() << "\n");
//   CFLog(VERBOSE, "Global ID: " << states[1]->getGlobalID() << ", parUpdatable: " << states[1]->isParUpdatable() << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void SetQuadLagrangeP1FluxReconstructionP1StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
