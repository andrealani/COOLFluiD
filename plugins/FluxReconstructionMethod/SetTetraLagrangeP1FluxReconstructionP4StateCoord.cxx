#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetTetraLagrangeP1FluxReconstructionP4StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTetraLagrangeP1FluxReconstructionP4StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetTetraLagrangeP1FluxReconstructionP4StateCoord("TetraLagrangeP1FluxReconstructionP4");

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP1FluxReconstructionP4StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 35);
  cf_assert(nodes.size() == 4);

  // assign nodes to the states
 
}

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP1FluxReconstructionP4StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
