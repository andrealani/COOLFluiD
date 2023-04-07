#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetTetraLagrangeP2FluxReconstructionP4StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTetraLagrangeP2FluxReconstructionP4StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetTetraLagrangeP2FluxReconstructionP4StateCoord("TetraLagrangeP2FluxReconstructionP4");

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP2FluxReconstructionP4StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 35);
  cf_assert(nodes.size() == 10);

  // assign nodes to the states
 
}

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP2FluxReconstructionP4StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
