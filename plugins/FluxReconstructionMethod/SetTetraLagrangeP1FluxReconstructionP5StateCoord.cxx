#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetTetraLagrangeP1FluxReconstructionP5StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTetraLagrangeP1FluxReconstructionP5StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetTetraLagrangeP1FluxReconstructionP5StateCoord("TetraLagrangeP1FluxReconstructionP5");

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP1FluxReconstructionP5StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 56);
  cf_assert(nodes.size() == 4);

  // assign nodes to the states
 
}

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP1FluxReconstructionP5StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
