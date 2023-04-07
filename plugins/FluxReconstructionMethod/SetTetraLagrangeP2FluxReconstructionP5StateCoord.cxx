#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetTetraLagrangeP2FluxReconstructionP5StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTetraLagrangeP2FluxReconstructionP5StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetTetraLagrangeP2FluxReconstructionP5StateCoord("TetraLagrangeP2FluxReconstructionP5");

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP2FluxReconstructionP5StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 56);
  cf_assert(nodes.size() == 10);

  // assign nodes to the states
 
}

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP2FluxReconstructionP5StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
