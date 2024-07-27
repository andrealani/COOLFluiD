#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetPrismLagrangeP1FluxReconstructionP4StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetPrismLagrangeP1FluxReconstructionP4StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetPrismLagrangeP1FluxReconstructionP4StateCoord("PrismLagrangeP1FluxReconstructionP4");

//////////////////////////////////////////////////////////////////////////////

void SetPrismLagrangeP1FluxReconstructionP4StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 75);
  cf_assert(nodes.size() == 6);

}

//////////////////////////////////////////////////////////////////////////////

void SetPrismLagrangeP1FluxReconstructionP4StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
