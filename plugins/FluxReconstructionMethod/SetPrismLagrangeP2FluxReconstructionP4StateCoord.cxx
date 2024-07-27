#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetPrismLagrangeP2FluxReconstructionP4StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetPrismLagrangeP2FluxReconstructionP4StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetPrismLagrangeP2FluxReconstructionP4StateCoord("PrismLagrangeP2FluxReconstructionP4");

//////////////////////////////////////////////////////////////////////////////

void SetPrismLagrangeP2FluxReconstructionP4StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 75);
  cf_assert(nodes.size() == 18);

}

//////////////////////////////////////////////////////////////////////////////

void SetPrismLagrangeP2FluxReconstructionP4StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
