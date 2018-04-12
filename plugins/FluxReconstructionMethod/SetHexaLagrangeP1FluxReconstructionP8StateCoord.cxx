#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetHexaLagrangeP1FluxReconstructionP8StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetHexaLagrangeP1FluxReconstructionP8StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetHexaLagrangeP1FluxReconstructionP8StateCoord("HexaLagrangeP1FluxReconstructionP8");

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1FluxReconstructionP8StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 729);
  cf_assert(nodes.size() == 8);
}

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1FluxReconstructionP8StateCoord::update(const vector<Framework::Node*>& nodes,
                                                     vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
