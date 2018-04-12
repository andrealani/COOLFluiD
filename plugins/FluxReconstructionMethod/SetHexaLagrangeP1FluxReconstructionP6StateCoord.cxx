#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetHexaLagrangeP1FluxReconstructionP6StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetHexaLagrangeP1FluxReconstructionP6StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetHexaLagrangeP1FluxReconstructionP6StateCoord("HexaLagrangeP1FluxReconstructionP6");

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1FluxReconstructionP6StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 343);
  cf_assert(nodes.size() == 8);
}

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1FluxReconstructionP6StateCoord::update(const vector<Framework::Node*>& nodes,
                                                     vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
