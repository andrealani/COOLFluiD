#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetHexaLagrangeP1FluxReconstructionP10StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetHexaLagrangeP1FluxReconstructionP10StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetHexaLagrangeP1FluxReconstructionP10StateCoord("HexaLagrangeP1FluxReconstructionP10");

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1FluxReconstructionP10StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 1331);
  cf_assert(nodes.size() == 8);
}

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1FluxReconstructionP10StateCoord::update(const vector<Framework::Node*>& nodes,
                                                     vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
