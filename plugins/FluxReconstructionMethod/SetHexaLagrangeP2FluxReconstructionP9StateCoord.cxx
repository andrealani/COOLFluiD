#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetHexaLagrangeP2FluxReconstructionP9StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetHexaLagrangeP2FluxReconstructionP9StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetHexaLagrangeP2FluxReconstructionP9StateCoord("HexaLagrangeP2FluxReconstructionP9");

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP2FluxReconstructionP9StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 1000);
  cf_assert(nodes.size() == 27);
}

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP2FluxReconstructionP9StateCoord::update(const vector<Framework::Node*>& nodes,
                                                     vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
