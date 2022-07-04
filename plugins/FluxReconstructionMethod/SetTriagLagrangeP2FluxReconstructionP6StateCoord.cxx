#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetTriagLagrangeP2FluxReconstructionP6StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTriagLagrangeP2FluxReconstructionP6StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetTriagLagrangeP2FluxReconstructionP6StateCoord("TriagLagrangeP2FluxReconstructionP6");

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP2FluxReconstructionP6StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 28);
  cf_assert(nodes.size() == 6);
}

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP2FluxReconstructionP6StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
