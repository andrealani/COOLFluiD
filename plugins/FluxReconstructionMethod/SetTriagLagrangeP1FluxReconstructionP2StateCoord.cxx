#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetTriagLagrangeP1FluxReconstructionP2StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTriagLagrangeP1FluxReconstructionP2StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetTriagLagrangeP1FluxReconstructionP2StateCoord("TriagLagrangeP1FluxReconstructionP2");

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP1FluxReconstructionP2StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP1FluxReconstructionP2StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
