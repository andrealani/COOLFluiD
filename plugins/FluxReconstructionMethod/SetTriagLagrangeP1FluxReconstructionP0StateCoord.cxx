#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetTriagLagrangeP1FluxReconstructionP0StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTriagLagrangeP1FluxReconstructionP0StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetTriagLagrangeP1FluxReconstructionP0StateCoord("TriagLagrangeP1FluxReconstructionP0");

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP1FluxReconstructionP0StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP1FluxReconstructionP0StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
