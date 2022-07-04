#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetTriagLagrangeP1FluxReconstructionP4StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTriagLagrangeP1FluxReconstructionP4StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetTriagLagrangeP1FluxReconstructionP4StateCoord("TriagLagrangeP1FluxReconstructionP4");

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP1FluxReconstructionP4StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP1FluxReconstructionP4StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
