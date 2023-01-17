#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetTetraLagrangeP2FluxReconstructionP1StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTetraLagrangeP2FluxReconstructionP1StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetTetraLagrangeP2FluxReconstructionP1StateCoord("TetraLagrangeP2FluxReconstructionP1");

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP2FluxReconstructionP1StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 4);
  cf_assert(nodes.size() == 10);

  // assign nodes to the states

}

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP2FluxReconstructionP1StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
