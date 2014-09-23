#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SetQuadLagrangeP2SpectralFDP5StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetQuadLagrangeP2SpectralFDP5StateCoord,
               Framework::SetElementStateCoord,
               SpectralFDModule>
SetQuadLagrangeP2SpectralFDP5StateCoord("QuadLagrangeP2SpectralFDP5");

//////////////////////////////////////////////////////////////////////////////

void SetQuadLagrangeP2SpectralFDP5StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                                          vector<Framework::State*>& states)
{
  cf_assert(states.size() == 36);
  cf_assert(nodes.size() == 9);
}

//////////////////////////////////////////////////////////////////////////////

void SetQuadLagrangeP2SpectralFDP5StateCoord::update(const vector<Framework::Node*>& nodes,
                                                     vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
