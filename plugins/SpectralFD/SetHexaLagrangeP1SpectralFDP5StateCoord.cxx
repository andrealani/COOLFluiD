#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SetHexaLagrangeP1SpectralFDP5StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetHexaLagrangeP1SpectralFDP5StateCoord,
               Framework::SetElementStateCoord,
               SpectralFDModule>
SetHexaLagrangeP1SpectralFDP5StateCoord("HexaLagrangeP1SpectralFDP5");

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1SpectralFDP5StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 216);
  cf_assert(nodes.size() == 8);
}

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP1SpectralFDP5StateCoord::update(const vector<Framework::Node*>& nodes,
                                                     vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
