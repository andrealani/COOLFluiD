#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SetHexaLagrangeP2SpectralFDP4StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetHexaLagrangeP2SpectralFDP4StateCoord,
               Framework::SetElementStateCoord,
               SpectralFDModule>
SetHexaLagrangeP2SpectralFDP4StateCoord("HexaLagrangeP2SpectralFDP4");

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP2SpectralFDP4StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 125);
  cf_assert(nodes.size() == 27);
}

//////////////////////////////////////////////////////////////////////////////

void SetHexaLagrangeP2SpectralFDP4StateCoord::update(const vector<Framework::Node*>& nodes,
                                                     vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
