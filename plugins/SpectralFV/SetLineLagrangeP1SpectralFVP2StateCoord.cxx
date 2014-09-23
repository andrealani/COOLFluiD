#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/SetLineLagrangeP1SpectralFVP2StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetLineLagrangeP1SpectralFVP2StateCoord,
               Framework::SetElementStateCoord,
               SpectralFVModule>
SetLineLagrangeP1SpectralFVP2StateCoord("LineLagrangeP1SpectralFVP2");

//////////////////////////////////////////////////////////////////////////////

void SetLineLagrangeP1SpectralFVP2StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

void SetLineLagrangeP1SpectralFVP2StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
