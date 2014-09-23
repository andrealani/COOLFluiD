#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/SetTriagLagrangeP1SpectralFVP2StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTriagLagrangeP1SpectralFVP2StateCoord,
               Framework::SetElementStateCoord,
               SpectralFVModule>
SetTriagLagrangeP1SpectralFVP2StateCoord("TriagLagrangeP1SpectralFVP2");

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP1SpectralFVP2StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP1SpectralFVP2StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
