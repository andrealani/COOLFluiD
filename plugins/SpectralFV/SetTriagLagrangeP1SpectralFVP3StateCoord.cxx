#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/SetTriagLagrangeP1SpectralFVP3StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTriagLagrangeP1SpectralFVP3StateCoord,
               Framework::SetElementStateCoord,
               SpectralFVModule>
SetTriagLagrangeP1SpectralFVP3StateCoord("TriagLagrangeP1SpectralFVP3");

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP1SpectralFVP3StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

void SetTriagLagrangeP1SpectralFVP3StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
