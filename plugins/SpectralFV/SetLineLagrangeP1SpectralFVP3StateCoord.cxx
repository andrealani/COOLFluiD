#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/SetLineLagrangeP1SpectralFVP3StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetLineLagrangeP1SpectralFVP3StateCoord,
               Framework::SetElementStateCoord,
               SpectralFVModule>
SetLineLagrangeP1SpectralFVP3StateCoord("LineLagrangeP1SpectralFVP3");

//////////////////////////////////////////////////////////////////////////////

void SetLineLagrangeP1SpectralFVP3StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

void SetLineLagrangeP1SpectralFVP3StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
