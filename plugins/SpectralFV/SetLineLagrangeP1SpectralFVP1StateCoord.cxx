#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/SetLineLagrangeP1SpectralFVP1StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetLineLagrangeP1SpectralFVP1StateCoord,
               Framework::SetElementStateCoord,
               SpectralFVModule>
SetLineLagrangeP1SpectralFVP1StateCoord("LineLagrangeP1SpectralFVP1");

//////////////////////////////////////////////////////////////////////////////

void SetLineLagrangeP1SpectralFVP1StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

void SetLineLagrangeP1SpectralFVP1StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
