#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/SetTetraLagrangeP1SpectralFVP2StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTetraLagrangeP1SpectralFVP2StateCoord,
               Framework::SetElementStateCoord,
               SpectralFVModule>
SetTetraLagrangeP1SpectralFVP2StateCoord("TetraLagrangeP1SpectralFVP2");

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP1SpectralFVP2StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP1SpectralFVP2StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
