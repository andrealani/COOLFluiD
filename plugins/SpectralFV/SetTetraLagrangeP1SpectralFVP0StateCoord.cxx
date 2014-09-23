#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/SetTetraLagrangeP1SpectralFVP0StateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetTetraLagrangeP1SpectralFVP0StateCoord,
               Framework::SetElementStateCoord,
               SpectralFVModule>
SetTetraLagrangeP1SpectralFVP0StateCoord("TetraLagrangeP1SpectralFVP0");

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP1SpectralFVP0StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

void SetTetraLagrangeP1SpectralFVP0StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
