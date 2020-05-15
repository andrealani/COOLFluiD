#include "GReKO.hh"
#include "NavierStokes2DGReKLogOPuvt.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/ShouldNotBeHereException.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::KOmega;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes2DGReKLogOPuvt, DiffusiveVarSet,
			    GReKOModule, 2>
ns2DGReKLogOPuvtProvider("NavierStokes2DGReKLogOPuvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGReKLogOPuvt::NavierStokes2DGReKLogOPuvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokesKLogOmegaPvt<NavierStokesKLogOmegaVarSet<NavierStokes2DVarSet, 0> >(name, model),
  _unperturbedFluxGa(),
  _unperturbedFluxRe()
{
  const CFuint nbTurbEquations = _eulerModel->getNbScalarVars(0);
  
  vector<std::string> names(4 + nbTurbEquations);
  names[0] = "p";
  names[1] = "u";
  names[2] = "v";
  names[3] = "T";
  
  // Names for turbulent variables
  names[4] = "K";
  names[5] = "Omega";
  
  // Names for transition onset variables
  names[6] = "Ga";
  names[7] = "Re";
  
  setVarNames(names);
  setModelCoefficients();
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGReKLogOPuvt::~NavierStokes2DGReKLogOPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& NavierStokes2DGReKLogOPuvt::getFlux(const RealVector& values,
                                              const vector<RealVector*>& gradients,
                                              const RealVector& normal,
                                              const CFreal& radius)
{
  // here nsData is filled in
  computeTransportProperties(values, gradients, normal);
  computeStressTensor(values, gradients, radius);
  
  const RealVector& gradT     = *gradients[3];
  const RealVector& gradK     = *gradients[4];
  const RealVector& gradOmega = *gradients[5];
  const RealVector& gradGa    = *gradients[6];
  const RealVector& gradRe    = *gradients[7];

  const CFreal avU = values[1];
  const CFreal avV = values[2];

  RealVector& nsData = getModel().getPhysicalData();
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];

  const CFreal qFlux = -getModel().getCoeffQ()*nsData[NSTurbTerm::LAMBDA]*
    (gradT[XX]*nx + gradT[YY]*ny);
  
  const CFreal tauXX = _tau(XX, XX);
  const CFreal tauXY = _tau(XX, YY);
  const CFreal tauYY = _tau(YY, YY);
  _flux[1] = tauXX*nx + tauXY*ny;
  _flux[2] = tauXY*nx + tauYY*ny;
  _flux[3] = (tauXX*avU + tauXY*avV)*nx + (tauXY*avU + tauYY*avV)*ny - qFlux;
  
  //diff flux corresponding to K
  computeBlendingCoefFromGradientVars(values,gradK, gradOmega);
  CFreal coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  CFreal coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]*getSigmaK();

  _flux[4] = (coeffTauMu1 + coeffTauMu2)*(gradK[XX]*nx + gradK[YY]*ny);
  //diff flux corresponding to Omega
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]*getSigmaOmega();

  _flux[5] = (coeffTauMu1 + coeffTauMu2)*(gradOmega[XX]*nx + gradOmega[YY]*ny);
  
  //diff flux corresponding to Gamma
  const CFreal sigmaf= 1.0;
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT]/sigmaf;

  _flux[6] = (coeffTauMu1 + coeffTauMu2)*(gradGa[XX]*nx + gradGa[YY]*ny);
  
  //diff flux corresponding to Re
  const CFreal sigmatheta= 2.0;
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT];

  _flux[7] = sigmatheta*(coeffTauMu1 + coeffTauMu2)*(gradRe[XX]*nx + gradRe[YY]*ny);
  //jacobian contribution
  if(_isPerturb)
  {
    if(_iPerturbVar != 4){
      _flux[4] = _unperturbedFluxK;
    }
    if(_iPerturbVar != 5){
      _flux[5] = _unperturbedFluxOmega;
    }
    if(_iPerturbVar != 6){
      _flux[6] = _unperturbedFluxGa;
    }
    if(_iPerturbVar != 7){
      _flux[7] = _unperturbedFluxRe;
    }
  }
  else
  {
    _unperturbedFluxK = _flux[4];
    _unperturbedFluxOmega = _flux[5];
    _unperturbedFluxGa    = _flux[6];
    _unperturbedFluxRe    = _flux[7];
  }
 //cout << _flux << endl;
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
