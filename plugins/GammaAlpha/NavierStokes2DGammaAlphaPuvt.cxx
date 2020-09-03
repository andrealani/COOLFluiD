#include "GammaAlpha.hh"
#include "NavierStokes2DGammaAlphaPuvt.hh"
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

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes2DGammaAlphaPuvt, DiffusiveVarSet,
			    GammaAlphaModule, 2>
ns2DGammaAlphaPuvtProvider("NavierStokes2DGammaAlphaPuvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGammaAlphaPuvt::NavierStokes2DGammaAlphaPuvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokesKLogOmegaPvt<NavierStokesKLogOmegaVarSet<NavierStokes2DVarSet, 0> >(name, model),
  _unperturbedFluxGa(),
  _unperturbedFluxAlpha()
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
  names[6] = "gamma";
  names[7] = "alpha";
  
  setVarNames(names);
  setModelCoefficients();
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGammaAlphaPuvt::~NavierStokes2DGammaAlphaPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& NavierStokes2DGammaAlphaPuvt::getFlux(const RealVector& values,
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
  const RealVector& gradAlpha    = *gradients[7];

  const CFreal p = max(values[0],1.0e-5);
  const CFreal avU = values[1];
  const CFreal avV = values[2];
  const CFreal gamma = min(max(values[6],0.01),0.99);
  //const CFreal gamma = min(max(values[6],0.0),1.0);

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
  
  // compute local freestream values with isentropic assumption
  const CFreal MInf = _eulerModel->getMachInf();
  const CFreal pInf = _eulerModel->getStaticPressInf();//getPressInf();
  const CFreal gammaIsentropic = _eulerModel->getGamma();
  const CFreal uInf = _eulerModel->getVelInf();
  const CFreal rhoInf = (MInf*MInf)/(uInf*uInf)*gammaIsentropic*pInf;
  
  const CFreal pTotTerm = 1.0+(gammaIsentropic-1.0)/2.0*MInf*MInf;
  const CFreal pTotExponent = gammaIsentropic/(gammaIsentropic-1.0);
  const CFreal pTotalInf = pow(pTotTerm,pTotExponent)*p;
    
  const CFreal R = _eulerModel->getR();
  
  const CFreal rhoInfLocal = rhoInf*pow(p/pInf,1.0/gammaIsentropic);
  const CFreal MInfLocal = sqrt(2.0/(gammaIsentropic-1.0)*(pow(pTotalInf/p,1.0/pTotExponent)-1.0));
  const CFreal uInfLocal = sqrt(gammaIsentropic*p/rhoInfLocal)*MInfLocal;
  const CFreal TInfLocal = p/(rhoInfLocal*R);
  
  RealVector copyState = values;
  copyState[3] = TInfLocal;
  
  const CFreal muInfLocal = getLaminarDynViscosityFromGradientVars(copyState);
  
  //const CFreal distance = std::max(_wallDistance, 1.0e-12);
  
  const CFreal fMuGamma = 1.0-exp(-256.0*(_wallDistance*uInfLocal*rhoInfLocal/muInfLocal)*(_wallDistance*uInfLocal*rhoInfLocal/muInfLocal));
  const CFreal fMMuGamma = (1.0+0.26*(gammaIsentropic-1.0)/2.0*MInfLocal*MInfLocal)*sqrt(1+0.38*pow(MInfLocal,0.6));
 //if (_wallDistance < 1.0e-4) CFLog(INFO, "wallDist: " << _wallDistance << "\n"); 
  //const CFreal gammaLim = min(max(0.01,gamma),0.99);
  const CFreal muGamma = 0.57*pow(-log(1-gamma),-5.0/6.0*(1.0-gamma))*fMuGamma*fMMuGamma*getModel().getCoeffTau()*nsData[NSTurbTerm::MU];

  _flux[6] = (muGamma)*(gradGa[XX]*nx + gradGa[YY]*ny);
  
  //diff flux corresponding to alpha
  const CFreal sigmatheta= 2.0;
  coeffTauMu1 = getModel().getCoeffTau()*nsData[NSTurbTerm::MU];
  coeffTauMu2 = getModel().getCoeffTau()*nsData[NSTurbTerm::MUT];

  _flux[7] = sigmatheta*(coeffTauMu1 + coeffTauMu2)*(gradAlpha[XX]*nx + gradAlpha[YY]*ny);
  
  //if (values[6]<0) CFLog(INFO, "Flux update: " << _flux << "\n");
  
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
      _flux[7] = _unperturbedFluxAlpha;
    }
  }
  else
  {
    _unperturbedFluxK = _flux[4];
    _unperturbedFluxOmega = _flux[5];
    _unperturbedFluxGa    = _flux[6];
    _unperturbedFluxAlpha    = _flux[7];
  }
 //cout << _flux << endl;
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
