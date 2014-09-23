#include "NavierStokesGReKO2DSourceTerm_Lang.hh"
#include "GReKO/NavierStokes2DGReKOPuvt.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "FiniteVolumeGReKO/FiniteVolumeGReKO.hh"
#include "Framework/SubSystemStatus.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "MathTools/MathConsts.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::GReKO;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

////////////////////////////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NavierStokesGReKO2DSourceTerm_Lang,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeGReKOModule>
NavierStokesGReKO2DSourceTerm_LangFVMCCProvider("NavierStokesGReKO2DSourceTerm_Lang");

///////////////////////////////////////////////////////////////////////////////////////////////////
void NavierStokesGReKO2DSourceTerm_Lang::defineConfigOptions(Config::OptionList& options)
{
 options.addConfigOption< bool >("SSTV","True for SST with Vorticity source term");
 options.addConfigOption< bool >("SSTsust","True for SST with  sustaining terms");
 options.addConfigOption< CFreal >("Kinf","K at the farfield");
 options.addConfigOption< CFreal >("Omegainf","Omega at the farfield");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

NavierStokesGReKO2DSourceTerm_Lang::NavierStokesGReKO2DSourceTerm_Lang(const std::string& name) :
  NavierStokesKOmega2DSourceTerm(name),
  _diffVarSetGReKO(CFNULL),
  _Rethetat(),
  _Rethetac(),
  _Flength(),
  _Vorticity(),
  _Strain(),
  _PGrad(true)

{ 
 addConfigOptionsTo(this);
  _SST_V = false;
  setParameter("SSTV",&_SST_V);
  _SST_sust = false;
  setParameter("SSTsust",&_SST_sust);
  _kamb = 100. ;
  setParameter("Kinf",&_kamb);
  _omegaamb = 0.1;
  setParameter("Omegainf",&_omegaamb);
}

/////////////////////////////////////////////////////////////////////////////////////////////////

NavierStokesGReKO2DSourceTerm_Lang::~NavierStokesGReKO2DSourceTerm_Lang()
{
  for(CFuint iGrad = 0; iGrad < _gradients.size(); iGrad++)
  {
    deletePtr(_gradients[iGrad]);
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////

void NavierStokesGReKO2DSourceTerm_Lang::SetDiffVarset()
{
  _diffVarSetGReKO = getMethodData().getDiffusiveVar().d_castTo<NavierStokes2DGReKOPuvt>();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

CFreal  NavierStokesGReKO2DSourceTerm_Lang::GetStrain(CFreal& VoverRadius)
{
  const CFreal gradU_X = (*(_gradients[1]))[XX];
  const CFreal gradU_Y = (*(_gradients[1]))[YY];
  const CFreal gradV_X = (*(_gradients[2]))[XX];
  const CFreal gradV_Y = (*(_gradients[2]))[YY];
  const CFreal gradSum = (gradU_Y+ gradV_X);
  const CFreal Strain = std::pow(gradU_X,2.)+ 0.5*std::pow(gradSum,2.)+ std::pow(gradV_Y,2.) + VoverRadius ;

  // const CFreal strain1     = 0.5*((*(_gradients[2]))[XX] + (*(_gradients[1]))[YY]);
  // const CFreal strain     = std::sqrt(2*strain1*strain1);

  // _Strain = std::sqrt(2.*gradSum*gradSum);
   _Strain = std::sqrt(2.*Strain);
 return _Strain;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
CFreal  NavierStokesGReKO2DSourceTerm_Lang::GetVorticity()
{
   const CFreal Vorticity1 = ((*(_gradients[2]))[XX] - (*(_gradients[1]))[YY]);
   const CFreal Vorticity2 = 0.5*Vorticity1*Vorticity1;
   _Vorticity =  std::sqrt(2*Vorticity2);
    return _Vorticity;
}
//////////////////////////////////////////////////////////////////////////////////////////////////

void NavierStokesGReKO2DSourceTerm_Lang::computeSource(Framework::GeometricEntity *const element,
						   RealVector& source,
						   RealMatrix& jacobian)
{
  DataHandle<CFreal> normals  = this->socket_normals.getDataHandle();
  DataHandle<CFint> isOutward = this->socket_isOutward.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
 
 
  cf_assert(_varSet.isNotNull());

  // Set the physical data for the cell considered
  State *const currState = element->getState(0);
  _varSet->computePhysicalData(*currState, _physicalData);

  // fill in the nodal states
  const vector<Node*>* const nodes = element->getNodes();
  const CFuint nbNodesInElem = nodes->size();
  _states.clear();
  for (CFuint i = 0; i < nbNodesInElem; ++i) {
    _states.push_back(&_nstates[(*nodes)[i]->getLocalID()]);
  }

  //From now on, we will use the gradient vars
  _diffVarSetGReKO->setGradientVars(_states, _values, _states.size());

  const vector<GeometricEntity*>& faces = *element->getNeighborGeos();
  cf_assert(faces.size() == nbNodesInElem);
  const CFuint elemID = element->getID();

  // compute the gradients by applying Green Gauss in the
  // cell d's
  for(CFuint iGrad = 0; iGrad < _gradients.size(); iGrad++)
  {
    *(_gradients[iGrad]) = 0.0;
  }

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  for (CFuint i = 0; i < nbNodesInElem; ++i) {
    // get the face normal
    const CFuint faceID = faces[i]->getID();
    const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
    CFreal nx = normals[startID];
    CFreal ny = normals[startID + 1];
    if (static_cast<CFuint>( isOutward[faceID]) != elemID) {
      nx *= -1.;
      ny *= -1.;
    }

    if (i < (nbNodesInElem - 1))
    {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	(*(_gradients[iEq]))[XX] += nx*(_values(iEq,i) + _values(iEq,i+1));
	(*(_gradients[iEq]))[YY] += ny*(_values(iEq,i) + _values(iEq,i+1));
      }
    }
    else {
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	(*(_gradients[iEq]))[XX] += nx*(_values(iEq,i) + _values(iEq,0));
	(*(_gradients[iEq]))[YY] += ny*(_values(iEq,i) + _values(iEq,0));
      }
    }
  }




  *(_gradients[0]) *= 0.5/volumes[elemID];
  *(_gradients[1]) *= 0.5/volumes[elemID];
  *(_gradients[2]) *= 0.5/volumes[elemID];
  *(_gradients[3]) *= 0.5/volumes[elemID];
  *(_gradients[4]) *= 0.5/volumes[elemID];
  *(_gradients[5]) *= 0.5/volumes[elemID];
  *(_gradients[6]) *= 0.5/volumes[elemID];
  *(_gradients[7]) *= 0.5/volumes[elemID];

// compute PUVTGReKO by averaging the nodes
// NO!!! If we do it like that we nearly certainly
// get negative values!!!
// So we just take the state value
  const CFuint iK = _varSet->getModel()->getFirstScalarVar(0);
  
  _avState[0] = _physicalData[EulerTerm::P];
  _avState[1] = _physicalData[EulerTerm::VX];
  _avState[2] = _physicalData[EulerTerm::VY];
  _avState[3] = _physicalData[EulerTerm::T];
  _avState[4] = _physicalData[iK];
  _avState[5] = _physicalData[iK+1];
  _avState[6] = _physicalData[iK+2];
  _avState[7] = _physicalData[iK+3];

  double avK     = _physicalData[iK];
  double avOmega = _physicalData[iK+1];
  double avGa    = _physicalData[iK+2];
  double avRe    = _physicalData[iK+3];

  CFreal avV     = _physicalData[EulerTerm::V];
  
  const CFreal rho = _diffVarSetGReKO->getDensity(_avState);

  ///Get the wall distance
  const CFreal avDist = _wallDistance[currState->getLocalID()];

  ///Set the wall distance before computing the turbulent viscosity
  _diffVarSetGReKO->setWallDistance(avDist);

   //Compute Radius for Axisymmetric case and v/Radius
  const CFreal  Radius = (this->getMethodData().isAxisymmetric()) ? (currState->getCoordinates())[YY] : 1.; 
        CFreal  vOverRadius = 0 ; 
  if( (Radius > MathTools::MathConsts::CFrealEps()) && (this->getMethodData().isAxisymmetric())) {
             vOverRadius = _avState[2]/Radius;
 }

  CFreal mut = _diffVarSetGReKO->getTurbDynViscosityFromGradientVars(_avState, _gradients);
   CFreal mu = _diffVarSetGReKO->getLaminarDynViscosityFromGradientVars(_avState);
  _diffVarSetGReKO->computeBlendingCoefFromGradientVars(_avState, *(_gradients[4]), *(_gradients[5]));

  //The MOdified Blending Function F1 
  const CFreal F1org  = _diffVarSetGReKO->getBlendingCoefficientF1();
  const CFreal Ry     = (rho*avDist * sqrt(avK) )/(mu);   
  const CFreal F3bis  = std::pow(Ry/120.0,8);
  const CFreal F3     = exp(-F3bis);
  const CFreal blendingCoefF1     = std::max(F1org,F3);
  
  //const CFreal blendingCoefF1 = _diffVarSetGReKO->getBlendingCoefficientF1();
 
  const CFreal sigmaOmega2 = _diffVarSetGReKO->getSigmaOmega2();


  ///Compute Reynolds stress tensor
  const CFreal twoThirdRhoK = (2./3.)*(avK * rho);
  const CFreal coeffTauMu = _diffVarSetGReKO->getModel().getCoeffTau();
  
  // Get Vorticity
   GetVorticity();

   ///Compute the blending function Fthetat
  const CFreal  Rew         = (rho*avDist * avDist * avOmega)/(mu);   
  const CFreal  Fwake1      = (1e-5 * Rew)*(1e-5 * Rew);   
  const CFreal  Fwake       = exp(-Fwake1);
  const CFreal  thetaBL     = (avRe*mu)/(rho*avV);
  const CFreal  deltaBL     = (0.5*15*thetaBL);
  const CFreal  delta       = (50 * _Vorticity * avDist * deltaBL)/(avV);
  const CFreal  CoefFtheta0 = (avDist/delta)*(avDist/delta)*(avDist/delta)*(avDist/delta);
  const CFreal  CoefFtheta1 = exp(-CoefFtheta0);
  const CFreal  Ftheta1     = Fwake * CoefFtheta1;
  const CFreal  ce2         = 50;
  //const CFreal  overce2     = 1/50;
  //const CFreal  Ftheta2     = (avGa-overce2)/(1.0-overce2);
  const CFreal  Ftheta3     = 1-(((ce2*avGa-1.0)/(ce2-1.0))*((ce2*avGa-1.0)/(ce2-1.0)));
  const CFreal  Ftheta4     = std::max(Ftheta1,Ftheta3);
  const CFreal   Fthetat     = std::min(Ftheta4,1.0);


   
  //The variables needed for the  production term of Re
  const CFreal cthetat   = 0.03;
  const CFreal t         = (500 * mu )/(rho * avV * avV);
   cf_assert(avV >0.);   
      CFreal Tu= 100 * (std::sqrt(2*avK/3))/(avV);


   _Rethetat = (_PGrad) ? GetRethetat(Tu) : GetRethetatwithPressureGradient(mu,Tu); 

   //Compute Flength
    GetFlength(avRe);
   //Compute _Retheta_C
   GetRethetac(avRe);
  //Compute Strain 
    GetStrain(vOverRadius); 
  //ompute Gasep
   const CFreal Rt         = (rho*avK)/(mu*avOmega);
   const CFreal Freattach0 = exp(-Rt/20);
   const CFreal Freattach  = std::pow(Freattach0,4);
   const CFreal Rev        = (rho*avDist*avDist*_Strain)/(mu);
   const CFreal Gasep1     = ((Rev)/((3.235 *  _Rethetac)))-1;
   const CFreal Gasep2     = std::max(0.,Gasep1);
   const CFreal Gasep3     = 2.0*Gasep2*Freattach;
   const CFreal Gasep4     = std::min(Gasep3,2.0);
   const CFreal Gasep      = Gasep4*Fthetat;
   

  ///Gaeff
   const CFreal Gaeff  = std::max(avGa,Gasep);
   

 
  // The Onset function of the  production term of the intermittency Ga

  const CFreal  Fonset1 = (Rev )/(2.93*_Rethetac);
  const CFreal  Fonset2 = std::pow(Fonset1,4);
  const CFreal  Fonset3 = std::max(Fonset1,Fonset2);
  const CFreal  Fonset4 = std::min(Fonset3,2.0);
  const CFreal  Fonset6 = 1-((Rt/2.5)*(Rt/2.5)*(Rt/2.5)) ;
  const CFreal  Fonset7 = std::max(Fonset6,0.);
  const CFreal  Fonset8 = (Fonset4 -Fonset7);
  const CFreal  Fonset  = std::max(Fonset8,0.);

   
  // computeProductionTerm(twoThirdRhoK,mut,vOverRadius,_KProdTerm);
  ///The Modified Production term: k
// const CFreal twoThirdDivV = (2./3.)*((*(_gradients[1]))[XX] + (*(_gradients[2]))[YY]); 
// CFreal tauXX = coeffTauMu * (mut * (2.*(*(_gradients[1]))[XX] - twoThirdDivV)) - twoThirdRhoK;
//  CFreal tauXY = coeffTauMu * (mut * ((*(_gradients[1]))[YY] + (*(_gradients[2]))[XX]));
//  CFreal tauYY = coeffTauMu * (mut * (2.*(*(_gradients[2]))[YY] - twoThirdDivV)) - twoThirdRhoK;
  

  CFreal prodTerm_k = Radius*(Gaeff*_KProdTerm);
 // CFreal prodTerm_k1 = Gaeff*(tauXX*(*(_gradients[1]))[XX] + tauYY*(*(_gradients[2]))[YY] + tauXY * ((*(_gradients[1]))[YY] + (*(_gradients[2]))[XX]));

  ///Production term: Omega
  CFreal prodTerm_Omega = (_diffVarSetGReKO->getGammaCoef()*rho/mut)*_KProdTerm;

  ///This is used in (BSL,SST), not for normal kOmega
  const CFreal overOmega = 1./avOmega;
  prodTerm_Omega += (1. - blendingCoefF1) * 2. * rho * overOmega * sigmaOmega2
                    * ((*(_gradients[4]))[XX]*(*(_gradients[5]))[XX] + (*(_gradients[4]))[YY]*(*(_gradients[5]))[YY]);

  prodTerm_Omega *= Radius; 


  ///The Modified Destruction term: k
  ///CoeffDk This coefficient is used in the destruciton term related to k: 
    const CFreal CoeffDk1  = std::max(Gaeff,0.1);
    const CFreal CoeffDk   = std::min(CoeffDk1,1.0);
    CFreal destructionTerm_k = Radius*CoeffDk*((-1.) * rho * avOmega * avK * _diffVarSetGReKO->getBetaStar(_avState));
  // CFreal destructionTerm_k = ((-1.) * rho * avOmega * avK * _diffVarSetGReKO->getBetaStar(_avState));

  ///Destruction term: Omega
  CFreal destructionTerm_Omega = Radius* (-1.) * rho * avOmega * avOmega * _diffVarSetGReKO->getBeta(_avState);

  //Make sure negative values dont propagate...
  destructionTerm_k     = min(0., destructionTerm_k);
  destructionTerm_Omega = min(0., destructionTerm_Omega);
  prodTerm_k            = max(0., prodTerm_k);
  prodTerm_Omega        = max(0., prodTerm_Omega);



  // This trick was used by W. Dieudonne in Euphoria (with 20.)
  prodTerm_k = min(10.*fabs(destructionTerm_k), prodTerm_k);
  prodTerm_Omega = min(10.*fabs(destructionTerm_Omega), prodTerm_Omega);

   // The production term of the intermittency Ga
  const CFreal ca1       = 2.0;
  const CFreal ce1       = 1.0;
  const CFreal  ca2 =  0.06;
  const CFreal GaFonset1 = avGa * Fonset;
  const CFreal GaFonset  = std::pow(GaFonset1,0.5);

    CFreal prodTerm_Ga = _Flength * ca1 * rho * _Strain * GaFonset * (1.0 - ce1*avGa);
      
       //prodTerm_Ga *=Radius;
  // CFreal prodTerm_Ga = _Flength * ca1 * rho * strain * GaFonset * 1.0 +  ca2 * rho *  vorticity * avGa * Fturb * (ce2*avGa);

  // The production term of  Re
    CFreal prodTerm_Re = cthetat * (rho/t) * (_Rethetat - avRe) * (1.0 - Fthetat);
     //prodTerm_Re *=Radius; 
   // CFreal prodTerm_Re = cthetat * (rho/t) * (_Rethetat) * (1.0 - Fthetat);

  //The variables needed for the  Destruction term of Ga   
  const CFreal  Fturb1 =  exp(-Rt/4); 
  const CFreal  Fturb =  std::pow(Fturb1,4); 
  //Destruction term of the intermittency Ga
      CFreal  destructionTerm_Ga  = (-1.0) *ca2 * rho *  _Vorticity * avGa * Fturb * (ce2*avGa - 1);
  
         destructionTerm_Ga *=Radius;
  //Destruction term of the intermittency Re
     CFreal destructionTerm_Re = 0;



  ///Make sure negative values dont propagate
   prodTerm_Ga        = max(0., prodTerm_Ga);
   prodTerm_Re        = max(0., prodTerm_Re);
   destructionTerm_Ga = min(0., destructionTerm_Ga);
   destructionTerm_Re = min(0., prodTerm_Re);

  
  CFreal Tau_tt = (-2./3.)*coeffTauMu*((*(_gradients[1]))[XX] + (*(_gradients[2]))[YY] - 2*vOverRadius);

  //Computation of the source term
  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = (this->getMethodData().isAxisymmetric()) ? _avState[0]-Tau_tt : 0.0 ;


  //What we do with the source term depends if
  //we are computing the jacobian or not
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();
  if(isPerturb)
  {
    /// Compute the jacobian contribution
    // only perturb the negative part of the source term
    if(iPerturbVar == 4)
    {
      source[4] = destructionTerm_k;
      source[4] += _unperturbedPositivePart[0];
    }
    else
    {
      source[4] = _unperturbedNegativePart[0];
      source[4] += _unperturbedPositivePart[0];
    }

    if(iPerturbVar == 5)
    {
      source[5] = destructionTerm_Omega;
      source[5] += _unperturbedPositivePart[1];
    }
    else
    {
      source[5] = _unperturbedNegativePart[1];
      source[5] += _unperturbedPositivePart[1];
    }


    if(iPerturbVar == 6)
    {
      source[6] = destructionTerm_Ga;
      source[6] += _unperturbedPositivePart[2];
    }
    else
    {
      source[6] = _unperturbedNegativePart[2];
      source[6] += _unperturbedPositivePart[2];
    }
    if(iPerturbVar == 7)
    {
      source[7] = destructionTerm_Re;
      source[7] += _unperturbedPositivePart[3];
    }
    else
    {
      source[7] = _unperturbedNegativePart[3];
      source[7] += _unperturbedPositivePart[3];
    }


  }
  else
  {
    /// Compute the rhs contribution
    // and Store the unperturbed source terms
    source[4] = prodTerm_k;
    source[4] += destructionTerm_k;
    _unperturbedPositivePart[0] = prodTerm_k;
    _unperturbedNegativePart[0] = destructionTerm_k;

    source[5] = prodTerm_Omega;
    source[5] += destructionTerm_Omega;
    _unperturbedPositivePart[1] = prodTerm_Omega;
    _unperturbedNegativePart[1] = destructionTerm_Omega;
    
    source[6] = prodTerm_Ga;
    source[6] += destructionTerm_Ga;
    _unperturbedPositivePart[2] = prodTerm_Ga;
    _unperturbedNegativePart[2] = destructionTerm_Ga;
    
    source[7] = prodTerm_Re;
    source[7] += destructionTerm_Re;
    _unperturbedPositivePart[3] = prodTerm_Re;
    _unperturbedNegativePart[3] = destructionTerm_Re;
  }

  ///Finally multiply by the cell volume
  source *= volumes[elemID];
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CFreal NavierStokesGReKO2DSourceTerm_Lang::GetRethetac(CFreal& Retheta)
{	
        if (Retheta <= 1860){
            _Rethetac  = Retheta - (396.035*1e-2 -120.656*1e-4*Retheta)+(868.230*1e-6)*Retheta*Retheta 
                          - 696.506*1e-9*Retheta*Retheta*Retheta + 174.105*1e-12*Retheta*Retheta*Retheta*Retheta;
           }
        else {
            _Rethetac = Retheta - 593.11 + (Retheta - 1870.0)*0.482;
              }
  return _Rethetac;
}           

/////////////////////////////////////////////////////////////////////////////////////////////////////
CFreal NavierStokesGReKO2DSourceTerm_Lang::GetFlength(CFreal& Retheta)
{
       if (Retheta <=400){  
            _Flength = 398.189*1e-1 -119.270*1e-4*Retheta -132.567*1e-6*Retheta*Retheta;   
          }
      else if ((Retheta>=400 ) && (Retheta < 596)) {
            _Flength = 263.404 - 123.939*1e-2*Retheta + 194.548*1.e-5*Retheta*Retheta - 101.695*1e-8*Retheta*Retheta*Retheta;
              }
        else if ((Retheta>=596 ) && (Retheta < 1200)) {
            _Flength = 0.5-(Retheta - 596.0)*3.0*1e-4;
              }
        else {
            _Flength = 0.3188;
          }
return _Flength; 

}
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CFreal NavierStokesGReKO2DSourceTerm_Lang::GetRethetat(const CFreal& Tu) 
{
   cf_assert(Tu > 0);   
  const CFreal overTu    = 1/Tu;
           if (Tu<=1.3) {
                _Rethetat = (1173.51-589.428*Tu + 0.2196*overTu*overTu);
		 }
    	  else {
  		const CFreal lamco5   = Tu - 0.5658;
  	 	const CFreal pwtu   = -0.671;
  	 	_Rethetat = 331.5*std::pow(lamco5,pwtu);
         	 }
             return _Rethetat;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CFreal NavierStokesGReKO2DSourceTerm_Lang::GetRethetatwithPressureGradient(CFreal& Viscosity, const CFreal& Tu)
{
        const CFreal overTu    = 1/Tu;
        const CFreal avV = _physicalData[EulerTerm::V];   //AvrageSpeeed;   
        const CFreal mu = Viscosity;
        const CFreal rho = _diffVarSetGReKO->getDensity(_avState); ;  
        const CFreal avu     = _physicalData[EulerTerm::VX];
        const CFreal avv     = _physicalData[EulerTerm::VY];        
        const CFreal VSound  = _physicalData[EulerTerm::A];
  	vector<CFreal> Theta(2);
        
         CFreal  Flambdak = 0;
         CFreal  Flambdakprime = 0;
         if(Tu <=1.3){
    	Theta[0]            = (mu/(rho*avV))*(1173.51-589.428*Tu + 0.2196*overTu*overTu);
 	//cout << "thetaN" << Theta[0] <<endl; 
  	}
  	else {
   	const CFreal lamco5   = Tu - 0.5658;
   	const CFreal pwtu     = -0.671;
    	Theta[0]            = (mu/(rho*avV))*331.5*std::pow(lamco5,pwtu);
 	//cout << "thetaP" << Theta[0] << endl; 
   	}
  	const CFuint MAXITER   = 10;
  	const CFreal TOL       = 1e-6;
  	CFreal Restheta        =  Theta[0]*TOL;
 	for (CFuint i = 0; i < MAXITER; ++i)
 	{
	  //cout << "ITER" << i  << endl;
 	 //The variables needed for the calculayion of Re_thetat
  	const CFreal overU     = 1./avV;
  	const CFreal overmu    = 1./mu;
  	const CFreal MACH      = avV/VSound;
  	const CFreal Fc1       = 1+16*(2*MACH)*(2*MACH)*(2*MACH)*(2*MACH);
  	const CFreal dUdx   	 = avV * (avu* (*(_gradients[1]))[XX]  + avv * (*(_gradients[2]))[XX]);
  	const CFreal dUdy      = avV * (avu* (*(_gradients[1]))[YY]  + avv * (*(_gradients[2]))[YY]);
  	const CFreal dUds      =  overU * (avu * dUdx + avv *dUdy);
  	const CFreal Fc        = std::pow(Fc1,0.75);
   	CFreal lambda1  = rho * Theta[0]* Theta[0] * overmu * dUds*Fc;
   	CFreal lambda2 = std::max(lambda1,-0.1);
   	CFreal lambda =  std::min(lambda2,0.1);
 	 const CFreal lambdaprime = 2 * rho * Theta[0] * overmu * dUds;
  	 if (lambda <=0) {
   	const CFreal lamco1   = (Tu/1.5);
   	const CFreal lamco2   = -1.0*std::pow(lamco1,1.5);
   	const CFreal FlambN   =  -12.986 * lambda - 123.66 * lambda*lambda -405.689 * lambda*lambda*lambda;
   	const CFreal FlambprimeN   = -12.986 - 2 * 123.66 * lambda*lambdaprime - 3* 405.689 * lambda*lambda*lambdaprime;
   	 Flambdak      = 1 - FlambN * std::exp(lamco2);
    	Flambdakprime =  -1.0*FlambprimeN * std::exp(lamco2);
 	 }
   	else {
   	const CFreal lamco3   = -1.0*(Tu/0.5);
   	const CFreal lamco4   = -35.0*lambda;
   	const CFreal FlambP   = 0.275*(1-std::exp(lamco4))*std::exp(lamco3);
   	const CFreal FlambprimeP   = 0.275*(1-(-35.0*lambdaprime)*std::exp(lamco4))*std::exp(lamco3); 
    	Flambdak= 1 + FlambP;
    	Flambdakprime= FlambprimeP;
   	}
  	if (Tu<=1.3) {
    	 CFreal Rethetat0 = (1173.51-589.428*Tu + 0.2196*overTu*overTu);
   	// const CFreal Rethetatprime = Rethetat0 *Flambdakprime;
         _Rethetat = Rethetat0 * Flambdak;
  	}
   	else {
   	const CFreal lamco5   = Tu - 0.5658;
   	const CFreal pwtu   = -0.671;
   	const CFreal Rethetat0 = 331.5*std::pow(lamco5,pwtu);
  	// const CFreal Rethetatprime = Rethetat0 * Flambdakprime;
         _Rethetat = Rethetat0 * Flambdak;
 	 }

   	const CFreal Rethetatprime =  (_Rethetat* Flambdakprime)/Flambdak;

     
   	const CFreal MainF = _Rethetat -(rho*avV)*Theta[0]/mu;
   	const CFreal MainFprime = Rethetatprime -(rho*avV)/mu;
    
    	Theta[1] = Theta[0] - MainF/MainFprime;
  	if ( std::abs(Theta[0]-Theta[1]) <= Restheta ) break;
    	Theta[0]= Theta[1];
   	//cout << "diff   " << Theta[0]-Theta[1] << endl;
   	// if ( i == MAXITER -1 ) cout << "GREK0: Newton convergence was  reached! at    " << i<< endl;
    	//if ( i == MAXITER - 1) cout << "GREK0: Newton convergence was  reached!" << endl;
  	}
 
 
  	//Compute Flength
 	//const CFreal FLENGTH1  = 0.1*exp(-0.022*avRe+12)+0.45;
      	// _Flength  = 8.5e7*std::pow(_Rethetat,-3)*(1+0.25*MACH*MACH*MACH*MACH);
//	}

return _Rethetat;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
