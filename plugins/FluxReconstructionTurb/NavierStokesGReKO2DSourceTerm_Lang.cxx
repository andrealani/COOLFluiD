#include "NavierStokesGReKO2DSourceTerm_Lang.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "FluxReconstructionTurb/FluxReconstructionTurb.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluxReconstructionTurb/KOmega2DSourceTerm.hh"


#include "Framework/MethodStrategyProvider.hh"
#include "MathTools/MathConsts.hh"
#include "KOmega/NavierStokesKOmegaVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::KOmega;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

////////////////////////////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NavierStokesGReKO2DSourceTerm_Lang, FluxReconstructionSolverData, FluxReconstructionKOmegaModule>
NavierStokesGReKO2DSourceTerm_LangFRProvider("NavierStokesGReKO2DSourceTerm_Lang");

///////////////////////////////////////////////////////////////////////////////////////////////////

void NavierStokesGReKO2DSourceTerm_Lang::defineConfigOptions(Config::OptionList& options)
{
 options.addConfigOption< bool >("SSTV","True for SST with Vorticity source term");
 options.addConfigOption< bool >("SSTsust","True for SST with  sustaining terms");
 options.addConfigOption< CFreal >("Kinf","K at the farfield");
 options.addConfigOption< CFreal >("Omegainf","Omega at the farfield");
 options.addConfigOption< bool >("PGrad","pressure Gradient");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

NavierStokesGReKO2DSourceTerm_Lang::NavierStokesGReKO2DSourceTerm_Lang(const std::string& name) :
  KOmega2DSourceTerm(name),
  m_avState(),
  m_Rethetat(),
  m_Rethetac(),
  m_Flength(),
  m_Vorticity(),
  m_Strain()

{ 
 addConfigOptionsTo(this);
  m_SST_V = false;
  setParameter("SSTV",&_SST_V);
  m_SST_sust = false;
  setParameter("SSTsust",&_SST_sust);
  m_kamb = 100. ;
  setParameter("Kinf",&_kamb);
  m_omegaamb = 0.1;
  setParameter("Omegainf",&_omegaamb);
  m_PGrad = false;
  setParameter("PGrad",&_PGrad);
}

/////////////////////////////////////////////////////////////////////////////////////////////////

NavierStokesGReKO2DSourceTerm_Lang::~NavierStokesGReKO2DSourceTerm_Lang()
{
  for(CFuint iGrad = 0; iGrad < m_gradients.size(); iGrad++)
  {
    deletePtr(m_gradients[iGrad]);
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////

void NavierStokesGReKO2DSourceTerm_Lang::setup()
{
  KOmega2DSourceTerm::setup();
  
  m_avState.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////////////////////////

CFreal  NavierStokesGReKO2DSourceTerm_Lang::GetStrain(CFreal& VoverRadius)
{
  const CFreal gradU_X = (*(m_gradients[1]))[XX];
  const CFreal gradU_Y = (*(m_gradients[1]))[YY];
  const CFreal gradV_X = (*(m_gradients[2]))[XX];
  const CFreal gradV_Y = (*(m_gradients[2]))[YY];
  const CFreal gradSum = (gradU_Y+ gradV_X);
  const CFreal Strain = std::pow(gradU_X,2.)+ 0.5*std::pow(gradSum,2.)+ std::pow(gradV_Y,2.) + VoverRadius ;
   m_Strain = std::sqrt(2.*Strain);
 return m_Strain;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

CFreal  NavierStokesGReKO2DSourceTerm_Lang::GetVorticity()
{
   const CFreal Vorticity1 = ((*(m_gradients[2]))[XX] - (*(m_gradients[1]))[YY]);
   const CFreal Vorticity2 = 0.5*Vorticity1*Vorticity1;
   m_Vorticity =  std::sqrt(2*Vorticity2);
    return m_Vorticity;
}
//////////////////////////////////////////////////////////////////////////////////////////////////

void NavierStokesGReKO2DSourceTerm_Lang::computeSource(Framework::GeometricEntity *const element,
						   RealVector& source,
						   RealMatrix& jacobian)
{
    PreparecomputeSource(element);
// compute PUVTGReKO by averaging the nodes
// NO!!! If we do it like that we nearly certainly
// get negative values!!!
// So we just take the state value
  const CFuint iK = m_varSet->getModel()->getFirstScalarVar(0);
  m_avState[0] = m_physicalData[EulerTerm::P];
  m_avState[1] = m_physicalData[EulerTerm::VX];
  m_avState[2] = m_physicalData[EulerTerm::VY];
  m_avState[3] = m_physicalData[EulerTerm::T];
  m_avState[4] = m_physicalData[iK];
  m_avState[5] = m_physicalData[iK+1];
  m_avState[6] = m_physicalData[iK+2];
  m_avState[7] = m_physicalData[iK+3];

  CFreal avV     = m_physicalData[EulerTerm::V];
  double avK     = m_physicalData[iK];
  double avOmega = m_physicalData[iK+1];
  double avGa    = m_physicalData[iK+2];
  double avRe    = m_physicalData[iK+3];
  const CFreal rho = m_diffVarSet->getDensity(m_avState);

  ///Get the wall distance
  m_diffVarSet->setWallDistance(m_avDist);
   CFreal mut = m_diffVarSet->getTurbDynViscosityFromGradientVars(m_avState, m_gradients);
   CFreal mu = m_diffVarSet->getLaminarDynViscosityFromGradientVars(m_avState);
  m_diffVarSet->computeBlendingCoefFromGradientVars(m_avState, *(m_gradients[4]), *(m_gradients[5]));

  //The MOdified Blending Function F1 
  //const CFreal F1org  = _diffVarSet->getBlendingCoefficientF1();
  //const CFreal Ry     = (rho*_avDist * sqrt(avK) )/(mu);   
  //const CFreal F3bis  = std::pow(Ry/120.0,8);
  //const CFreal F3     = exp(-F3bis);
  //const CFreal blendingCoefF1     = std::max(F1org,F3);
  
  //const CFreal blendingCoefF1 = _diffVarSet->getBlendingCoefficientF1();
 
  //const CFreal sigmaOmega2 = _diffVarSet->getSigmaOmega2();


  ///Compute Reynolds stress tensor
  //const CFreal twoThirdRhoK = (2./3.)*(avK * rho);
  //const CFreal coeffTauMu = _diffVarSet->getModel().getCoeffTau();
  
  // Get Vorticity
   GetVorticity();

   ///Compute the blending function Fthetat
  const CFreal  Rew         = (rho*m_avDist * m_avDist * avOmega)/(mu);   
  const CFreal  Fwake1      = (1e-5 * Rew)*(1e-5 * Rew);   
  const CFreal  Fwake       = exp(-Fwake1);
  const CFreal  thetaBL     = (avRe*mu)/(rho*avV);
  const CFreal  deltaBL     = (0.5*15*thetaBL);
  const CFreal  delta       = (50 * m_Vorticity * m_avDist * deltaBL)/(avV);
  const CFreal  CoefFtheta0 = (m_avDist/delta)*(m_avDist/delta)*(m_avDist/delta)*(m_avDist/delta);
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


   m_Rethetat = (!m_PGrad) ? GetRethetat(Tu) : GetRethetatwithPressureGradient(mu,Tu); 

   //Compute Flength
    GetFlength(avRe);
   //Compute _Retheta_C
   GetRethetac(avRe);
  //Compute Strain 
    GetStrain(_vOverRadius); 
  //ompute Gasep
   const CFreal Rt         = (rho*avK)/(mu*avOmega);
   const CFreal Freattach0 = exp(-Rt/20);
   const CFreal Freattach  = std::pow(Freattach0,4);
   const CFreal Rev        = (rho*m_avDist*m_avDist*m_Strain)/(mu);
   const CFreal Gasep1     = ((Rev)/((3.235 *  m_Rethetac)))-1;
   const CFreal Gasep2     = std::max(0.,Gasep1);
   const CFreal Gasep3     = 2.0*Gasep2*Freattach;
   const CFreal Gasep4     = std::min(Gasep3,2.0);
   const CFreal Gasep      = Gasep4*Fthetat;
   

  ///Gaeff
   const CFreal Gaeff  = std::max(avGa,Gasep);
   
 
  // The Onset function of the  production term of the intermittency Ga

  const CFreal  Fonset1 = (Rev )/(2.93*m_Rethetac);
  const CFreal  Fonset2 = std::pow(Fonset1,4);
  const CFreal  Fonset3 = std::max(Fonset1,Fonset2);
  const CFreal  Fonset4 = std::min_(Fonset3,2.0);
  const CFreal  Fonset6 = 1-((Rt/2.5)*(Rt/2.5)*(Rt/2.5)) ;
  const CFreal  Fonset7 = std::max(Fonset6,0.);
  const CFreal  Fonset8 = (Fonset4 -Fonset7);
  const CFreal  Fonset  = std::max(Fonset8,0.);

   

  ///The Modified Production  term: k
    ///Gaeff This coefficient is used in the destruciton term related to k: 
    computeProductionTerm(_avState, Gaeff,mut, m_prodTerm_k,m_prodTerm_Omega);
    
  ///The Modified Destruction term: k
  ///CoeffDk This coefficient is used in the destruciton term related to k: 
    const CFreal CoeffDk1  = std::max(Gaeff,0.1);
    const CFreal CoeffDk   = std::min(CoeffDk1,1.0);
    computeDestructionTerm(_avState, CoeffDk,m_destructionTerm_k, m_destructionTerm_Omega);
     
    //Limit the production terms
    LimitProductionTerms();

   // The production term of the intermittency Ga
  const CFreal ca1       = 2.0;
  const CFreal ce1       = 1.0;
  const CFreal  ca2 =  0.06;
  const CFreal GaFonset1 = avGa * Fonset;
  const CFreal GaFonset  = std::pow(GaFonset1,0.5);

    CFreal prodTerm_Ga = m_Flength * ca1 * rho * m_Strain * GaFonset * (1.0 - ce1*avGa);
  
// The production term of  Re
    CFreal prodTerm_Re = cthetat * (rho/t) * (m_Rethetat - avRe) * (1.0 - Fthetat);

  //The variables needed for the  Destruction term of Ga   
  const CFreal  Fturb1 =  exp(-Rt/4); 
  const CFreal  Fturb =  std::pow(Fturb1,4); 
  //Destruction term of the intermittency Ga
      CFreal  destructionTerm_Ga  = (-1.0) *ca2 * rho *  m_Vorticity * avGa * Fturb * (ce2*avGa - 1);
  
         destructionTerm_Ga *= m_Radius;
  //Destruction term of the intermittency Re
     CFreal destructionTerm_Re = 0;

  ///Make sure negative values dont propagate
   prodTerm_Ga        = max(0., prodTerm_Ga);
   prodTerm_Re        = max(0., prodTerm_Re);
   destructionTerm_Ga = min(0., destructionTerm_Ga);
   destructionTerm_Re = min(0., prodTerm_Re);

  

  //Computation of the source term
  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = (m_isAxisymmetric) ? GetNSSourceTerm() : 0.0 ;
  source[3] = 0.0;


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
      source[4] = m_destructionTerm_k;
      source[4] += m_unperturbedPositivePart[0];
    }
    else
    {
      source[4] = m_unperturbedNegativePart[0];
      source[4] += m_unperturbedPositivePart[0];
    }

    if(iPerturbVar == 5)
    {
      source[5] = m_destructionTerm_Omega;
      source[5] += m_unperturbedPositivePart[1];
    }
    else
    {
      source[5] = m_unperturbedNegativePart[1];
      source[5] += m_unperturbedPositivePart[1];
    }


    if(iPerturbVar == 6)
    {
      source[6] = destructionTerm_Ga;
      source[6] += m_unperturbedPositivePart[2];
    }
    else
    {
      source[6] = m_unperturbedNegativePart[2];
      source[6] += m_unperturbedPositivePart[2];
    }
    if(iPerturbVar == 7)
    {
      source[7] = destructionTerm_Re;
      source[7] += m_unperturbedPositivePart[3];
    }
    else
    {
      source[7] = m_unperturbedNegativePart[3];
      source[7] += m_unperturbedPositivePart[3];
    }


  }
  else
  {
    /// Compute the rhs contribution
    // and Store the unperturbed source terms
    source[4] = m_prodTerm_k;
    source[4] += m_destructionTerm_k;
    m_unperturbedPositivePart[0] = m_prodTerm_k;
    m_unperturbedNegativePart[0] = m_destructionTerm_k;

    source[5] = m_prodTerm_Omega;
    source[5] += m_destructionTerm_Omega;
    m_unperturbedPositivePart[1] = m_prodTerm_Omega;
    m_unperturbedNegativePart[1] = m_destructionTerm_Omega;
    
    source[6] = prodTerm_Ga;
    source[6] += destructionTerm_Ga;
    m_unperturbedPositivePart[2] = prodTerm_Ga;
    m_unperturbedNegativePart[2] = destructionTerm_Ga;
    
    source[7] = prodTerm_Re;
    source[7] += destructionTerm_Re;
    m_unperturbedPositivePart[3] = prodTerm_Re;
    m_unperturbedNegativePart[3] = destructionTerm_Re;
  }

  ///Finally multiply by the cell volume
  source *= m_volumes_elemID;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CFreal NavierStokesGReKO2DSourceTerm_Lang::GetRethetac(CFreal& Retheta)
{	
        if (Retheta <= 1860){
            m_Rethetac  = Retheta - (396.035*1e-2 -120.656*1e-4*Retheta)+(868.230*1e-6)*Retheta*Retheta 
                          - 696.506*1e-9*Retheta*Retheta*Retheta + 174.105*1e-12*Retheta*Retheta*Retheta*Retheta;
           }
        else {
            m_Rethetac = Retheta - 593.11 + (Retheta - 1870.0)*0.482;
              }
  return m_Rethetac;
}           

/////////////////////////////////////////////////////////////////////////////////////////////////////
CFreal NavierStokesGReKO2DSourceTerm_Lang::GetFlength(CFreal& Retheta)
{
       if (Retheta <=400){  
            m_Flength = 398.189*1e-1 -119.270*1e-4*Retheta -132.567*1e-6*Retheta*Retheta;   
          }
      else if ((Retheta>=400 ) && (Retheta < 596)) {
            m_Flength = 263.404 - 123.939*1e-2*Retheta + 194.548*1.e-5*Retheta*Retheta - 101.695*1e-8*Retheta*Retheta*Retheta;
              }
        else if ((Retheta>=596 ) && (Retheta < 1200)) {
            m_Flength = 0.5-(Retheta - 596.0)*3.0*1e-4;
              }
        else {
            m_Flength = 0.3188;
          }
return m_Flength; 

}
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CFreal NavierStokesGReKO2DSourceTerm_Lang::GetRethetat(const CFreal& Tu) 
{
   cf_assert(Tu > 0);   
  const CFreal overTu    = 1/Tu;
           if (Tu<=1.3) {
                m_Rethetat = (1173.51-589.428*Tu + 0.2196*overTu*overTu);
		 }
    	  else {
  		const CFreal lamco5   = Tu - 0.5658;
  	 	const CFreal pwtu   = -0.671;
  	 	m_Rethetat = 331.5*std::pow(lamco5,pwtu);
         	 }
             return m_Rethetat;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CFreal NavierStokesGReKO2DSourceTerm_Lang::GetLambda(CFreal& Lambda, CFreal& Theta, CFreal& Viscosity)
{
        const CFreal avV = m_physicalData[EulerTerm::V];   //AvrageSpeeed;   
        const CFreal mu = Viscosity;
        const CFreal rho = m_diffVarSet->getDensity(m_avState); 
        const CFreal rhoovermu = rho /mu;   
        const CFreal avu     = m_physicalData[EulerTerm::VX];
        const CFreal avv     = m_physicalData[EulerTerm::VY];        
  	const CFreal overU     = 1./avV;
  	const CFreal dUdx   	 = avV * (avu* (*(m_gradients[1]))[XX]  + avv * (*(m_gradients[2]))[XX]);
  	const CFreal dUdy      = avV * (avu* (*(m_gradients[1]))[YY]  + avv * (*(m_gradients[2]))[YY]);
  	const CFreal dUds      =  overU * (avu * dUdx + avv *dUdy);
        const CFreal theta_sq  = Theta * Theta;
        const CFreal lambda0 =  rhoovermu * theta_sq * dUds;

   	CFreal lambda1 = std::max(lambda0,-0.1);
   	   Lambda =  std::min(lambda1,0.1);

 return Lambda;
}        
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CFreal NavierStokesGReKO2DSourceTerm_Lang::GetFlambda(CFreal& Lambda,const CFreal& Tu,CFreal& Flambda,CFreal& Theta,bool Prime )
{
       const  CFreal lambdaprime         = Lambda*Theta;
  if (Lambda > 0) {
   	const CFreal lamco1        = (Tu/1.5);
   	const CFreal lamco2        = -1.0*std::pow(lamco1,1.5);
   	const CFreal Flamb        =  -12.986 * Lambda - 123.66 * Lambda*Lambda - 405.689 * Lambda*Lambda*Lambda;
   	const CFreal Flambprime   =  -12.986 * lambdaprime  - 2 * 123.66 * Lambda*lambdaprime - 3 * 405.689 * Lambda*Lambda*lambdaprime;
   	             Flambda      = (Prime)? 1 - (Flamb * std::exp(lamco2)): -1.0*Flambprime * std::exp(lamco2);
  }
 else {
       	const CFreal lamco3   = -1.0*(Tu/0.5);
   	const CFreal lamco4   = -35.0*Lambda;
   	const CFreal FlambP   = 0.275*(1-std::exp(lamco4))*std::exp(lamco3);
   	const CFreal FlambprimeP   = 0.275*(1-((-35.0*lambdaprime)*std::exp(lamco4)))*std::exp(lamco3); 
    	            Flambda        = (Prime)? 1 + FlambP : FlambprimeP;
  }

return Flambda;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CFreal NavierStokesGReKO2DSourceTerm_Lang::GetRethetatwithPressureGradient(CFreal& Viscosity, const CFreal& Tu)
{

        const CFreal avV = m_physicalData[EulerTerm::V];   //AvrageSpeeed;     
        const CFreal mu = Viscosity;                                          
        const CFreal rho = m_diffVarSet->getDensity(m_avState);                 
  	const CFreal overTu   = 1./Tu; 
        vector<CFreal> Theta(2);
        CFreal Lambda = 0;
        CFreal Flambda = 0;
        CFreal FlambdaPrime = 0;
         if(Tu <=1.3){
    	     Theta[0] = (mu/(rho*avV))*(1173.51-589.428*Tu + 0.2196*overTu*overTu);
  	   }
  	else {
   	     const CFreal lamco5   = Tu - 0.5658;
   	     const CFreal pwtu     = -0.671;
    	     Theta[0]            = (mu/(rho*avV))*331.5*std::pow(lamco5,pwtu);
   	   }


  	const CFuint MAXITER   = 10;
  	const CFreal TOL       = 1e-6;
 	for (CFuint iter = 0; iter < MAXITER; ++iter)
 	{
  	CFreal Restheta        =  std::abs(Theta[0]*TOL);
	  //cout << "ITER" << i  << endl;
 	 //The variables needed for the calculation of Re_thetat
  	GetLambda(Lambda, Theta[0], Viscosity);
        GetFlambda(Lambda,Tu,Flambda,Theta[0],true);  
        GetFlambda(Lambda,Tu,FlambdaPrime,Theta[0],false);  
          
         if (Tu<=1.3) {
         	 CFreal Rethetat0 = (1173.51-589.428*Tu + 0.2196*overTu*overTu);
                 m_Rethetat = Rethetat0 * Flambda;
         }
   	else {
     	      const CFreal lamco5   = Tu - 0.5658;
   	      const CFreal pwtu   = -0.671;
              const CFreal Rethetat0 = 331.5*std::pow(lamco5,pwtu);
              m_Rethetat = Rethetat0 * Flambda;
 	 }

   	const CFreal Rethetatprime =  (m_Rethetat* FlambdaPrime)/Flambda;

     
   	const CFreal MainF = m_Rethetat -(rho*avV)*Theta[0]/mu;
   	const CFreal MainFprime = Rethetatprime -(rho*avV)/mu;
    
    	Theta[1] = Theta[0] - MainF/MainFprime;
   	//cout.precision(20); cout << "diff   " << Theta[1]/Theta[0] << endl;
  	if ( std::abs(Theta[0]-Theta[1]) <= Restheta ) break;
    	Theta[0]= Theta[1];
  	}

return m_Rethetat;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
