#include "SA3DSourceTerm.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "FiniteVolumeTurb/FiniteVolumeSA.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "NavierStokes/EulerVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::SA;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<SA3DSourceTerm,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeSAModule>
SA3DSourceTermFVMCCProvider("SA3DSourceTerm");

//////////////////////////////////////////////////////////////////////////////

void SA3DSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("CompressibilityCorrectionTerm","Add the extra destruction term (Default = False)");
}

//////////////////////////////////////////////////////////////////////////////

SA3DSourceTerm::SA3DSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _varSet(CFNULL),
  _diffVarSet(CFNULL),
  _temp(),
  _avState(),
  _physicalData(),
  _nstates(CFNULL),
  _wallDistance(CFNULL),
  _values()
{
  addConfigOptionsTo(this);
  
  _CompTerm = false;
  setParameter("CompressibilityCorrectionTerm",&_CompTerm);

}

//////////////////////////////////////////////////////////////////////////////

SA3DSourceTerm::~SA3DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void SA3DSourceTerm::setup()
{
  CFAUTOTRACE;

  ComputeSourceTermFVMCC::setup();  
  
  _varSet = getMethodData().getUpdateVar().d_castTo<EulerVarSet>();
  _diffVarSet = getMethodData().getDiffusiveVar().d_castTo<NavierStokes3DSA>();
  cf_assert(_diffVarSet.isNotNull());
  
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());
  _avState.resize(PhysicalModelStack::getActive()->getNbEq());
  
  cf_assert(_varSet.isNotNull());
  _varSet->getModel()->resizePhysicalData(_physicalData);
  
  _nstates = _sockets.getSocketSink<RealVector>("nstates")->getDataHandle();
  _wallDistance = _sockets.getSocketSink<CFreal>("wallDistance")->getDataHandle();

  SafePtr<DerivativeComputer> derComput =
    this->getMethodData().getDerivativeComputer();
  const CFuint nbNodesInControlVolume =
    derComput->getMaxNbVerticesInControlVolume();

  _values.resize(PhysicalModelStack::getActive()->getNbEq(), nbNodesInControlVolume);
}

//////////////////////////////////////////////////////////////////////////////
      
void SA3DSourceTerm::computeSource(Framework::GeometricEntity *const element,
				   RealVector& source,
				   RealMatrix& jacobian)
{

  ///Reset the source term
  source = 0.0;

  const EquationSubSysDescriptor& eqData =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint nbEqs = eqData.getNbEqsSS();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();

  //jacobian contribution
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar(); // whats this
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  // I change this!!!!!!!!!!!!!!!!!!!!!
  
  if((nbEqs != 5) && (isPerturb && iPerturbVar != 5))
  {
    source[5] = _unperturbedNegativePart;
    source[5] += _unperturbedPositivePart;
  }
  const CFuint elemID = element->getID();

  //if we solve using coupling, then no contribution when solving the first block (normal NS)
  if((nbEqs != 5) && ( (isPerturb && iPerturbVar == 5) || (!isPerturb)))
  {
    DataHandle<CFreal> normals = this->socket_normals.getDataHandle(); // what we need the normals???
    DataHandle<CFint> isOutward = this->socket_isOutward.getDataHandle(); // what is this???
    
    cf_assert(_varSet.isNotNull());

    //Set the physical data for the cell considered
    State *const currState = element->getState(0); // what we take from this the rho/P or the first cell/BC
    _varSet->computePhysicalData(*currState, _physicalData);

    //fill in the nodal states
    const CFreal R = _varSet->getModel()->getR();

    //compute the gradients by applying Green Gauss in the
    //cell d's
    CFreal dNutildX = 0.0;
    CFreal dNutildY = 0.0;
    CFreal dNutildZ = 0.0; 
    CFreal dRhodX = 0.0;
    CFreal dRhodY = 0.0;
    CFreal dRhodZ = 0.0; 
	   dVdX = 0.0;
	   dUdY = 0.0;
	   dUdZ = 0.0; 
	   dWdX = 0.0; 
	   dVdZ = 0.0; 
	   dWdY = 0.0; 
    
    // these terms are added for the comp model
	   dUdX = 0.0; 
	   dVdY = 0.0; 
	   dWdZ = 0.0; 
    
    if (this->m_useGradientLS && this->m_gradientsExist) { // what does this check
      const CFuint pID = 0;
      const CFuint uID = 1;
      const CFuint vID = 2;
      const CFuint wID = 3; 
      const CFuint TID = 4;
      const CFuint NutilID = 5; 
      const CFuint start = elemID*totalNbEqs;
      
      dNutildX = this->m_ux[start+NutilID]; 
      dNutildY = this->m_uy[start+NutilID];
      dNutildZ = this->m_uz[start+NutilID]; 
      dVdX = this->m_ux[start+vID];
      dUdY = this->m_uy[start+uID];
      dUdZ = this->m_uz[start+uID]; 
      dWdX = this->m_ux[start+wID]; 
      dVdZ = this->m_uz[start+vID]; 
      dWdY = this->m_uy[start+wID]; 
      
      // these terms are added for the comp model
      dUdX = this->m_ux[start+uID]; 
      dVdY = this->m_uy[start+vID]; 
      dWdZ = this->m_uz[start+wID]; 
      
      // AL: here we assume to have one single species
      // p=rho*R*T  =>  dP=dRho*R*T+rho*R*dT  =>  dRho=(dP-rho*R*dT)/(R*T)
      //const CFreal avRho = _physicalData[EulerTerm::RHO]; // which rho should I use <------------------------------------
      const CFreal avRhoR = _physicalData[EulerTerm::RHO]*R;
      const CFreal RT = R*_physicalData[EulerTerm::T];
      dRhodX = (this->m_ux[start+pID] - avRhoR*this->m_ux[start+TID])/RT; // calculation of the grad(rho)
      dRhodY = (this->m_uy[start+pID] - avRhoR*this->m_uy[start+TID])/RT;
      dRhodZ = (this->m_uz[start+pID] - avRhoR*this->m_uz[start+TID])/RT; // I add this !!!!!!!!!!!!!!
    } 
    else { 
     {
  CFLog(INFO, "Only implemented for Least Square Reconstructor \n");
  throw Common::NotImplementedException (FromHere(),"For !(this->m_useGradientLS && this->m_gradientsExist) not implemented...");
  
}
    }
    
    
    const CFuint iNutil = _varSet->getModel()->getFirstScalarVar(0);
    
    _avState[0] = _physicalData[EulerTerm::P];
    _avState[1] = _physicalData[EulerTerm::VX];
    _avState[2] = _physicalData[EulerTerm::VY];
    _avState[3] = _physicalData[EulerTerm::VZ];
    _avState[4] = _physicalData[EulerTerm::T];
    _avState[5] = _physicalData[iNutil];
    
    const CFreal mu = _diffVarSet->getLaminarDynViscosityFromGradientVars(_avState);
    const CFreal rho = _diffVarSet->getDensity(_avState);
    NIU = mu / rho;
    NIUtilda = _avState[5];
    
    //make sure we don't have negative values of niutilda
    ///@Attention CG: According to Spalart "cliping updates prevents the convergence of discrete PDE residuals and hampers efforts
    //to quantify discrete truncation and solution errors"
    ///@see Allmaras, S. R., Johnson, F. T., and Spalart, P. R., "Modifications and Clarifications for the Implementation 
    ///of the Spalart-Allmaras Turbulence Model," ICCFD7-1902, 7th International Conference on Computational Fluid Dynamics, Big Island, Hawaii, 9-13 July 2012. 
    NIUtilda = max(0.,NIUtilda); // to avoid this use of the SA-negative model. See the above paper.

    // To prevent division by zero (see Ashford, G., An Unstructured Grid
    // Generation and Adaptive Solution Technique for High-Reynolds-Number Compressible Flows, Ph.D. Thesis, University of Michigan 1996.)
    // Page 155
    CFreal Qsi = NIUtilda / NIU;

    //constants of the SA turbulence model:
    const CFreal Cb1 = 0.1355;
    const CFreal Cb2 = 0.622;
    const CFreal sigma = 2.0/3.0;
    const CFreal kappa = 0.41;
    const CFreal Cw1 = ( Cb1 / (kappa*kappa) ) + (( 1. + Cb2)/sigma );
    const CFreal Cw2 = 0.3;
    const CFreal Cw3 = 2.0;
    const CFreal Cv1 = 7.1;
    
    const CFreal Cv2 = 0.7;// these two parameters are used for the modified Stilda
    const CFreal Cv3 = 0.9; ///@see Allmaras, S. R., Johnson, F. T., and Spalart, P. R., "Modifications and Clarifications for the Implementation 
    ///of the Spalart-Allmaras Turbulence Model," ICCFD7-1902, 7th International Conference on Computational Fluid Dynamics, Big Island, Hawaii, 9-13 July 2012. 

    ///Original definition of the vorticity (can be negative):
         const CFreal fv1 = Qsi*Qsi*Qsi / (Qsi*Qsi*Qsi + Cv1*Cv1*Cv1);
         const CFreal fv2 = 1. - ( Qsi /(1. + (Qsi * fv1)));
	 
    /// Calculation of the NIUturbulent
	 NIUturbulent = fv1*NIUtilda;
	 
	 const CFreal dw = SA3DSourceTerm::getDistance(element);
	 
	 const CFreal Nitiloverkapa2d2 = NIUtilda/( kappa*kappa*dw*dw);

	 //const CFreal Nitiloverkapa2d2 = NIUtilda/( kappa*kappa*_d*_d);
	 
         const CFreal S = std::sqrt((dVdX - dUdY)*(dVdX - dUdY) + 
				    (dUdZ - dWdX)*(dUdZ - dWdX) +
				    (dVdZ - dWdY)*(dVdZ - dWdY)); // definition of the 3D vorticity magnitude
	 
	 const CFreal Soverbar = Nitiloverkapa2d2*fv2;
	 
	 CFreal Stilda = 0.; // definition and initialization of Stilda
	 
	 //Preventing Negative Values of Modified Vorticity Stilda
	 ///@see Allmaras, S. R., Johnson, F. T., and Spalart, P. R., "Modifications and Clarifications for the Implementation 
	 ///of the Spalart-Allmaras Turbulence Model," ICCFD7-1902, 7th International Conference on Computational Fluid Dynamics, Big Island, Hawaii, 9-13 July 2012.
	 if (Soverbar >= -(Cv2*S))
	 {
	     Stilda = S + Soverbar;
	 }
	 
	 else
	 {
	     Stilda = S + (S*(Cv2*Cv2*S + Cv3*Soverbar))/((Cv3 - 2.0*Cv2)*S - Soverbar);
	 }
	 
    //Clip the vorticity
    //     Stilda = max(Stilda,1.e-10);
    
    ///@Attention SA - fv3 model should not be used
    /// Modified vorticity definition (always positive)
    // see Ashford, G., An Unstructured Grid Generation and Adaptive Solution
    // Technique for High-Reynolds-Number Compressible Flows, Ph.D. Thesis, University of Michigan 1996.)
    // Page 155
    //const CFreal fv1 = Qsi*Qsi*Qsi / (Qsi*Qsi*Qsi + Cv1*Cv1*Cv1);
    //const CFreal fv2 = pow((1. + (Qsi/Cv2)),-3);
    //const CFreal fv3 = ((1. + (Qsi * fv1))*(1. - fv2))/Qsi; <------------------------------------
    //const CFreal S = fabs(dVdX - dUdY);
    //const CFreal Stilda = (S * fv3) + (( NIUtilda/( kappa*kappa*d*d)) * fv2);
    
    const CFreal rlim = 10.0; // definition of the first SA model
    
    CFreal r = Nitiloverkapa2d2/Stilda;
    
    r = min(r, rlim);
    
    //for the adimensionalization of the maximum value...
    //    RealVector& refDataConv = *_varSet->getModel()->getReferencePhysicalData();
    //     const CFreal rhoRef = refDataConv[EulerTerm::RHO];
    //     const CFreal Uref = refDataConv[EulerTerm::V];
    //     const CFreal Lref = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
    //     CFreal rRef = 1./(rhoRef * Uref*Uref * Lref * Lref);
    //CFreal rRef = 1.;
    
    ///Clip r to 10.0 (original article)
    // r = min(r, 10.*rRef);
    ///Clip r to 2.0 (G. Ashford)
    //r = min(r, 2.*rRef);
    
    const CFreal g = r + Cw2 * ( pow(r,6) - r);
    const CFreal g6 = g*g*g*g*g*g;
    const CFreal Cw3_6 = Cw3*Cw3*Cw3*Cw3*Cw3*Cw3;
    const CFreal sixth = 1./6.;
    const CFreal fw = g * pow( (1. + Cw3_6) / (g6 + Cw3_6) ,sixth);

    /// the below used only for SA - la model
    //const CFreal DU = 0.; // velocity difference point vs.trip point
    //unused // const CFreal DX = 1.; // grid spacing at the trip point
    //unused // const CFreal wt = 1.; // wall vorticity at the trip point
    //unused // const CFreal dt = 1.; // distance from point to trip point
    //unused // const CFreal gt = min(0.1,DU / (wt * DX));
    
    //trip terms
    //     const CFreal ft1 = Ct1 * gt * exp( -1.0 * Ct2 * ((wt*wt)/(DU*DU)) * (d*d + gt*gt*dt*dt) );
    //     const CFreal ft2 = Ct3 * exp(-1.0 * Ct4 * Qsi*Qsi);
    //const CFreal ft1 = 0.0; //override trip term
    //const CFreal ft2 = 0.0; //override trip term
    
    /*
    const CFreal positivePart1 = _diffVarSet->getModel().getCoeffTau() * ( Cb2 / sigma ) * (pow(dKdX,2) + pow(dKdY,2));
    const CFreal positivePart2 = Cb1 * ( 1. - ft2 ) * Stilda * NIUtilda;
    const CFreal positivePart3 = (1./_diffVarSet->getModel().getCoeffTau()) * ft1 * pow(DU,2);
    const CFreal positivePart = (positivePart1 + positivePart2 + positivePart3)*rho;
    */

    /*
    const CFreal negativePart1 = - _diffVarSet->getModel().getCoeffTau() * ( Cw1 * fw ) * ((NIUtilda*NIUtilda) / (d*d));
    const CFreal negativePart2 = _diffVarSet->getModel().getCoeffTau() * ((Cb1/(kappa*kappa))*ft2) * ((NIUtilda*NIUtilda) / (d*d));
    CFreal negativePart = (negativePart1 + negativePart2) * rho;
    
    // correction for the introduction of rho in the convective flux : G
    const CFreal G = ( 1. / sigma ) * (NIU + NIUtilda) * ( dKdX + dKdY ) * ( dRhodX + dRhodY );
    const CFreal adimCoef = _diffVarSet->getModel().getCoeffTau();
    negativePart -= G*adimCoef;
    */
    
    ///@details the notations given from Spalart are used in this implementation. That is
    //P := Production term
    //D := Destruction term
    //T := extra trip term
    
     // it calls the function of the RANS or for the DES modes
	 _d = this->getDistance(element);
    
    const CFreal adimCoef = _diffVarSet->getModel().getCoeffTau();
    
    const CFreal P = Cb1 * Stilda * NIUtilda; // production term
    
    CFreal D = adimCoef * ( Cw1 * fw ) * ((NIUtilda*NIUtilda) / (_d*_d)); // destruction term
    
    // tranform the model into SA - noft2 - comp by adding the extra destruction term
    // It improves the performance of the model in compressible mixing layers
    ///@see the site http://turbmodels.larc.nasa.gov/spalart.html#qcr2000
    if (_CompTerm)
    {
      const CFreal C5 = 3.5;
      
      const CFreal spSound = _physicalData[EulerTerm::A];
      
      const CFreal extraDest = (C5*NIUtilda*NIUtilda/(spSound*spSound))*(SA3DSourceTerm::compSumOfVelocityGrads());
      
      D += extraDest;
    }
    
    const CFreal nonConsDiffTerm = adimCoef * ( Cb2 / sigma ) * (pow(dNutildX,2) + pow(dNutildY,2) + pow(dNutildZ,2));
    
    // correction for the introduction of rho in the convective flux : G
    const CFreal G = adimCoef * ( 1. / sigma ) * (NIU + NIUtilda) * ( dNutildX + dNutildY + dNutildZ ) * ( dRhodX + dRhodY + dRhodZ );
    
    const CFreal negativePart = - (D * rho + G);
    const CFreal positivePart = (P + nonConsDiffTerm )*rho;
    
    
    if(isPerturb)
    {
      cf_assert(iPerturbVar == 5);
      source[5] = negativePart;
      source[5] += _unperturbedPositivePart;
    }
    else
    {
      //rhs contribution
      source[5] = negativePart;
      source[5] += positivePart;
      _unperturbedPositivePart = positivePart;
      _unperturbedNegativePart = negativePart;
    }
  }

  source *= volumes[elemID];

}

///////////////////////////////////////////////////////////////////////////////


///@Attention CG: this method should be overrided from the DES models
CFreal SA3DSourceTerm::getDistance (Framework::GeometricEntity *const element)
{
  //Set the physical data for the cell considered
    State *const currState = element->getState(0); 
    
  CFreal d = _wallDistance[currState->getLocalID()];
    d = max(d,1.e-10);
    
    return d;
}

//////////////////////////////////////////////////////////////////////////////

// CG: this method compute the sum of the velocity gradients. It is needed for the comp model
// and for the DDES nad IDDES modes
 CFreal SA3DSourceTerm::compSumOfVelocityGrads ()
{
  
  CFreal Sum = (dUdX*dUdX + dVdX*dVdX + dWdX*dWdX +
		dUdY*dUdY + dVdY*dVdY + dWdY*dWdY + 
		dUdZ*dUdZ + dVdZ*dVdZ + dWdZ*dWdZ );
  
  return Sum;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
