#include "SA2DSourceTerm.hh"
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

MethodStrategyProvider<SA2DSourceTerm,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeSAModule>
SA2DSourceTermFVMCCProvider("SA2DSourceTerm");

//////////////////////////////////////////////////////////////////////////////

SA2DSourceTerm::SA2DSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _varSet(CFNULL),
  _diffVarSet(CFNULL),
  _temp(),
  _avState(),
  _physicalData(),
  _nstates(CFNULL),
  _wallDistance(CFNULL),
  _values(),
  _states(),
  _rho()
{
}

//////////////////////////////////////////////////////////////////////////////

SA2DSourceTerm::~SA2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void SA2DSourceTerm::setup()
{
  CFAUTOTRACE;

  ComputeSourceTermFVMCC::setup();  
  
  _varSet = getMethodData().getUpdateVar().d_castTo<EulerVarSet>();
  _diffVarSet = getMethodData().getDiffusiveVar().d_castTo<NavierStokes2DSA>();
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
  _states.reserve(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////
      
void SA2DSourceTerm::computeSource(Framework::GeometricEntity *const element,
				   RealVector& source,
				   RealMatrix& jacobian)
{

  ///Reset the source term
  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = 0.0;
  source[4] = 0.0;

  const EquationSubSysDescriptor& eqData =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint nbEqs = eqData.getNbEqsSS();

  //jacobian contribution
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  if((nbEqs != 4) && (isPerturb && iPerturbVar != 4))
  {
    source[4] = _unperturbedNegativePart;
    source[4] += _unperturbedPositivePart;
  }
  const CFuint elemID = element->getID();

  //if we solve using coupling, then no contribution when solving the first block (normal NS)
  if((nbEqs != 4) && ( (isPerturb && iPerturbVar == 4) || (!isPerturb)))
  {
    DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
    DataHandle<CFint> isOutward = this->socket_isOutward.getDataHandle();
    
    cf_assert(_varSet.isNotNull());

    //Set the physical data for the cell considered
    State *const currState = element->getState(0);
    _varSet->computePhysicalData(*currState, _physicalData);

    //fill in the nodal states
    const vector<Node*>* const nodes = element->getNodes();
    const CFuint nbNodesInElem = nodes->size();
    _states.clear();
    for (CFuint i = 0; i < nbNodesInElem; ++i) {
      _states.push_back(&_nstates[(*nodes)[i]->getLocalID()]);
    }

    _rho.resize(nbNodesInElem);

    _diffVarSet->setGradientVars(_states, _values, _states.size());

    const CFreal R = _varSet->getModel()->getR();
    for(CFuint iNode = 0; iNode < nbNodesInElem ; iNode++)
    {
      _rho[iNode] = _values(0,iNode)/(R*_values(3,iNode));
    }

    const vector<GeometricEntity*>& faces = *element->getNeighborGeos();
    cf_assert(faces.size() == nbNodesInElem);

    //compute the gradients by applying Green Gauss in the
    //cell d's
    CFreal dKdX = 0.0;
    CFreal dKdY = 0.0;

    CFreal dRhodX = 0.0;
    CFreal dRhodY = 0.0;

    CFreal dVdX = 0.0;
    CFreal dUdY = 0.0;

    for (CFuint i = 0; i < nbNodesInElem; ++i) {
      //get the face normal
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
        dKdX += nx*(_values(4,i) + _values(4,i+1));
        dKdY += ny*(_values(4,i) + _values(4,i+1));

        dRhodX += nx*(_rho[i] + _rho[i+1]);
        dRhodY += ny*(_rho[i] + _rho[i+1]);

        dVdX += nx*(_values(2,i) + _values(2,i+1));
        dUdY += ny*(_values(1,i) + _values(1,i+1));
      }
      else {
        dKdX += nx*(_values(4,i) + _values(4,0));
        dKdY += ny*(_values(4,i) + _values(4,0));

        dRhodX += nx*(_rho[i] + _rho[0]);
        dRhodY += ny*(_rho[i] + _rho[0]);

        dVdX += nx*(_values(2,i) + _values(2,0));
        dUdY += ny*(_values(1,i) + _values(1,0));
      }
    }
    dKdX *= 0.5/volumes[elemID];
    dKdY *= 0.5/volumes[elemID];
    dRhodX *= 0.5/volumes[elemID];
    dRhodY *= 0.5/volumes[elemID];
    dUdY *= 0.5/volumes[elemID];
    dVdX *= 0.5/volumes[elemID];

    //compute PUVTK by averaging the nodes (since it is not computed in the
    //physical data, would be generally wrong computing it from the equation
    //of perfect gases ... think about the LTE case)
//     CFreal avP = 0.;
//     CFreal avU = 0.;
//     CFreal avV = 0.;
//     CFreal avT = 0.;
//     CFreal avK = 0.;
     CFreal avDist = _wallDistance[currState->getLocalID()];
//     for (CFuint i = 0; i < nbNodesInElem; ++i) {
//       avP += _p[i];
//       avU += _u[i];
//       avV += _v[i];
//       avT += _T[i];
//       avK += _K[i];
//     }
//     avP /= nbNodesInElem;
//     avU /= nbNodesInElem;
//     avV /= nbNodesInElem;
//     avT /= nbNodesInElem;
//     avK /= nbNodesInElem;
//
//     _avState[0] = avP;
//     _avState[1] = avU;
//     _avState[2] = avV;
//     _avState[3] = avT;
//     _avState[4] = avK;

  const CFuint iK = _varSet->getModel()->getFirstScalarVar(0);

  _avState[0] = _physicalData[EulerTerm::P];
  _avState[1] = _physicalData[EulerTerm::VX];
  _avState[2] = _physicalData[EulerTerm::VY];
  _avState[3] = _physicalData[EulerTerm::T];
  _avState[4] = _physicalData[iK];



    const CFreal mu = _diffVarSet->getLaminarDynViscosityFromGradientVars(_avState);
    const CFreal rho = _diffVarSet->getDensity(_avState);
    const CFreal NIU = mu / rho;
    CFreal NIUtilda = _avState[4];

    //make sure we don't have negative values of niutilda
    NIUtilda = max(0.,NIUtilda);

    // To prevent division by zero (see Ashford, G., An Unstructured Grid
    // Generation and Adaptive Solution Technique for High-Reynolds-Number Compressible Flows, Ph.D. Thesis, University of Michigan 1996.)
    // Page 155
    CFreal Qsi = NIUtilda / NIU;
    Qsi = max (Qsi, 0.001);

    //constants of the SA turbulence model:
    const CFreal Cb1 = 0.1355;
    const CFreal Cb2 = 0.622;
    const CFreal sigma = 2.0/3.0;
    const CFreal kappa = 0.41;
    const CFreal Cw1 = ( Cb1 / (kappa*kappa) ) + (( 1. + Cb2)/sigma );
    const CFreal Cw2 = 0.3;
    const CFreal Cw3 = 2.0;
    const CFreal Cv1 = 7.1;
    const CFreal Cv2 = 5.0;
    //unused // const CFreal Ct1 = 1.0;
    //unused // const CFreal Ct2 = 2.0;
    //unused // const CFreal Ct3 = 1.2;
    //unused // const CFreal Ct4 = 0.5;

    //Clip the distance from the cell to the closest wall
    CFreal d = avDist;
    d = max(d,1.e-10);

    ///Original definition of the vorticity (can be negative):
//     const CFreal fv1 = Qsi*Qsi*Qsi / (Qsi*Qsi*Qsi + Cv1*Cv1*Cv1);
//     const CFreal fv2 = 1. - ( Qsi /(1. + (Qsi * fv1)));
//     const CFreal S = fabs(dVdX - dUdY);
//     CFreal Stilda = S + (( NIUtilda/( kappa*kappa*d*d)) * fv2);
    //Clip the vorticity
//     Stilda = max(Stilda,1.e-10);

    /// Modified vorticity definition (always positive)
    // see Ashford, G., An Unstructured Grid Generation and Adaptive Solution
    // Technique for High-Reynolds-Number Compressible Flows, Ph.D. Thesis, University of Michigan 1996.)
    // Page 155
    const CFreal fv1 = Qsi*Qsi*Qsi / (Qsi*Qsi*Qsi + Cv1*Cv1*Cv1);
    const CFreal fv2 = pow((1. + (Qsi/Cv2)),-3);
    const CFreal fv3 = ((1. + (Qsi * fv1))*(1. - fv2))/Qsi;
    const CFreal S = fabs(dVdX - dUdY);
    const CFreal Stilda = (S * fv3) + (( NIUtilda/( kappa*kappa*d*d)) * fv2);

    CFreal r = NIUtilda / (Stilda * kappa * kappa * d * d );

    //for the adimensionalization of the maximum value...
//    RealVector& refDataConv = *_varSet->getModel()->getReferencePhysicalData();
//     const CFreal rhoRef = refDataConv[EulerTerm::RHO];
//     const CFreal Uref = refDataConv[EulerTerm::V];
//     const CFreal Lref = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
//     CFreal rRef = 1./(rhoRef * Uref*Uref * Lref * Lref);
    CFreal rRef = 1.;

    ///Clip r to 10.0 (original article)
    // r = min(r, 10.*rRef);
    ///Clip r to 2.0 (G. Ashford)
    r = min(r, 2.*rRef);

    const CFreal g = r + Cw2 * ( pow(r,6) - r);
    const CFreal g6 = g*g*g*g*g*g;
    const CFreal Cw3_6 = Cw3*Cw3*Cw3*Cw3*Cw3*Cw3;
    const CFreal sixth = 1./6.;
    const CFreal fw = g * pow( (1. + Cw3_6) / (g6 + Cw3_6) ,sixth);

    const CFreal DU = 0.; // velocity diference point vs.trip point
    //unused // const CFreal DX = 1.; // grid spacing at the trip point
    //unused // const CFreal wt = 1.; // wall vorticity at the trip point
    //unused // const CFreal dt = 1.; // distance from point to trip point
    //unused // const CFreal gt = min(0.1,DU / (wt * DX));

    //trip terms
//     const CFreal ft1 = Ct1 * gt * exp( -1.0 * Ct2 * ((wt*wt)/(DU*DU)) * (d*d + gt*gt*dt*dt) );
//     const CFreal ft2 = Ct3 * exp(-1.0 * Ct4 * Qsi*Qsi);
    const CFreal ft1 = 0.0; //override trip term
    const CFreal ft2 = 0.0; //override trip term

    const CFreal positivePart1 = _diffVarSet->getModel().getCoeffTau() * ( Cb2 / sigma ) * (pow(dKdX,2) + pow(dKdY,2));
    const CFreal positivePart2 = Cb1 * ( 1. - ft2 ) * Stilda * NIUtilda;
    const CFreal positivePart3 = (1./_diffVarSet->getModel().getCoeffTau()) * ft1 * pow(DU,2);

    CFreal positivePart = (positivePart1 + positivePart2 + positivePart3)*rho;

    const CFreal negativePart1 = - _diffVarSet->getModel().getCoeffTau() * ( Cw1 * fw ) * ((NIUtilda*NIUtilda) / (d*d));
    const CFreal negativePart2 = _diffVarSet->getModel().getCoeffTau() * ((Cb1/(kappa*kappa))*ft2) * ((NIUtilda*NIUtilda) / (d*d));

    CFreal negativePart = (negativePart1 + negativePart2) * rho;

    //correction for the introduction of rho in the convective flux : G
    const CFreal G = ( 1. / sigma ) * (NIU + NIUtilda) * ( dKdX + dKdY ) * ( dRhodX + dRhodY );
    const CFreal adimCoef = _diffVarSet->getModel().getCoeffTau();
    negativePart -= G*adimCoef;

    if(isPerturb)
    {
      cf_assert(iPerturbVar == 4);
      source[4] = negativePart;
      source[4] += _unperturbedPositivePart;
    }
    else
    {
      //rhs contribution
      source[4] = negativePart;
      source[4] += positivePart;
      _unperturbedPositivePart = positivePart;
      _unperturbedNegativePart = negativePart;
    }
  }

  source *= volumes[elemID];

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
