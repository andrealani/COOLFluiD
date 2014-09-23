#include "Framework/MethodStrategyProvider.hh"
#include "Framework/DataHandle.hh"
#include "Framework/GeometricEntity.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/MHD3DProjectionHeatingSourceTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<MHD3DProjectionHeatingSourceTerm,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mHD3DProjectionHeatingSTFVMCCProvider("MHD3DProjectionHeatingST");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionHeatingSourceTerm::MHD3DProjectionHeatingSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _varSet(CFNULL),
  _physicalData(),
  _dataLeftState(),
  _dataRightState()
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionHeatingSourceTerm::~MHD3DProjectionHeatingSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionHeatingSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  cf_assert(_varSet.isNotNull());
  
  _varSet->getModel()->resizePhysicalData(_physicalData);
  _varSet->getModel()->resizePhysicalData(_dataLeftState);
  _varSet->getModel()->resizePhysicalData(_dataRightState);
  
  DataHandle<State*, GLOBAL> state = _globalSockets.getSocketSink<State*>("states")->getDataHandle();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionHeatingSourceTerm::computeSource(Framework::GeometricEntity *const element,
					  RealVector& source,
					  RealMatrix& jacobian)
{
  cf_assert(_varSet.isNotNull());
  
  DataHandle<State*, GLOBAL> state = _globalSockets.getSocketSink<State*>("states")->getDataHandle();
  DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
  DataHandle<CFint> isOutward = this->socket_isOutward.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  // this source term is for MHD flows
  const vector<State*>* const states = element->getStates();
  const CFuint elementID = element->getID();

  // all elements in FVM should have only one state
  cf_assert(states->size() == 1);

  State *const currState = (*states)[0];

  RealVector stateCoordsCartesian = currState->getCoordinates();
  RealVector stateCoordsSpherical(PhysicalModelStack::getActive()->getDim());
  RealMatrix carSphTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim()),
        sphCarTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim());

  _varSet->setTransformationMatrices(stateCoordsCartesian,stateCoordsSpherical,carSphTransMat,sphCarTransMat);

  const CFreal r = stateCoordsSpherical[0];
  const CFreal theta = stateCoordsSpherical[1];

  const CFreal th1 = (MathTools::MathConsts::CFrealPi()/180.0)*17.5;
  const CFreal th2 = (MathTools::MathConsts::CFrealPi()/180.0)*61.5;

  // critical angle
  CFreal critTheta = 0.;
  if (r <= 7.0) {
	critTheta = asin(sqrt(sin(th1)*sin(th1)+cos(th1)*cos(th1)*(r-1.0)/8.0));
  }
  if ((r > 7.0) && (r < 30.0)) {
	critTheta = asin(sqrt(sin(th2)*sin(th2)+cos(th2)*cos(th2)*(r-7.0)/40.0));
  }

  // volumetric heating amplitude (J/(kg*s*K))
  const CFreal q0 = 100.0;

  // reference length = RSun
  const CFreal lRef = _varSet->getLRef();

  // target temperature and the heating scale height
  CFreal T0 = 0.;
  CFreal sigma0 = 0.;
  if (fabs(0.5*MathTools::MathConsts::CFrealPi()-theta) < critTheta) {
    T0 = 1.5e6;
    sigma0 = 4.5*lRef;
  }
  if (fabs(0.5*MathTools::MathConsts::CFrealPi()-theta) > critTheta) {
    T0 = 2.63e6;
    sigma0 = 4.5*lRef*(2.0-(sin(theta)*sin(theta)/(sin(critTheta)*sin(critTheta))));
  }
  
  _varSet->computePhysicalData(*currState, _physicalData);

  // Boltzmann constant
  const CFreal k = 1.3806503e-23;

  // mass of proton
  const CFreal mp = 1.67262158e-27;

  // mass of electron
  const CFreal me = 9.10938188e-31;

  const CFreal nRef = _varSet->getNRef();
  const CFreal TRef = _varSet->getTRef(); 

  const CFreal rhoRef = nRef*(mp+me);
  const CFreal vRef = sqrt(2.0*k*TRef/mp);

  const CFreal T = _physicalData[MHDTerm::P]*(mp+me)*vRef*vRef/(_physicalData[MHDTerm::RHO]*k);

  // dimensional volumetric heating term (J/(m^3*s))
  const CFreal Q = _physicalData[MHDTerm::RHO]*rhoRef*q0*(T0-T)*exp(-(r-1.0)*(r-1.0)/(sigma0*sigma0));

  const CFreal nondimconst = lRef/(rhoRef*vRef*vRef*vRef);

  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = 0.0;
  source[4] = 0.0;
  source[5] = 0.0;
  source[6] = 0.0;
  source[7] = Q*nondimconst*volumes[elementID];
  source[8] = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
