#include "Framework/MethodStrategyProvider.hh"
#include "Framework/DataHandle.hh"
#include "Framework/GeometricEntity.hh"
#include "MHD/MHD3DVarSet.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/MHD3DCoriolisCentrifugalForcesSourceTerm.hh"
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

MethodStrategyProvider<MHD3DCoriolisCentrifugalForcesSourceTerm,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mHD3DCoriolisCentrifugalForcesSTFVMCCProvider("MHD3DCoriolisCentrifugalForcesST");

//////////////////////////////////////////////////////////////////////////////

void MHD3DCoriolisCentrifugalForcesSourceTerm::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("OmegaX","x-component of the angular velocity of the external object to which the reference frame is attached.");
   options.addConfigOption< CFreal >("OmegaY","y-component of the angular velocity of the external object to which the reference frame is attached.");
   options.addConfigOption< CFreal >("OmegaZ","z-component of the angular velocity of the external object to which the reference frame is attached.");
}

//////////////////////////////////////////////////////////////////////////////

MHD3DCoriolisCentrifugalForcesSourceTerm::MHD3DCoriolisCentrifugalForcesSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _varSet(CFNULL),
  _physicalData(),
  _dataLeftState(),
  _dataRightState()
{
  addConfigOptionsTo(this);

  // the x-component of the angular velocity of the Sun at its radius on the equator is assigned by default 
  _OmegaX = 0.0;
  setParameter("OmegaX",&_OmegaX);

  // the y-component of the angular velocity of the Sun at its radius on the equator is assigned by default 
  _OmegaY = 0.0;
  setParameter("OmegaY",&_OmegaY);

  // the z-component of the angular velocity of the Sun at its radius on the equator is assigned by default 
  _OmegaZ = 2.0*MathTools::MathConsts::CFrealPi()/(27.2753*24.0*3600.0);
  setParameter("OmegaZ",&_OmegaZ);
}

//////////////////////////////////////////////////////////////////////////////

MHD3DCoriolisCentrifugalForcesSourceTerm::~MHD3DCoriolisCentrifugalForcesSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DCoriolisCentrifugalForcesSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DVarSet>();
  cf_assert(_varSet.isNotNull());
  
  _varSet->getModel()->resizePhysicalData(_physicalData);
  _varSet->getModel()->resizePhysicalData(_dataLeftState);
  _varSet->getModel()->resizePhysicalData(_dataRightState);
  
  DataHandle<State*, GLOBAL> state = _globalSockets.getSocketSink<State*>("states")->getDataHandle();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DCoriolisCentrifugalForcesSourceTerm::computeSource(Framework::GeometricEntity *const element,
					  RealVector& source,
					  RealMatrix& jacobian)
{
  cf_assert(_varSet.isNotNull());
  
  DataHandle<State*, GLOBAL> state = _globalSockets.getSocketSink<State*>("states")->getDataHandle();
  DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
  DataHandle<CFint> isOutward =  this->socket_isOutward.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  // this source term is for MHD flows
  const vector<State*>* const states = element->getStates();
  const CFuint elementID = element->getID();

  // all elements in FVM should have only one state
  cf_assert(states->size() == 1);

  State *const currState = (*states)[0];

  RealVector stateCoord = currState->getCoordinates();

  const CFreal xCoord = stateCoord[0];
  const CFreal yCoord = stateCoord[1];
  /// !!! TBD: The centrifugal force is computed by assuming that the Sun (or the rotating object) rotates around the z-axis. No tilt is allowed yet (i.e. no non-zero OmegaX or OmegaY). However, the formulation below is general and rotation axis tilting will be implemented in the future when necessary.
  //const CFreal zCoord = stateCoord[2];
  const CFreal zCoord = 0.0;

  _varSet->computePhysicalData(*currState, _physicalData);

  // magnetic permeability at vacuum
  const CFreal mu0 = 4.e-7*MathTools::MathConsts::CFrealPi();

  // mass of proton
  const CFreal mp = 1.67262158e-27;

  // mass of electron
  const CFreal me = 9.10938188e-31;

  const CFreal nRef = _varSet->getNRef();
  const CFreal BRef = _varSet->getBRef();
  const CFreal lRef = _varSet->getLRef();

  const CFreal rhoRef = nRef*(mp+me);
  const CFreal vRef = BRef / sqrt(mu0*rhoRef); 

  // non-dimensionalize the angular velocity components
  const CFreal OmegaX = _OmegaX*lRef/vRef;
  const CFreal OmegaY = _OmegaY*lRef/vRef;
  const CFreal OmegaZ = _OmegaZ*lRef/vRef;

  const CFreal coriolisX = -2.0*_physicalData[MHDTerm::RHO]*
		(OmegaY*_physicalData[MHDTerm::VZ]-OmegaZ*_physicalData[MHDTerm::VY]);
  const CFreal coriolisY = -2.0*_physicalData[MHDTerm::RHO]*
		(OmegaZ*_physicalData[MHDTerm::VX]-OmegaX*_physicalData[MHDTerm::VZ]);
  const CFreal coriolisZ = -2.0*_physicalData[MHDTerm::RHO]*
                (OmegaX*_physicalData[MHDTerm::VY]-OmegaY*_physicalData[MHDTerm::VX]);
  const CFreal centrifugalX = -_physicalData[MHDTerm::RHO]*
		(OmegaX*OmegaY*yCoord - OmegaY*OmegaY*xCoord - OmegaZ*OmegaZ*xCoord + OmegaX*OmegaZ*zCoord);
  const CFreal centrifugalY = -_physicalData[MHDTerm::RHO]*
                (OmegaY*OmegaZ*zCoord - OmegaZ*OmegaZ*yCoord - OmegaX*OmegaX*yCoord + OmegaX*OmegaY*xCoord);
  const CFreal centrifugalZ = -_physicalData[MHDTerm::RHO]*
                (OmegaX*OmegaZ*xCoord - OmegaX*OmegaX*zCoord - OmegaY*OmegaY*zCoord + OmegaY*OmegaZ*yCoord);

  const CFreal Vdotcentrifugal = _physicalData[MHDTerm::VX]*centrifugalX+_physicalData[MHDTerm::VY]*centrifugalY
                        +_physicalData[MHDTerm::VZ]*centrifugalZ;

  source[0] = 0.0;
  source[1] = (coriolisX+centrifugalX)*volumes[elementID];
  source[2] = (coriolisY+centrifugalY)*volumes[elementID];
  source[3] = (coriolisZ+centrifugalZ)*volumes[elementID];
  source[4] = 0.0;
  source[5] = 0.0;
  source[6] = 0.0;
  source[7] = Vdotcentrifugal*volumes[elementID];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
