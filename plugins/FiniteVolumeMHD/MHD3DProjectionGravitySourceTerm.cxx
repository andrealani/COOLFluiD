#include "Framework/DataHandle.hh"
#include "Framework/GeometricEntity.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/MHD3DProjectionGravitySourceTerm.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<MHD3DProjectionGravitySourceTerm,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mHD3DProjectionGravitySTFVMCCProvider("MHD3DProjectionGravityST");

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionGravitySourceTerm::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("xMassCenter","x-coordinate of the center of mass of the external object for which the gravitational force will be computed.");
   options.addConfigOption< CFreal >("yMassCenter","y-coordinate of the center of mass of the external object for which the gravitational force will be computed.");
   options.addConfigOption< CFreal >("zMassCenter","z-coordinate of the center of mass of the external object for which the gravitational force will be computed.");
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionGravitySourceTerm::MHD3DProjectionGravitySourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _varSet(CFNULL),
  _physicalData()
{
  addConfigOptionsTo(this);

  _xMassCenter = 0.0;
  setParameter("xMassCenter",&_xMassCenter);

  _yMassCenter = 0.0;
  setParameter("yMassCenter",&_yMassCenter);

  _zMassCenter = 0.0;
  setParameter("zMassCenter",&_zMassCenter);

}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionGravitySourceTerm::~MHD3DProjectionGravitySourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionGravitySourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  cf_assert(_varSet.isNotNull());

  _varSet->getModel()->resizePhysicalData(_physicalData);

  DataHandle<State*, GLOBAL> state = _globalSockets.getSocketSink<State*>("states")->getDataHandle();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionGravitySourceTerm::computeSource(GeometricEntity *const element,
 					  RealVector& source,
					  RealMatrix& jacobian)
{
 cf_assert(_varSet.isNotNull());

 DataHandle<State*, GLOBAL> state = _globalSockets.getSocketSink<State*>("states")->getDataHandle();
 DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
 // this source term is for MHD flows
 const vector<State*>* const states = element->getStates();
 const CFuint elementID = element->getID();

 // all elements in FVM should have only one state
 cf_assert(states->size() == 1);
 
 State *const currState = (*states)[0];

 _varSet->computePhysicalData(*currState, _physicalData);
   
 const RealVector& stateCoords = currState->getCoordinates();
 const CFreal xCoord = stateCoords[0];
 const CFreal yCoord = stateCoords[1];
 const CFreal zCoord = stateCoords[2];

 // distance from the center of mass of the external object by default assumed to be at the origin
 const CFreal distMass = sqrt((xCoord-_xMassCenter)*(xCoord-_xMassCenter)
	+(yCoord-_yMassCenter)*(yCoord-_yMassCenter)
	+(zCoord-_zMassCenter)*(zCoord-_zMassCenter));
 
 // gravitational constant
 const CFreal G = 6.67384e-11;

 // Boltzmann constant
 const CFreal k = 1.3806503e-23;

 // mass of proton
 const CFreal mp = 1.67262158e-27;

 const CFreal TRef = _varSet->getTRef();
 const CFreal lRef = _varSet->getLRef();

 const CFreal vRef = sqrt(2.0*k*TRef/mp);

 // mass of the external object
 const CFreal mass = _varSet->getMass();

 const CFreal gx = -G*mass*(xCoord-_xMassCenter)/(lRef*lRef*distMass*distMass*distMass);
 const CFreal gy = -G*mass*(yCoord-_yMassCenter)/(lRef*lRef*distMass*distMass*distMass);
 const CFreal gz = -G*mass*(zCoord-_zMassCenter)/(lRef*lRef*distMass*distMass*distMass);

 const CFreal nondimconst = lRef/(vRef*vRef);

 const CFreal Vdotg = _physicalData[MHDTerm::VX]*gx
	+ _physicalData[MHDTerm::VY]*gy
	+ _physicalData[MHDTerm::VZ]*gz;

 source[0] = 0.0;
 source[1] = _physicalData[MHDTerm::RHO]*(gx*nondimconst)*volumes[elementID];
 source[2] = _physicalData[MHDTerm::RHO]*(gy*nondimconst)*volumes[elementID];
 source[3] = _physicalData[MHDTerm::RHO]*(gz*nondimconst)*volumes[elementID];
 source[4] = 0.0;
 source[5] = 0.0;
 source[6] = 0.0;
 source[7] = _physicalData[MHDTerm::RHO]*(Vdotg*nondimconst)*volumes[elementID];
 source[8] = 0.0;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
