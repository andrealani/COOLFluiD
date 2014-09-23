#include "Framework/MethodStrategyProvider.hh"
#include "Framework/DataHandle.hh"
#include "Framework/GeometricEntity.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/MHD3DProjectionIoTorusSourceTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<MHD3DProjectionIoTorusSourceTerm,
		       CellCenterFVMData,
		       Framework::ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mHD3DProjectionIoTorusSTFVMCCProvider("MHD3DProjectionIoTorusST");

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionIoTorusSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("distanceTorusCenterFromOrigin","Distance of center of the torus from the origin.");
  options.addConfigOption< CFreal >("radiusTorus","Radius of torus.");
  options.addConfigOption< CFreal >("velTangentialTorus","Corotation velocity of the torus (tangential).");
  options.addConfigOption< CFreal >("dRhodtTorusPlasma","Rate of change of plasma density in the torus.");
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionIoTorusSourceTerm::MHD3DProjectionIoTorusSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _varSet(CFNULL),
  _stateCoord(CFNULL)
{
  addConfigOptionsTo(this);

  _a0 = 5.9;
  setParameter("distanceTorusCenterFromOrigin",&_a0);

  _radiusTorus = 1.0;
  setParameter("radiusTorus",&_radiusTorus);

  _velTangentialTorus = 1.0;
  setParameter("velTangentialTorus",&_velTangentialTorus);

  _dRhodtTorusPlasma = 1.0;
  setParameter("dRhodtTorusPlasma",&_dRhodtTorusPlasma);
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionIoTorusSourceTerm::~MHD3DProjectionIoTorusSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionIoTorusSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  cf_assert(_varSet.isNotNull());

  _stateCoord.resize(Framework::PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionIoTorusSourceTerm::computeSource(Framework::GeometricEntity *const element,
						     RealVector& source,
						     RealMatrix& jacobian)
{
  cf_assert(_varSet.isNotNull());

  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  // this source term is for MHD flows
  const vector<State*>* const states = element->getStates();
  const CFuint elementID = element->getID();

  // all elements in FVM should have only one state
  cf_assert(states->size() == 1);

  State *const currState = (*states)[0];

  _stateCoord = currState->getCoordinates();
  const CFreal mX = _varSet->getModel()->getMX();
  const CFreal mY = _varSet->getModel()->getMY();
  const CFreal mZ = _varSet->getModel()->getMZ();

  CFreal uTorusPlasma = 0.0, vTorusPlasma = 0.0, wTorusPlasma = 0.0, rhoChangeTorusPlasma = 0.0;

  // one case is implemented for testing purposes

  if ((mX == 0.0) && (mY == 0.0) && (mZ != 0.0)) {
    const CFreal a = sqrt(_stateCoord[0]*_stateCoord[0] +
		  	  _stateCoord[1]*_stateCoord[1]);
    const CFreal pi = 3.141592653;
    const CFreal r = sqrt((a-_a0)*(a-_a0) + _stateCoord[2]*_stateCoord[2]);
    if (r <= _radiusTorus) {
      uTorusPlasma = _velTangentialTorus*_stateCoord[1]/a;
      vTorusPlasma = -_velTangentialTorus*_stateCoord[0]/a;
      rhoChangeTorusPlasma = _dRhodtTorusPlasma*0.5*(1.0-cos(pi*(_radiusTorus-r)/_radiusTorus));
    }
  }
  
  // const CFreal velocityTorusPlasmaSq = uTorusPlasma*uTorusPlasma +
  //  vTorusPlasma*vTorusPlasma +
  //  wTorusPlasma*wTorusPlasma;
  
  source[0] = rhoChangeTorusPlasma*volumes[elementID];
  source[1] = rhoChangeTorusPlasma*uTorusPlasma*volumes[elementID];
  source[2] = rhoChangeTorusPlasma*vTorusPlasma*volumes[elementID];
  source[3] = rhoChangeTorusPlasma*wTorusPlasma*volumes[elementID];
  source[4] = 0.0;
  source[5] = 0.0;
  source[6] = 0.0;
  // TO BE DONE ! the rate of change pressure due to the extra plasma, dpExtra/dt, should be calculated if possible, make it a user-defined parameter in the CFcase file
  //source[7] = dpExtradt/gammaMinus1 + 0.5*rhoChangeTorusPlasma*velocityTorusPlasmaSq*volumes[elementID];
  source[8] = 0.0;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
