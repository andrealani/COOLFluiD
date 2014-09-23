#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubOutletEuler2D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubOutletEuler2D, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
subOutletEuler2DFVMCCProvider("SubOutletEuler2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("P","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler2D::SubOutletEuler2D(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _pressure = 1.0;
   setParameter("P",&_pressure);
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler2D::~SubOutletEuler2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2D::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;

  const CFreal u = _dataInnerState[EulerTerm::VX];
  const CFreal v = _dataInnerState[EulerTerm::VY];
  // unused // const CFreal rho = _dataInnerState[EulerTerm::RHO];
  const CFreal vn = u*nx + v*ny;
  // unused // const CFreal pInnerState = _dataInnerState[EulerTerm::P];
  const CFreal aInnerState = _dataInnerState[EulerTerm::A];
  const CFreal machInner = vn / aInnerState;

  // depending if the outlet is subsonic or supersonic
  // we impose pressure or not

  // supersonic outlet case
  if (machInner >= 1.0) {
    (*ghostState)[0] = (*innerState)[0];
    (*ghostState)[1] = (*innerState)[1];
    (*ghostState)[2] = (*innerState)[2];
    (*ghostState)[3] = (*innerState)[3];
  }
  // subsonic outlet case
  else{
    _dataGhostState[EulerTerm::VX] = _dataInnerState[EulerTerm::VX];
    _dataGhostState[EulerTerm::VY] = _dataInnerState[EulerTerm::VY];
    _dataGhostState[EulerTerm::V] = _dataInnerState[EulerTerm::V];
    _dataGhostState[EulerTerm::RHO] = _dataInnerState[EulerTerm::RHO];
    _dataGhostState[EulerTerm::P] = 2.0*_pressure - _dataInnerState[EulerTerm::P];
    _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
                                    + 0.5*_dataGhostState[EulerTerm::RHO]*
                                    _dataGhostState[EulerTerm::V]*
                                    _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
    _dataGhostState[EulerTerm::A] = sqrt(_varSet->getModel()->getGamma()*
                                        _dataGhostState[EulerTerm::P]/
                                        _dataGhostState[EulerTerm::RHO]);
     
     _dataGhostState[EulerTerm::T] = _dataInnerState[EulerTerm::T];

    _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
  }

}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2D::setup()
{

  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

  _pressure /= _varSet->getModel()->getPressRef();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
