#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallAdiabaticNSTurb3D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "NavierStokes/NSTurbTerm.hh"
#include "Framework/SubSystemStatus.hh"
#include "NavierStokes/NavierStokes3DVarSet.hh"

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

MethodCommandProvider<NoSlipWallAdiabaticNSTurb3D, CellCenterFVMData,
		      FiniteVolumeNavierStokesModule>
NoSlipWallAdiabaticNSTurb3DFVMCCProvider("NoSlipWallAdiabaticNSTurb3DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNSTurb3D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("xWallVelocity","X-component of a velocity vector of the wall.");
   options.addConfigOption< CFreal >("yWallVelocity","Y-component of a velocity vector of the wall.");
   options.addConfigOption< CFreal >("zWallVelocity","Z-component of a velocity vector of the wall.");
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallAdiabaticNSTurb3D::NoSlipWallAdiabaticNSTurb3D(const std::string& name) :
  FVMCC_BC(name),
  _varSetTurb(CFNULL),
  _diffVarTurb(CFNULL),
  _dataInnerState(),
  _dataGhostState(),
  _xWallVelocity(0.),
  _yWallVelocity(0.),
  _zWallVelocity(0.),
  socket_wallDistance("wallDistance")
{
   addConfigOptionsTo(this);
   setParameter("xWallVelocity",&_xWallVelocity);
   setParameter("yWallVelocity",&_yWallVelocity);
   setParameter("zWallVelocity",&_zWallVelocity);
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallAdiabaticNSTurb3D::~NoSlipWallAdiabaticNSTurb3D()
{
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
NoSlipWallAdiabaticNSTurb3D::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = FVMCC_BC::needsSockets();

  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNSTurb3D::setup()
{
  CFAUTOTRACE;
  
  FVMCC_BC::setup();

  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb3DVarSet>();
  _varSetTurb->getModel()->resizePhysicalData(_dataInnerState);
  _varSetTurb->getModel()->resizePhysicalData(_dataGhostState);
  
  _diffVarTurb = getMethodData().getDiffusiveVar().d_castTo<DiffTurb3DVarSet>();

  _xWallVelocity /= _varSetTurb->getModel()->getVelRef();
  _yWallVelocity /= _varSetTurb->getModel()->getVelRef();
  _zWallVelocity /= _varSetTurb->getModel()->getVelRef(); 
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNSTurb3D::setGhostState(GeometricEntity *const face)
{

  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  _varSetTurb->computePhysicalData(*innerState, _dataInnerState);

  const CFreal R = _varSetTurb->getModel()->getR();
  const CFreal ghostT = _dataInnerState[EulerTerm::P]/
    (R*_dataInnerState[EulerTerm::RHO]);
  const CFreal ghostP = _dataInnerState[EulerTerm::P];
  const CFreal gamma = _varSetTurb->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma - 1.);

  _dataGhostState[EulerTerm::RHO] = _dataInnerState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::VX] = 2.*_xWallVelocity - _dataInnerState[EulerTerm::VX];
  _dataGhostState[EulerTerm::VY] = 2.*_yWallVelocity - _dataInnerState[EulerTerm::VY];
  _dataGhostState[EulerTerm::VZ] = 2.*_zWallVelocity - _dataInnerState[EulerTerm::VZ];
  _dataGhostState[EulerTerm::V] = _dataInnerState[EulerTerm::V];
  _dataGhostState[EulerTerm::P] = ghostP;
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*ghostP +
			           0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/
                                   _dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*ghostP/_dataGhostState[EulerTerm::RHO]);

  _dataGhostState[EulerTerm::T] = ghostT;
  _dataGhostState[EulerTerm::E] = _dataGhostState[EulerTerm::H] -
     (_dataGhostState[EulerTerm::P]/_dataGhostState[EulerTerm::RHO]);

  const CFuint iK = _varSetTurb->getModel()->getFirstScalarVar(0);
  const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);

  //k=0 at the wall
  //(in reality, put a very small value to avoid negative values close to the wall)
  const CFreal kWall = 1.e-20;
  _dataGhostState[iK] = kWall - _dataInnerState[iK];

  //At the wall, according to Menter's definition, Omega_w = 10*((6*nu)/(beta1 * y0 *y0))
  if(nbTurbVars == 2){

      ///Compute y0getDynViscosity
      DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();
      CFreal y0 = wallDistance[innerState->getLocalID()];

      ///This is a temporary fix for the flat plate because, at initialization, the wall Distance is not known!!!!
      ///@todo remove this line...
      if(SubSystemStatusStack::getActive()->getNbIter() == 0) y0 = innerState->getCoordinates()[YY];

      //avoid too small distances
      y0 = max(y0, 10.e-10);

      CFreal pdim = _dataInnerState[EulerTerm::P] * _varSetTurb->getModel()->getPressRef();
      CFreal Tdim = _dataInnerState[EulerTerm::T] * _varSetTurb->getModel()->getTempRef();
      
      CFreal mu = _diffVarTurb->getModel().getDynViscosityDim(pdim, Tdim)/
                      _diffVarTurb->getModel().getReferencePhysicalData()[NSTurbTerm::MU];

      CFreal nu = mu / _dataInnerState[EulerTerm::RHO];

      //this is not the best, but it avoids having to code another BC! because I
      //would have to dynamic cast to the KOmega varset to get the beta1
      CFreal beta1 = 0.075;

      ///@todo here should this be adimensionalized (by the distance)???
      //Menter's definition
      CFreal omegaWall = (10. * 6. * nu) / (beta1 * y0 * y0);

      //Wilcox's definition
      //CFreal omegaWall = (1. * 6. * nu) / (beta1 * y0 * y0);

      _dataGhostState[iK + 1] = 2.0*omegaWall - _dataInnerState[iK + 1];
  }

  _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
