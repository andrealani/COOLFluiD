#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallAdiabaticNSTurb2D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "NavierStokes/NSTurbTerm.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NoSlipWallAdiabaticNSTurb2D, CellCenterFVMData,
		      FiniteVolumeNavierStokesModule>
NoSlipWallAdiabaticNSTurb2DFVMCCProvider("NoSlipWallAdiabaticNSTurb2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNSTurb2D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("yWallVelocity","Y-component of a velocity vector of the wall.");
   options.addConfigOption< CFreal >("xWallVelocity","X-component of a velocity vector of the wall.");
   options.addConfigOption< CFreal >("KWall","Wall value for turbulent intensity");
   options.addConfigOption< CFreal >("KGhostMin","Minimum turb. intensity in the ghost state");
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallAdiabaticNSTurb2D::NoSlipWallAdiabaticNSTurb2D(const std::string& name) :
  FVMCC_BC(name),
  _varSetTurb(CFNULL),
  _diffVarTurb(CFNULL),
  _xWallVelocity(0.),
  _yWallVelocity(0.),
  m_kID(4),
  m_ghostK(0.)
{
  addConfigOptionsTo(this);
  setParameter("xWallVelocity",&_xWallVelocity);
  setParameter("yWallVelocity",&_yWallVelocity);

  m_wallK = 1.e-8;
  setParameter("KWall",&m_wallK);

  m_ghostKMin = 0.0;
  setParameter("KGhostMin",&m_ghostKMin);
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallAdiabaticNSTurb2D::~NoSlipWallAdiabaticNSTurb2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNSTurb2D::setup()
{
  CFAUTOTRACE;
  
  FVMCC_BC::setup();
  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb2DVarSet>();
  _diffVarTurb = getMethodData().getDiffusiveVar().d_castTo<DiffTurb2DVarSet>();
  
  _xWallVelocity /= _varSetTurb->getModel()->getVelRef();
  _yWallVelocity /= _varSetTurb->getModel()->getVelRef();

  cf_assert(m_wallK >= 0.0);
}
      
//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNSTurb2D::setGhostState(GeometricEntity *const face)
{
  computeGhostPosition(face);
  setGhostStateImpl(*face->getState(0), *face->getState(1));
}


//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNSTurb2D::setGhostStateImpl
(const Framework::State& innerState, Framework::State& ghostState)
{
  // here a fix is needed in order to have always m_ghostK > 0
  // dynamic relocation of the ghost state: the position of the
  // ghost state is locally changed, and the BC is imposed
  // using a weighted average of ghost state (in the new location)
  // and inner state
  const CFreal innerK = innerState[m_kID];
  this->repositionNode(innerK, m_ghostK, m_wallK, m_ghostKMin);
  
  const CFreal innerP = innerState[0];
  const CFreal innerT = innerState[3];
  const CFreal innerDensity = innerP / (_varSetTurb->getModel()->getR() * innerT);
  //   const CFreal innerSoundSpeed = sqrt(_varSetTurb->getModel()->getR()*_varSetTurb->getModel()->getGamma()*innerT);
  //   const CFreal innerVel = sqrt(innerU*innerU + innerV*innerV);
  //   const CFreal normalVel = m_faceNormal[0]*innerU + m_faceNormal[1]*innerV;
  
  //Modification of the pressure to avoid pressure oscillations at the wall
  //see article from M-S Liou, A sequel to AUSM: AUSM+, JCP 129 364-382 (1996)
  //Here, the velocity is not the normal velocity because it relates to the u+-a eigen values
//   CFreal pressureCorrection = sqrt(innerVel*innerVel) * innerDensity * innerSoundSpeed;
  CFreal wallPressure = innerP;
  //   if(normalVel < 0.) wallPressure += pressureCorrection ;
  //   else wallPressure -= pressureCorrection;

  cf_assert(innerP > 0.);
  cf_assert(innerT > 0.);

  // reset the ghost node with the new position
  ghostState.getCoordinates() = m_tempGhostNode;
  
  // here we assume that we are in PUVTKOmega
  ghostState[0] = wallPressure;
  this->linearInterpolate(innerState[1], _xWallVelocity, ghostState[1]);
  this->linearInterpolate(innerState[2], _yWallVelocity, ghostState[2]);
  ghostState[3] = innerState[3];
  ghostState[4] = m_ghostK;
  
  const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);
  if(nbTurbVars == 2){

      //Compute distance to innerstate
      CFreal y0 = m_drXiXw;

      //avoid too small distances
      y0 = max(y0, 10.e-10);

      const CFreal pdim =  innerState[0] * _varSetTurb->getModel()->getPressRef();
      const CFreal Tdim =  innerState[3] * _varSetTurb->getModel()->getTempRef();

      const CFreal mu = _diffVarTurb->getModel().getDynViscosityDim(pdim, Tdim)/
                      (_diffVarTurb->getModel().getReferencePhysicalData())[NSTurbTerm::MU];

      CFreal nu = mu / innerDensity;

      //this is not the best, but it avoids having to code another BC! because I
      //would have to dynamic cast to the KOmega varset to get the beta1
      const CFreal beta1 = 0.075;

      ///@todo here should this be adimensionalized (by the distance)???
      //Menter's definition
      const CFreal omegaWall = (10. * 6. * nu) / (beta1 * y0 * y0);

      //Wilcox's definition
      //CFreal omegaWall = (1. * 6. * nu) / (beta1 * y0 * y0);
      this->linearInterpolate(innerState[5], omegaWall, ghostState[5]);

      if(ghostState[5] < 0.) {
       ghostState[5] = innerState[5];
      }
   }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
