#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubInletNSTurb3DTtPtRadialAlpha.hh"
#include "Framework/MethodCommandProvider.hh"
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

MethodCommandProvider<SubInletNSTurb3DTtPtRadialAlpha, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
subInletNSTurb3DTtPtRadialAlphaFVMCCProvider("SubInletNSTurb3DTtPtRadialAlphaFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb3DTtPtRadialAlpha::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("Ptot","total pressure");
   options.addConfigOption< bool   >("Swirl","contains swirling velocity component");
   options.addConfigOption< CFreal >("SwirlIntensity","swirling intensity 0->zero u-velocity, 1->only u-velocity");
}

//////////////////////////////////////////////////////////////////////////////

SubInletNSTurb3DTtPtRadialAlpha::SubInletNSTurb3DTtPtRadialAlpha(const std::string& name) :
  FVMCC_BC(name),
  _varSetTurb(CFNULL),
  _diffVarTurb(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _tTotal = 0.0;
   setParameter("Ttot",&_tTotal);

   _pTotal = 0.0;
   setParameter("Ptot",&_pTotal);

   _isSwirling = false;
   setParameter("Swirl",&_isSwirling);

   _swirl = 0.0;
   setParameter("SwirlIntensity",&_swirl);
}

//////////////////////////////////////////////////////////////////////////////

SubInletNSTurb3DTtPtRadialAlpha::~SubInletNSTurb3DTtPtRadialAlpha()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb3DTtPtRadialAlpha::setGhostState(GeometricEntity *const face)
{

  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = normals[startID + 2];
  // using a normal of interior cell as the ghost cell => reversed sign
  const CFreal FaceLength = -1.*sqrt(nx*nx + ny*ny + nz*nz);
  // unit normal vector of the ghost cell, i.e. pointing towards the domain interior
  nx /= FaceLength;
  ny /= FaceLength;
  nz /= FaceLength;

if( fabs( sqrt(nx*nx+ny*ny+nz*nz)-1. ) > 0.00001 ) CFout << "\n SOMETHING IS WRONG !!!" ;

  // set the physical data starting from the inner state
  _varSetTurb->computePhysicalData(*innerState, _dataInnerState);

  // physical constants
  const CFreal gamma = _varSetTurb->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal R = _varSetTurb->getModel()->getR();

  // inner state
  const CFreal u = _dataInnerState[EulerTerm::VX];
  const CFreal v = _dataInnerState[EulerTerm::VY];
  const CFreal w = _dataInnerState[EulerTerm::VZ];
  const CFreal u2v2w2 = u*u + v*v + w*w;
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];
  const CFreal tInnerState = _dataInnerState[EulerTerm::P] / (R*_dataInnerState[EulerTerm::RHO]);
  const CFreal machInner = sqrt( u2v2w2 / (gamma*R*tInnerState) );
  const CFreal coefficient = 1.0 + 0.5*gammaMinus1*machInner*machInner;
  const CFreal tTotalInner = tInnerState*coefficient;
  const CFreal coeffPow = pow(coefficient, gamma/gammaMinus1);
  const CFreal pTotalInner = pInnerState*coeffPow;

  // ghost state
  const CFreal tTotalGhost = 2.0*_tTotal - tTotalInner; // total temperature
  const CFreal pTotalGhost = 2.0*_pTotal - pTotalInner; // total pressure
  const CFreal machGhost = machInner;
  const CFreal tGhost = tTotalGhost / coefficient; // temperature

  // physical data of the ghost state
  _dataGhostState[EulerTerm::P]   = pTotalGhost/coeffPow;
  _dataGhostState[EulerTerm::RHO] = _dataGhostState[EulerTerm::P]/(R*tGhost);
  _dataGhostState[EulerTerm::A]   = sqrt(gamma*_dataGhostState[EulerTerm::P]/_dataGhostState[EulerTerm::RHO]);
  _dataGhostState[EulerTerm::T] = tGhost;
  _dataGhostState[EulerTerm::VX]  = nx * machGhost * _dataGhostState[EulerTerm::A];
  _dataGhostState[EulerTerm::VY]  = ny * machGhost * _dataGhostState[EulerTerm::A];
  _dataGhostState[EulerTerm::VZ]  = nz * machGhost * _dataGhostState[EulerTerm::A];

  if( _isSwirling == true ){
    cf_assert( fabs( _swirl ) < 1.0 );

    RealVector& coord = ghostState->getCoordinates();
    const CFreal coorX = coord[0];
    const CFreal coorY = coord[1];
    const CFreal coorZ = coord[2];

    const CFreal A = coorY*coorY + coorZ*coorZ;
    const CFreal B = 2.0*_swirl*coorX*coorZ*u;
    const CFreal C = _swirl*_swirl*u*u*(coorX*coorX + coorY*coorY) - coorY*coorY*u2v2w2;
    const CFreal Det = B*B - 4.0*A*C;
    cf_assert( Det >= 0. );

    if( coorY > 0. ){ _dataGhostState[EulerTerm::VZ] = ( -B + sqrt(Det) ) / ( 2.0*A ); }
    else            { _dataGhostState[EulerTerm::VZ] = ( -B - sqrt(Det) ) / ( 2.0*A ); }

    _dataGhostState[EulerTerm::VY] = -1.0 * ( coorZ * _dataGhostState[EulerTerm::VZ] + _swirl*coorX*u ) / coorY;
    _dataGhostState[EulerTerm::VX] *= _swirl ;
  }

  _dataGhostState[EulerTerm::V]   = sqrt(_dataGhostState[EulerTerm::VX]*_dataGhostState[EulerTerm::VX] +
                                         _dataGhostState[EulerTerm::VY]*_dataGhostState[EulerTerm::VY] +
                                         _dataGhostState[EulerTerm::VZ]*_dataGhostState[EulerTerm::VZ] );
  _dataGhostState[EulerTerm::H]   = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
                                    + 0.5*_dataGhostState[EulerTerm::RHO]*
                                    _dataGhostState[EulerTerm::V]*_dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];

    //Values taken from: F.R. Menter: Two-Eq Eddy-Viscosity Turbulence Models for Engineering Applications (Aug 1994)
    const CFreal L = PhysicalModelStack::getActive()->getImplementor()->getRefLength() ;
    const CFreal Uinf = _dataGhostState[EulerTerm::V] ;
    //Wilcox BC
    //   _turbVars[1] = Uinf/L;
    //Menter BC
    _turbVars[1] = 10.*Uinf/L;

    const CFreal pdim = _dataGhostState[EulerTerm::P] * _varSetTurb->getModel()->getPressRef();
    const CFreal Tdim = _dataGhostState[EulerTerm::P] / ( _dataGhostState[EulerTerm::RHO] * R ) * _varSetTurb->getModel()->getTempRef();
    
    const CFreal muInf = _diffVarTurb->getModel().getDynViscosityDim(pdim, Tdim)/
      (_diffVarTurb->getModel().getReferencePhysicalData())[NSTurbTerm::MU];
        
    //upper bound
    const CFreal muTurbInf = muInf/100.;
    //lower bound
    //const CFreal muTurbInf = muInf/100000.;

    _turbVars[0] = _turbVars[1] * muTurbInf ;

  const CFuint firstTurbVar = _varSetTurb->getModel()->getFirstScalarVar(0);

    _dataGhostState[firstTurbVar]   = _turbVars[0];
    _dataGhostState[firstTurbVar+1] = _turbVars[1];

  // set the ghost state starting from the physical data
  _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);

}

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb3DTtPtRadialAlpha::setup()
{
  FVMCC_BC::setup();

  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb3DVarSet>();
  cf_assert(_varSetTurb.isNotNull());

  _varSetTurb->getModel()->resizePhysicalData(_dataInnerState);
  _varSetTurb->getModel()->resizePhysicalData(_dataGhostState);

  if(_turbVars.size() == 0){
    _turbVars.resize(_varSetTurb->getModel()->getNbScalarVars(0));
  }

  _diffVarTurb = getMethodData().getDiffusiveVar().d_castTo<DiffTurb3DVarSet>();
  cf_assert(_diffVarTurb.isNotNull());

  _pTotal/=_varSetTurb->getModel()->getPressRef();
  _tTotal/=_varSetTurb->getModel()->getTempRef();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
