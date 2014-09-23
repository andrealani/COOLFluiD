#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNavierStokes/SubInletNSTurb2DTtPtRadialAlpha.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

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

MethodCommandProvider<SubInletNSTurb2DTtPtRadialAlpha, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
subInletNSTurb2DTtPtRadialAlphaFVMCCProvider("SubInletNSTurb2DTtPtRadialAlphaFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb2DTtPtRadialAlpha::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("Ptot","total pressure");
}

//////////////////////////////////////////////////////////////////////////////

SubInletNSTurb2DTtPtRadialAlpha::SubInletNSTurb2DTtPtRadialAlpha(const std::string& name) :
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
}

//////////////////////////////////////////////////////////////////////////////

SubInletNSTurb2DTtPtRadialAlpha::~SubInletNSTurb2DTtPtRadialAlpha()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb2DTtPtRadialAlpha::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  // using a normal of interior cell as the ghost cell => reversed sign
  const CFreal FaceLength = -1.*sqrt(nx*nx + ny*ny);
  // unit normal vector of the ghost cell, i.e. pointing towards the domain interior
  nx /= FaceLength;
  ny /= FaceLength;

// if( fabs( sqrt(nx*nx+ny*ny)-1. ) > 0.00001 ) CFout << "\n SOMETHING IS WRONG !!!" ;

  // set the physical data starting from the inner state
  _varSetTurb->computePhysicalData(*innerState, _dataInnerState);

  const CFreal gamma = _varSetTurb->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal R = _varSetTurb->getModel()->getR();
  const CFreal u = _dataInnerState[EulerTerm::VX];
  const CFreal v = _dataInnerState[EulerTerm::VY];
  const CFreal uSqvSq = u*u + v*v;
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];
  const CFreal tInnerState = pInnerState / (R*_dataInnerState[EulerTerm::RHO]);
  const CFreal machInner = sqrt(uSqvSq/(gamma*R*tInnerState));
  const CFreal coefficient = 1.0 + 0.5*gammaMinus1*machInner*machInner;
  const CFreal tTotalInner = tInnerState*coefficient;
  const CFreal coeffPow = pow(coefficient, gamma/gammaMinus1);
  const CFreal pTotalInner = pInnerState*coeffPow;

  // ghost state quantities
  const CFreal tTotalGhost = 2.0*_tTotal - tTotalInner;
  const CFreal pTotalGhost = 2.0*_pTotal - pTotalInner;
  const CFreal machGhost = machInner;
  const CFreal tGhost = tTotalGhost / coefficient;

  // set the physical data for the ghost state
  _dataGhostState[EulerTerm::P]   = pTotalGhost/coeffPow;
  _dataGhostState[EulerTerm::RHO] = _dataGhostState[EulerTerm::P]/(R*tGhost);
  _dataGhostState[EulerTerm::A]   = sqrt(gamma*_dataGhostState[EulerTerm::P]/_dataGhostState[EulerTerm::RHO]);
  _dataGhostState[EulerTerm::T]   = tGhost;

  _dataGhostState[EulerTerm::VX]  = nx * machGhost * _dataGhostState[EulerTerm::A];
  _dataGhostState[EulerTerm::VY]  = ny * machGhost * _dataGhostState[EulerTerm::A];
  _dataGhostState[EulerTerm::V]   = sqrt(_dataGhostState[EulerTerm::VX]*_dataGhostState[EulerTerm::VX] +
				                                 _dataGhostState[EulerTerm::VY]*_dataGhostState[EulerTerm::VY] );
  _dataGhostState[EulerTerm::H]   = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				                            + 0.5*_dataGhostState[EulerTerm::RHO]*
				                            _dataGhostState[EulerTerm::V]*_dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];


    ///@todo this is only valid for k-Omega model
    //check the adimensionalisation - velocities are in (m/s), temperature in (Kelvin)

    //Values taken from: F.R. Menter: Two-Eq Eddy-Viscosity Turbulence Models for Engineering Applications (Aug 1994)
    const CFreal L = PhysicalModelStack::getActive()->getImplementor()->getRefLength() ;
    const CFreal Uinf = _dataGhostState[EulerTerm::V] ;
    //Wilcox BC
//     _turbVars[1] = Uinf/L;
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

// CF_DEBUG_OBJ( _dataGhostState );
  // set the ghost state starting from the physical data
  _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);
// CF_DEBUG_OBJ( *ghostState );
// CF_DEBUG_EXIT;
 }

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb2DTtPtRadialAlpha::setup()
{
  FVMCC_BC::setup();
 
  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb2DVarSet>();
  cf_assert(_varSetTurb.isNotNull());

  _varSetTurb->getModel()->resizePhysicalData(_dataInnerState);
  _varSetTurb->getModel()->resizePhysicalData(_dataGhostState);

  if(_turbVars.size() == 0){
    _turbVars.resize(_varSetTurb->getModel()->getNbScalarVars(0));
  }
// CFout << "\n _turbVars.size() =  " << _turbVars.size() << "\n _varSetTurb->getModel()->getNbScalarVars(0) =  " << _varSetTurb->getModel()->getNbScalarVars(0) ;
  //Check that the initial values for the turbulent variables have been set
//   cf_assert(_turbVars.size() == _varSetTurb->getModel()->getNbScalarVars(0));

  _diffVarTurb = getMethodData().getDiffusiveVar().d_castTo<DiffTurb2DVarSet>();
  cf_assert(_diffVarTurb.isNotNull());

  _pTotal/=_varSetTurb->getModel()->getPressRef();
  _tTotal/=_varSetTurb->getModel()->getTempRef();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
