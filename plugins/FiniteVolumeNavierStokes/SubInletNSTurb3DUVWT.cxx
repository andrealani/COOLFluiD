#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubInletNSTurb3DUVWT.hh"
#include "Framework/MethodCommandProvider.hh"

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

MethodCommandProvider<SubInletNSTurb3DUVWT, CellCenterFVMData, FiniteVolumeNavierStokesModule> subInletNSTurb3DUVWTFVMCCProvider("SubInletNSTurb3DUVWTFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb3DUVWT::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Vx","x velocity");
   options.addConfigOption< CFreal >("Vy","y velocity");
   options.addConfigOption< CFreal >("Vz","z velocity");
   options.addConfigOption< CFreal >("T","static temperature");
   options.addConfigOption< std::vector<CFreal> >("TurbVars","Freestream K, Omega values");
}

//////////////////////////////////////////////////////////////////////////////

SubInletNSTurb3DUVWT::SubInletNSTurb3DUVWT(const std::string& name) :
  FVMCC_BC(name),
  _varSetTurb(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _uinf = 0.0;
   setParameter("Vx",&_uinf);

  _vinf = 0.0;
   setParameter("Vy",&_vinf);

  _winf = 0.0;
   setParameter("Vz",&_winf);

  _temperature = 0.0;
   setParameter("T",&_temperature);

   _turbVars = vector<CFreal>();
   setParameter("TurbVars",&_turbVars);
}

//////////////////////////////////////////////////////////////////////////////

SubInletNSTurb3DUVWT::~SubInletNSTurb3DUVWT()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb3DUVWT::setGhostState(GeometricEntity *const face)
{

  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // set the physical data starting from the inner state
  _varSetTurb->computePhysicalData(*innerState, _dataInnerState);

  // physical constants
  const CFreal gamma = _varSetTurb->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
  const CFreal R = _varSetTurb->getModel()->getR();
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];

  _dataGhostState[EulerTerm::RHO] = pInnerState/(R*_temperature);
  _dataGhostState[EulerTerm::VX] = 2.0*_uinf - _dataInnerState[EulerTerm::VX];
  _dataGhostState[EulerTerm::VY] = 2.0*_vinf - _dataInnerState[EulerTerm::VY];
  _dataGhostState[EulerTerm::VZ] = 2.0*_winf - _dataInnerState[EulerTerm::VZ];
  _dataGhostState[EulerTerm::P] = pInnerState;
  _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
				       _dataGhostState[EulerTerm::VX] +
				       _dataGhostState[EulerTerm::VY]*
				       _dataGhostState[EulerTerm::VY] +
				       _dataGhostState[EulerTerm::VZ]*
				       _dataGhostState[EulerTerm::VZ]);

  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				   + 0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
				       _dataGhostState[EulerTerm::RHO]);

  _dataGhostState[EulerTerm::T] = _temperature;
  _dataGhostState[EulerTerm::E] = _dataGhostState[EulerTerm::H] -
     (_dataGhostState[EulerTerm::P]/_dataGhostState[EulerTerm::RHO]);

  const CFuint iK = _varSetTurb->getModel()->getFirstScalarVar(0);
  const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);
  for(CFuint iTurb = 0; iTurb < nbTurbVars; iTurb++)
  {
    _dataGhostState[iK + iTurb] = 2.0*_turbVars[iTurb] - _dataInnerState[iK + iTurb];
  }

  // set the ghost state starting from the physical data
  _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);

}

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb3DUVWT::setup()
{
  CFAUTOTRACE;

  FVMCC_BC::setup();

  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb3DVarSet>();
  _varSetTurb->getModel()->resizePhysicalData(_dataInnerState);
  _varSetTurb->getModel()->resizePhysicalData(_dataGhostState);

//   if(_turbVars.size() == 0){
// 
//     _turbVars.resize(_varSetTurb->getModel()->getNbScalarVars(0));
// 
//   }
// 
//   //This is for one equation model
//   if(_turbVars.size() == 1){
// 
//   }
// 
//   // for k-Omega model
//   else{
//     cf_assert (_turbVars.size() == 2);
// 
//     ///@todo this is only valid for k-Omega model
//     //velocities are in (m/s), temperature in (Kelvin)
// 
//     //Values taken from: F.R. Menter: Two-Eq Eddy-Viscosity Turbulence Models for Engineering Applications (Aug 1994)
//     const CFreal L = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
//     const CFreal Uinf = sqrt(_uInf*_uInf + _vInf*_vInf + _wInf*_wInf);
//     //Wilcox BC
//     _turbVars[1] = Uinf/L;
//     //Menter BC
//     //_turbVars[1] = 10.*Uinf/L;
// 
//     const CFreal pdim = _pressure * _varSetTurb->getModel()->getPressRef();
//     const CFreal Tdim = _temperature * _varSetTurb->getModel()->getTempRef();
// 
//     const CFreal muInf = _diffVarTurb->getModel().getDynViscosity(pdim, Tdim)/
//                             (*_diffVarTurb->getModel().getReferencePhysicalData())[NSTurbTerm::MU];
//     //upper bound
//     const CFreal muTurbInf = muInf/100.;
//     //lower bound
//     //const CFreal muTurbInf = muInf/100000.;
// 
//     _turbVars[0] = _turbVars[1] * muTurbInf ;
// 
//   }
// 
//   CFout << "Setting the turbulent values to :     k = " << _turbVars[0] << "\n                                  Omega = " << _turbVars[1] << "\n";

  //Check that the initial values for the turbulent variables have been set
  cf_assert (_turbVars.size() == _varSetTurb->getModel()->getNbScalarVars(0));

  _uinf /= _varSetTurb->getModel()->getVelRef();
  _vinf /= _varSetTurb->getModel()->getVelRef();
  _winf /= _varSetTurb->getModel()->getVelRef();
  _temperature /= _varSetTurb->getModel()->getTempRef();

  const CFuint firstScalarVar = _varSetTurb->getModel()->getFirstScalarVar(0);
  const CFuint nbScalarVars = _varSetTurb->getModel()->getNbScalarVars(0);

  const RealVector& refValues = _varSetTurb->getModel()->getReferencePhysicalData();
  for(CFuint iVar=0; iVar < nbScalarVars ;iVar++)
  {
    _turbVars[iVar] /= refValues[firstScalarVar + iVar];
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
