#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubOutletNSTurb3D.hh"
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

MethodCommandProvider<SubOutletNSTurb3D, CellCenterFVMData, FiniteVolumeNavierStokesModule> subOutletNSTurb3DFVMCCProvider("SubOutletNSTurb3DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubOutletNSTurb3D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("P","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletNSTurb3D::SubOutletNSTurb3D(const std::string& name) :
  FVMCC_BC(name),
  _varSetTurb(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _pressure = 1.0;
   setParameter("P",&_pressure);
}

//////////////////////////////////////////////////////////////////////////////

SubOutletNSTurb3D::~SubOutletNSTurb3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletNSTurb3D::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // set the physical data starting from the inner state
  _varSetTurb->computePhysicalData(*innerState, _dataInnerState);

  const CFreal R = _varSetTurb->getModel()->getR();
  const CFreal gamma = _varSetTurb->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);

  _dataGhostState[EulerTerm::VX] = _dataInnerState[EulerTerm::VX];
  _dataGhostState[EulerTerm::VY] = _dataInnerState[EulerTerm::VY];
  _dataGhostState[EulerTerm::VZ] = _dataInnerState[EulerTerm::VZ];
  _dataGhostState[EulerTerm::V]  = _dataInnerState[EulerTerm::V];
  _dataGhostState[EulerTerm::RHO] = _dataInnerState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::P] = 2.0*_pressure - _dataInnerState[EulerTerm::P];
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				   + 0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(_varSetTurb->getModel()->getGamma()*
				       _dataGhostState[EulerTerm::P]/
					 _dataGhostState[EulerTerm::RHO]);

  _dataGhostState[EulerTerm::T] = _dataGhostState[EulerTerm::P]/(R*_dataGhostState[EulerTerm::RHO]);
  _dataGhostState[EulerTerm::E] = _dataGhostState[EulerTerm::H] -
     (_dataGhostState[EulerTerm::P]/_dataGhostState[EulerTerm::RHO]);

  const CFuint iK = _varSetTurb->getModel()->getFirstScalarVar(0);
  const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);
  for(CFuint iTurb = 0; iTurb < nbTurbVars; iTurb++)
  {
    _dataGhostState[iK + iTurb] = _dataInnerState[iK + iTurb];
  }


  _varSetTurb->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletNSTurb3D::setup()
{
  CFAUTOTRACE;

  FVMCC_BC::setup();

  _varSetTurb = getMethodData().getUpdateVar().d_castTo<ConvTurb3DVarSet>();
  _varSetTurb->getModel()->resizePhysicalData(_dataInnerState);
  _varSetTurb->getModel()->resizePhysicalData(_dataGhostState);

  _pressure /= _varSetTurb->getModel()->getPressRef();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
