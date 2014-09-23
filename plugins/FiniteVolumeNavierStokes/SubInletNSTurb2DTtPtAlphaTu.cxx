#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNavierStokes/SubInletNSTurb2DTtPtAlphaTu.hh"
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

MethodCommandProvider<SubInletNSTurb2DTtPtAlphaTu, CellCenterFVMData, FiniteVolumeNavierStokesModule> 
subInletEuler2DTtPtAlphaFVMCCProvider("SubInletNSTurb2DTtPtAlphaTuFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb2DTtPtAlphaTu::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("Ptot","total pressure");
   options.addConfigOption< CFreal >("alpha","alpha");
}

//////////////////////////////////////////////////////////////////////////////

SubInletNSTurb2DTtPtAlphaTu::SubInletNSTurb2DTtPtAlphaTu(const std::string& name) :
  SubInletEuler2DTtPtAlpha(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _tTotal = 0.0;
   setParameter("Ttot",&_tTotal);

   _pTotal = 0.0;
   setParameter("Ptot",&_pTotal);

   _alpha = 0.0;
   setParameter("alpha",&_alpha);
}

//////////////////////////////////////////////////////////////////////////////

SubInletNSTurb2DTtPtAlphaTu::~SubInletNSTurb2DTtPtAlphaTu()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletNSTurb2DTtPtAlphaTu::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

   ComputeGhostState();
   
  




  // set the ghost state starting from the physical data
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
 }
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
