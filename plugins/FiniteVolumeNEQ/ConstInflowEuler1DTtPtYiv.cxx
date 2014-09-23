#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "FiniteVolumeNEQ/ConstInflowEuler1DTtPtYiv.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

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

MethodCommandProvider<ConstInflowEuler1DTtPtYiv,
                      CellCenterFVMData,
                       FiniteVolumeNEQModule>
constInflowEuler1DTtPtYivFVMCCProvider("ConstInflowEuler1DTtPtYivFVMCC");

//////////////////////////////////////////////////////////////////////////////

void ConstInflowEuler1DTtPtYiv::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("Ptot","total pressure");
   options.addConfigOption< std::vector<CFreal> >("Yi","mass fraction");
   options.addConfigOption< CFreal >("u","velocity");
}

//////////////////////////////////////////////////////////////////////////////

ConstInflowEuler1DTtPtYiv::ConstInflowEuler1DTtPtYiv(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _library(CFNULL),
  _dataInnerState()
{
   addConfigOptionsTo(this);
   _tTotal = 0.0;
   setParameter("Ttot",&_tTotal);

   _pTotal = 0.0;
   setParameter("Ptot",&_pTotal);

   _Yi = std::vector<CFreal>();
   setParameter("Yi",&_Yi);

   _u = 0.0;
   setParameter("u",&_u);

}

//////////////////////////////////////////////////////////////////////////////

ConstInflowEuler1DTtPtYiv::~ConstInflowEuler1DTtPtYiv()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConstInflowEuler1DTtPtYiv::setGhostState(GeometricEntity *const face)
{
  State& innerState = *face->getState(0);
  State& ghostState = *face->getState(1);

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(innerState, _dataInnerState);

  SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(library.isNotNull());

  // const CFreal Rgas = library->getRgas();
   const CFreal Rgas =_varSet->getModel()->getR();
   const CFuint nbSpecies = _varSet->getModel()->getNbScalarVars(0);

  const CFreal gamma = _dataInnerState[EulerTerm::GAMMA];
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal R = Rgas;
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];
  const CFreal tInnerState = pInnerState /(R*_dataInnerState[EulerTerm::RHO]);
  const CFreal uInnerState = _dataInnerState[EulerTerm::VX];
  const CFreal uSqvSq = uInnerState*uInnerState;
  const CFreal machInner = sqrt(uSqvSq)/(_dataInnerState[EulerTerm::A]);
  const CFreal coefficient = 1.0 + 0.5*gammaMinus1*machInner*machInner;
  const CFreal tTotalInner = tInnerState*coefficient;
  const CFreal coeffPow = pow(coefficient, gamma/gammaMinus1);
  const CFreal pTotalInner = pInnerState*coeffPow;


  // ghost state Total quantities and velocities
  const CFreal tTotalGhost = 2.0*_tTotal - tTotalInner;
  const CFreal pTotalGhost = 2.0*_pTotal - pTotalInner;
  const CFreal uGhost = 2.0*_u - uInnerState;
  const CFreal tGhost = tTotalGhost / coefficient;
  const CFreal p_Ghost = pTotalGhost/coeffPow;
  const CFreal Rho_Ghost =  p_Ghost/(R*tGhost);

  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal yIGhost = _Yi[i];
    ghostState[i] = yIGhost*Rho_Ghost;
  }
      ghostState[nbSpecies] = uGhost;
      ghostState[nbSpecies+1] = tGhost;

 }

//////////////////////////////////////////////////////////////////////////////

void ConstInflowEuler1DTtPtYiv::setup()
{
  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MultiScalarVarSet<Euler1DVarSet> >();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);

  _pTotal/=_varSet->getModel()->getPressRef();
  _tTotal/=_varSet->getModel()->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
