#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "FiniteVolumeNEQ/SubInletEuler1DTtPtYiTvTe.hh"
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

MethodCommandProvider<SubInletEuler1DTtPtYiTvTe,
                      CellCenterFVMData,
                       FiniteVolumeNEQModule>
subInletEuler1DTtPtYiTvTeFVMCCProvider("SubInletEuler1DTtPtYiTvTeFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler1DTtPtYiTvTe::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("Ptot","total pressure");
   options.addConfigOption< CFreal >("Tv","Vibrational Temperature");
   options.addConfigOption< CFreal >("Te","Electron Temperature");
   options.addConfigOption< std::vector<CFreal> >("Yi","mass fraction");
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler1DTtPtYiTvTe::SubInletEuler1DTtPtYiTvTe(const std::string& name) :
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

   _Tv = 0.0;
   setParameter("Tv",&_Tv);

   _Te = 0.0;
   setParameter("Te",&_Te); 

}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler1DTtPtYiTvTe::~SubInletEuler1DTtPtYiTvTe()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler1DTtPtYiTvTe::setGhostState(GeometricEntity *const face)
{
  State& innerState = *face->getState(0);
  State& ghostState = *face->getState(1);

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(innerState, _dataInnerState);

  SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(library.isNotNull());

  const CFreal Rgas = library->getRgas();
  const CFuint nbSpecies = _varSet->getModel()->getNbScalarVars(0);

  library-> getMolarMasses(_mm);
  CFreal mass = 0;  

  for (CFuint k = 0; k < nbSpecies; k++) {
      mass += _Yi[k]/_mm[k];
  }

  mass = 1.0/mass;
  const CFreal R = Rgas/mass; 
  const CFreal gamma = _dataInnerState[EulerTerm::GAMMA];
  const CFreal gammaMinus1 = gamma - 1.;  
  const CFreal u = _dataInnerState[EulerTerm::VX];
  const CFreal uSqvSq = u*u;
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];
  const CFreal tInnerState = pInnerState /(R*_dataInnerState[EulerTerm::RHO]);
  const CFreal machInner = sqrt(uSqvSq)/(_dataInnerState[EulerTerm::A]);
  const CFreal coefficient = 1.0 + 0.5*gammaMinus1*machInner*machInner;
  const CFreal tTotalInner = tInnerState*coefficient;
  const CFreal coeffPow = pow(coefficient, gamma/gammaMinus1);
  const CFreal pTotalInner = pInnerState*coeffPow;


  // ghost state Total quantities and velocities
  const CFreal tTotalGhost = 2.0*_tTotal - tTotalInner;
  const CFreal pTotalGhost = 2.0*_pTotal - pTotalInner;
  const CFreal machGhost = machInner;
  const CFreal tGhost = tTotalGhost / coefficient;
  
  const CFreal p_Ghost = pTotalGhost/coeffPow;
  const CFreal Rho_Ghost =  p_Ghost/(R*tGhost);
  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal yIGhost = _Yi[i];
    ghostState[i] = yIGhost*Rho_Ghost;
  }
    ghostState[nbSpecies]    = machGhost*sqrt(gamma*R*tGhost);
    ghostState[nbSpecies+1] = tGhost;
    ghostState[nbSpecies+2] = 2.0*_Tv-innerState[nbSpecies+2];   
    ghostState[nbSpecies+3] = 2.0*_Te-innerState[nbSpecies+3];

 }

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler1DTtPtYiTvTe::setup()
{
  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MultiScalarVarSet<Euler1DVarSet> >();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);

  _mm.resize(_varSet->getModel()->getNbScalarVars(0)); 
  _pTotal/=_varSet->getModel()->getPressRef();
  _tTotal/=_varSet->getModel()->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
