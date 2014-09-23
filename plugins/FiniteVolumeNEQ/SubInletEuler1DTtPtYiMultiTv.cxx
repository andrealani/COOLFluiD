#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "FiniteVolumeNEQ/SubInletEuler1DTtPtYiMultiTv.hh"
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

MethodCommandProvider<SubInletEuler1DTtPtYiMultiTv,
                      CellCenterFVMData,
                       FiniteVolumeNEQModule>
subInletEuler1DTtPtYiMultiTvFVMCCProvider("SubInletEuler1DTtPtYiMultiTvFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler1DTtPtYiMultiTv::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("Ptot","total pressure");
   options.addConfigOption< std::vector<CFreal> >("Tv","Vibrational Temperature");
   options.addConfigOption< std::vector<CFreal> >("Yi","mass fraction");
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler1DTtPtYiMultiTv::SubInletEuler1DTtPtYiMultiTv(const std::string& name) :
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

   _Tv = std::vector<CFreal>();
   setParameter("Tv",&_Tv); 

}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler1DTtPtYiMultiTv::~SubInletEuler1DTtPtYiMultiTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler1DTtPtYiMultiTv::setGhostState(GeometricEntity *const face)
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
  const CFuint nbEvEqs = _varSet->getModel()->getNbScalarVars(1);

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
  for (CFuint k = 0;k < nbEvEqs;k++ ){
      ghostState[nbSpecies+2+k] = 2.0*_Tv[k]-innerState[nbSpecies+2+k];
  }

 }

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler1DTtPtYiMultiTv::setup()
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
