#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNEQ/SubInletEuler2DTtPtAlphaYiMultiTv.hh"
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

MethodCommandProvider<SubInletEuler2DTtPtAlphaYiMultiTv,
                      CellCenterFVMData,
                       FiniteVolumeNEQModule>
subInletEuler2DTtPtAlphaYiMultiTvFVMCCProvider("SubInletEuler2DTtPtAlphaYiMultiTv");

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DTtPtAlphaYiMultiTv::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("Ptot","total pressure");
   options.addConfigOption< CFreal >("alpha","alpha");
   options.addConfigOption< std::vector<CFreal> >("Tv","Vibrational Temperature");
   options.addConfigOption< std::vector<CFreal> >("Yi","mass fraction");
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler2DTtPtAlphaYiMultiTv::SubInletEuler2DTtPtAlphaYiMultiTv(const std::string& name) :
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

   _alpha = 0.0;
   setParameter("alpha",&_alpha);   

   _Yi = std::vector<CFreal>();
   setParameter("Yi",&_Yi);

   _Tv = std::vector<CFreal>();
   setParameter("Tv",&_Tv);

}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler2DTtPtAlphaYiMultiTv::~SubInletEuler2DTtPtAlphaYiMultiTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DTtPtAlphaYiMultiTv::setGhostState(GeometricEntity *const face)
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
   const CFuint nbEvEqs = _varSet->getModel()->getNbScalarVars(1);

  const CFreal gamma = _dataInnerState[EulerTerm::GAMMA];
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal R = Rgas;
  const CFreal u = _dataInnerState[EulerTerm::VX];
  const CFreal v = _dataInnerState[EulerTerm::VY];
  const CFreal uSqvSq = u*u + v*v;
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];
  const CFreal tInnerState = pInnerState /(R*_dataInnerState[EulerTerm::RHO]);
  const CFreal machInner = sqrt(uSqvSq)/(_dataInnerState[EulerTerm::A]);
  const CFreal coefficient = 1.0 + 0.5*gammaMinus1*machInner*machInner;
  const CFreal tTotalInner = tInnerState*coefficient;
  const CFreal coeffPow = pow(coefficient, gamma/gammaMinus1);
  const CFreal pTotalInner = pInnerState*coeffPow;
  const CFreal tgAlphaInner = v/u;

  // ghost state Total quantities
  const CFreal tTotalGhost = 2.0*_tTotal - tTotalInner;
  const CFreal pTotalGhost = 2.0*_pTotal - pTotalInner;
  const CFreal tgAlphaGhost = 2.0*tan(_alpha) - tgAlphaInner;
  const CFreal machGhost = machInner;
  const CFreal tGhost = tTotalGhost / coefficient;


  const CFreal p_Ghost = pTotalGhost/coeffPow;
  const CFreal Rho_Ghost =  p_Ghost/(R*tGhost);
  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal yIGhost = _Yi[i];
    ghostState[i] = yIGhost*Rho_Ghost;
  }
      ghostState[nbSpecies]    = machGhost*sqrt(gamma*R*tGhost/
                                        (1.0 + tgAlphaGhost*tgAlphaGhost));
  
      ghostState[nbSpecies+1] = tgAlphaGhost*ghostState[nbSpecies];
      ghostState[nbSpecies+2] = tGhost;
   
      for (CFuint k = 0;k < nbEvEqs;k++ ){
      ghostState[nbSpecies+3+k] = 2.0*_Tv[k]-innerState[nbSpecies+3+k];
  }

 }

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DTtPtAlphaYiMultiTv::setup()
{
  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MultiScalarVarSet<Euler2DVarSet> >();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);

  _pTotal/=_varSet->getModel()->getPressRef();
  _tTotal/=_varSet->getModel()->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
