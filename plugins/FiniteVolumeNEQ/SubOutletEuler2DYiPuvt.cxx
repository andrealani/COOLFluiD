#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNEQ/SubOutletEuler2DYiPuvt.hh"
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

MethodCommandProvider<SubOutletEuler2DYiPuvt, CellCenterFVMData, FiniteVolumeNEQModule>
subOutletEuler2DYiPuvtFVMCCProvider("SubOutletEuler2DYiPuvt");


//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DYiPuvt::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("Pout","pressure compressible");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler2DYiPuvt::SubOutletEuler2DYiPuvt(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _library(CFNULL),
 _dataInnerState()

{
  addConfigOptionsTo(this);
  m_pressure = 0.0;
  setParameter("Pout",&m_pressure);
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler2DYiPuvt::~SubOutletEuler2DYiPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DYiPuvt::setGhostState(GeometricEntity *const face)
{
  State& innerState = *face->getState(0);
  State& ghostState = *face->getState(1);
  
  // set the physical data starting from the inner state
  _varSet->computePhysicalData(innerState, _dataInnerState);
  
  SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(library.isNotNull());
  
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];
  const CFreal Rgas = library->getRgas();
  const CFuint nbSpecies = _varSet->getModel()->getNbScalarVars(0);
  const CFuint Tghost =  _dataInnerState[EulerTerm::T];
  
 library-> getMolarMasses(_mm);
 
 const CFuint firstSpecies = _varSet->getModel()->getFirstScalarVar(0);
 CFreal mass = 0;
 for (CFuint k = 0; k < nbSpecies; k++) {
   _Yi[k] = _dataInnerState[firstSpecies + k];
   mass += _Yi[k]/_mm[k];
 }
 
 //cout << "mass" << "" <<  mass << endl;
 
 const CFuint sum = mass*Rgas*Tghost;
 
 // first variable is the 
 const CFreal pGhost    = 2.*m_pressure  - pInnerState;
 const CFreal Rho_Ghost = pGhost /sum; 
 
 for (CFuint i = 0; i < nbSpecies; ++i) {
   const CFreal yIGhost = _Yi[i];
   ghostState[i] = yIGhost*Rho_Ghost;
 }
 
 // following variables are exptrapolated
 ghostState[nbSpecies]    = _dataInnerState[EulerTerm::VX];
 ghostState[nbSpecies+1]  = _dataInnerState[EulerTerm::VY];
 ghostState[nbSpecies+2]  = _dataInnerState[EulerTerm::T];

 if (_library->getNbTempVib() == 1) {
   ghostState[nbSpecies+3]  = innerState[nbSpecies+3]; 
 }
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DYiPuvt::setup()
{
   //cout << "mm" <<  endl;
   FVMCC_BC::setup();
  _varSet = getMethodData().getUpdateVar().d_castTo<MultiScalarVarSet<Euler2DVarSet> >();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _mm.resize(_varSet->getModel()->getNbScalarVars(0));
  _Yi.resize(_varSet->getModel()->getNbScalarVars(0));
   //m_pressure/=_varSet->getModel()->getPressRef();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
