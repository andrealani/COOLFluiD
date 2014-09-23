#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "NavierStokes/Euler1DVarSet.hh"
#include "FiniteVolumeNEQ/SubOutletEuler1DT.hh"
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

MethodCommandProvider<SubOutletEuler1DT,
                      CellCenterFVMData,
                       FiniteVolumeNEQModule>
subOutletEuler1DTFVMCCProvider("SubOutletEuler1DTFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler1DT::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Tout","outlet pressure");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler1DT::SubOutletEuler1DT(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _library(CFNULL),
  _dataInnerState()
{
   addConfigOptionsTo(this);

   _tOut = 0.0;
   setParameter("Tout",&_tOut);
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler1DT::~SubOutletEuler1DT()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler1DT::setGhostState(GeometricEntity *const face)
{
  State& innerState = *face->getState(0);
  State& ghostState = *face->getState(1);
   
  // set the physical data starting from the inner state
  _varSet->computePhysicalData(innerState, _dataInnerState);

  SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(library.isNotNull());

  const CFuint nbSpecies = _varSet->getModel()->getNbScalarVars(0);
  const CFuint nbTv = _varSet->getModel()->getNbScalarVars(1);

  const CFreal uGhost = innerState[nbSpecies];
  const CFreal tInner = innerState[nbSpecies+1];
  const CFreal tGhost = 2.*_tOut - tInner;

  ghostState[nbSpecies] = uGhost;
  ghostState[nbSpecies+1] = tGhost;
 
  for (CFuint is = 0; is < nbSpecies; is++) {
       ghostState[is] = innerState[is]*tInner/tGhost;
  }
 
  if (nbTv != 0) { 
      for (CFuint k = 0; k < nbTv; k++) {
           ghostState[nbSpecies+2+k] = innerState[nbSpecies+2+k]; 
      }
  }
  
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler1DT::setup()
{
  FVMCC_BC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<MultiScalarVarSet<Euler1DVarSet> >();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);

  _tOut/=_varSet->getModel()->getTempRef();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
