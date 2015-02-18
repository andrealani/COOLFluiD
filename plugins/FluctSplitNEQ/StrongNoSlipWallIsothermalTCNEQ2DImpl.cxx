#include "FluctSplitNEQ/FluctSplitNEQ.hh"
#include "FluctSplitNEQ/StrongNoSlipWallIsothermalTCNEQ2DImpl.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFL.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::FluctSplit;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongNoSlipWallIsothermalTCNEQ2DImpl,
		      FluctuationSplitData,
		      FluctSplitNEQModule>
strongNoSlipWallIsothermalTCNEQ2DImplProvider("StrongNoSlipWallIsothermalTCNEQ2DImpl");

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalTCNEQ2DImpl::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("TWall","Wall temperature.");
  options.addConfigOption<CFuint, Config::DynamicOption<> >
    ("NbIterAdiabatic", "Number of iterations to run adiabatic");
}
      
//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallIsothermalTCNEQ2DImpl::StrongNoSlipWallIsothermalTCNEQ2DImpl(const std::string& name) :
  StrongImplBC(name)
{
  addConfigOptionsTo(this);
  
  _TWall = 0.0;
  setParameter("TWall",&_TWall);
  
  _nbIterAdiabatic = 0;
  setParameter("NbIterAdiabatic",&_nbIterAdiabatic);
}

//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallIsothermalTCNEQ2DImpl::~StrongNoSlipWallIsothermalTCNEQ2DImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalTCNEQ2DImpl::setup()
{
  StrongImplBC::setup();
  
  cf_assert(_TWall > 0.);
  
  SafePtr<MultiScalarTerm<EulerTerm> > term = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<MultiScalarTerm<EulerTerm> >();
  _TWall /= term->getTempRef();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbSpecies = term->getNbScalarVars(0);
  const CFuint nbSpPlus1 = nbSpecies+1; 
  const CFuint nbSpPlus2 = nbSpecies+2; 
  const CFuint nbSpPlus3 = nbSpecies+3;
  
  const CFuint nbTvs = term->getNbScalarVars(1);

  // apply always the initial condition
  std::vector< Common::SafePtr<TopologicalRegionSet> >& trsList = getTrsList();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< vector<bool> > discardTimeJacob = socket_discardTimeJacob.getDataHandle();
  
  for(CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> trs = trsList[iTRS];
   
    CFLogDebugMax("StrongNoSlipWallIsothermalTCNEQ2DImpl::setup() called for TRS: "
		  << trs->getName() << "\n");

    Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();
        
    for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
      const CFuint localStateID = (*statesIdx)[iState];
      State *const currState = states[localStateID];
      // fix the wall temperature (assume rho_i v T Tv variables)
      
      (*currState)[nbSpecies] = 0.0;
      (*currState)[nbSpPlus1] = 0.0;
      
      if (_nbIterAdiabatic == 0) {
	(*currState)[nbSpPlus2] = _TWall;
	if ( nbTvs > 0 ){
	  (*currState)[nbSpPlus3] = _TWall;
	}
      }
      
      /// @TODO AL: here check that the order in which BCs are setup
      /// is the order of their execution !!!
      discardTimeJacob[localStateID].resize(nbEqs);
      for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
	discardTimeJacob[localStateID][iVar] = (iVar >= nbSpecies) ? true : false;
      }
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalTCNEQ2DImpl::executeOnTrs()
{
  const CFuint nbIter = SubSystemStatusStack::getActive()->getNbIter();
  if (_nbIterAdiabatic == 0 || nbIter > _nbIterAdiabatic) {
    computeIsothermal();
  }
  else {
    computeAdiabatic();
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalTCNEQ2DImpl::computeIsothermal()
{
  CFLog(VERBOSE, "StrongNoSlipWallIsothermalTCNEQ2DImpl::computeIsothermal() START\n");
  
  SafePtr<LinearSystemSolver> lss = getMethodData().getLinearSystemSolver()[0];
  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();
  
  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "StrongNoSlipWallIsothermalTCNEQ2DImpl::execute() called for TRS: "
		 << trs->getName() << "\n");
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle< std::valarray<State*> > bStatesNeighbors =
    socket_bStatesNeighbors.getDataHandle();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal invcfl = 1./getMethodData().getCFL()->getCFLValue();
  const LSSIdxMapping& idxMapping =
    getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();
  
  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();
  SafePtr<MultiScalarTerm<EulerTerm> > term = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<MultiScalarTerm<EulerTerm> >();
  
  const CFuint nbSpecies = term->getNbScalarVars(0);
  const CFuint nbSpPlus1 = nbSpecies+1; 
  const CFuint nbSpPlus2 = nbSpecies+2; 
  const CFuint nbSpPlus3 = nbSpecies+3;
  const CFuint nbTvs = term->getNbScalarVars(1);
 
  _useForInitialization = getMethodData().isInitializationPhase();
  
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint localStateID = (*statesIdx)[iState];
    State *const currState = states[localStateID];

    if (!_useForInitialization) {
      if (currState->isParUpdatable()) {
	if (!isUpdated[localStateID]) {
	  //	  const CFreal coeff = max(updateCoeff[localStateID]*invcfl, invcfl);
	  const CFreal coeff = updateCoeff[localStateID]*invcfl;
	  
	  // velocity has to remain null
	  rhs(localStateID, nbSpecies, nbEqs) = 0.0;
	  rhs(localStateID, nbSpPlus1, nbEqs) = 0.0;
	  
	  rhs(localStateID, nbSpPlus2, nbEqs) = coeff*(_TWall - (*currState)[nbSpPlus2]);
	  
	  if ( nbTvs > 0 ){
	    rhs(localStateID, nbSpPlus3, nbEqs) = coeff*(_TWall - (*currState)[nbSpPlus3]);
	  }
	  
          const CFuint globalID = idxMapping.getColID(localStateID)*nbEqs;
          // here we have to know how many and which vertices
          // reference the current boundary node to avoid VERY
          // expensive reallocations!!!!
          const CFuint nbEntries = bStatesNeighbors[localStateID].size();
          cf_assert(nbEntries > 0);
	  
          for (CFuint i = 0; i < nbEntries; ++i) {
            const CFuint entryID = bStatesNeighbors[localStateID][i]->getLocalID();
            const CFuint uID  = globalID + nbSpecies;
            const CFuint vID  = globalID + nbSpPlus1;
	    const CFuint eID  = globalID + nbSpPlus2;
	    const CFuint evID = globalID + nbSpPlus3;
	    
	    if (entryID == localStateID) {
	      
	      // condition for the adiabatic wall
	      CFuint gID = globalID;
	      for (CFuint ib = 0; ib < nbEqs; ++ib, ++gID) {
		if (gID < uID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(eID, gID, 0.0);
		  if ( nbTvs > 0 ){
		    jacobMatrix->setValue(evID, gID, 0.0);
		  }
		  continue;
		}
		else if (gID == uID) {
		  jacobMatrix->setValue(uID, gID, coeff);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(eID, gID, 0.0);
		  if ( nbTvs > 0 ){
		    jacobMatrix->setValue(evID, gID, 0.0);
		  }
		  continue;
		}
		else if (gID == vID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, coeff);
		  jacobMatrix->setValue(eID, gID, 0.0);
		  if ( nbTvs > 0 ){
		    jacobMatrix->setValue(evID, gID, 0.0);
		  }
		  continue;
		}
		else if (gID == eID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(eID, gID, coeff);
		  if ( nbTvs > 0 ){
		    jacobMatrix->setValue(evID, gID, 0.0);
		  }
		}
		else if ( (gID == evID) &&  (nbTvs > 0) ) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(eID, gID, 0.0);
		  jacobMatrix->setValue(evID, gID, coeff);
		}
	      }
	    }
	    else {
	      CFuint globalColID =  idxMapping.getColID(entryID)*nbEqs;
	      for (CFuint ib = 0; ib < nbEqs; ++ib, ++globalColID) {
		jacobMatrix->setValue(uID, globalColID, 0.0);
		jacobMatrix->setValue(vID, globalColID, 0.0);
		jacobMatrix->setValue(eID, globalColID, 0.0);
		if ( nbTvs > 0 ){
		  jacobMatrix->setValue(evID, globalColID, 0.0);
		}
	      }
	    }
          }// for (CFuint i = 0; i < nbEntries; ++i) {
          isUpdated[localStateID] = true; // flagging is important!!!!!
        }
      }
    }
    else {
      (*currState)[nbSpecies] = 0.0;
      (*currState)[nbSpPlus1] = 0.0;
      (*currState)[nbSpPlus2] = _TWall;
      if ( nbTvs > 0 ){
        (*currState)[nbSpPlus3] = _TWall;
      }
    }
  }
  
  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();
  
  CFLog(VERBOSE, "StrongNoSlipWallIsothermalTCNEQ2DImpl::computeIsothermal() END\n");
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalTCNEQ2DImpl::computeAdiabatic()
{
  CFLog(VERBOSE, "StrongNoSlipWallIsothermalTCNEQ2DImpl::computeAdiabatic() START\n");
  
  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];
  
  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();
  
  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();
  
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "StrongNoSlipWallAdiabaticNS2DImpl::execute() called for TRS: "
        << trs->getName() << "\n");

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle< std::valarray<State*> > bStatesNeighbors =
    socket_bStatesNeighbors.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();
  const LSSIdxMapping& idxMapping =
    getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();

  _useForInitialization = getMethodData().isInitializationPhase();
  
  SafePtr<MultiScalarTerm<EulerTerm> > term = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<MultiScalarTerm<EulerTerm> >();
  const CFuint nbSpecies = term->getNbScalarVars(0);
  const CFuint nbSpPlus1 = nbSpecies+1;
  
  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint localStateID = (*statesIdx)[iState];
    State *const currState = states[localStateID];
    
    if (!_useForInitialization) {
      if (currState->isParUpdatable()) {
	if (!isUpdated[localStateID]) {
          const CFreal coeff = max(updateCoeff[localStateID]/cfl, 1./cfl);
	  
	  // velocity has to remain null
          rhs(localStateID, nbSpecies, nbEqs) = 0.0;
          rhs(localStateID, nbSpPlus1, nbEqs) = 0.0;
	  
          const CFuint globalID = idxMapping.getColID(localStateID)*nbEqs;
          // here we have to know how many and which vertices
          // reference the current boundary node to avoid VERY
          // expensive reallocations!!!!
          const CFuint nbEntries = bStatesNeighbors[localStateID].size();
          cf_assert(nbEntries > 0);

          for (CFuint i = 0; i < nbEntries; ++i) {
            const CFuint entryID = bStatesNeighbors[localStateID][i]->getLocalID();
            const CFuint uID = globalID + nbSpecies;
            const CFuint vID = globalID + nbSpPlus1;
            if (entryID == localStateID) {
	      CFuint gID = globalID;
	      for (CFuint ib = 0; ib < nbEqs; ++ib, ++gID) {
		if (gID == uID) {
		  jacobMatrix->setValue(uID, gID, coeff);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  continue;
		}
		else if (gID == vID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, coeff);
		  continue;
		}
		else {
		  cf_assert((gID != vID) && (gID != uID));
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, 0.0);
		}
	      }
            }
            else {
	      CFuint globalColID =  idxMapping.getColID(entryID)*nbEqs;
              for (CFuint ib = 0; ib < nbEqs; ++ib, ++globalColID) {
		jacobMatrix->setValue(uID, globalColID, 0.0);
		jacobMatrix->setValue(vID, globalColID, 0.0);
	      }
	    }
          }
          isUpdated[localStateID] = true; // flagging is important!!!!!
        }
      }
    }
    else {
      (*currState)[nbSpecies] = 0.0;
      (*currState)[nbSpPlus1] = 0.0;
    }
  }

  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly(); 
  
  CFLog(VERBOSE, "StrongNoSlipWallIsothermalTCNEQ2DImpl::computeAdiabatic() END\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
