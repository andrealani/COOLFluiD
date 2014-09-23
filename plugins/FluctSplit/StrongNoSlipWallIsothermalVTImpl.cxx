#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/StrongNoSlipWallIsothermalVTImpl.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFL.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongNoSlipWallIsothermalVTImpl,
		      FluctuationSplitData,
		      FluctSplitModule>
strongNoSlipWallIsothermalVTImplProvider("StrongNoSlipWallIsothermalVTImpl");

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalVTImpl::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >
    ("TWall","Wall temperature.");
}

//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallIsothermalVTImpl::StrongNoSlipWallIsothermalVTImpl(const std::string& name) :
  StrongImplBC(name)
{
   addConfigOptionsTo(this);

   _TWall = 0.0;
   setParameter("TWall",&_TWall);
}

//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallIsothermalVTImpl::~StrongNoSlipWallIsothermalVTImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalVTImpl::setup()
{
  CFAUTOTRACE;
  
  StrongImplBC::setup();
  
  cf_assert(_TWall > 0.);
  
  SafePtr<EulerTerm> term = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<EulerTerm>();
  _TWall /= term->getTempRef();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbSpecies = 1;
  const CFuint nbSpPlus1 = nbSpecies+1; 
  const CFuint nbSpPlus2 = nbSpecies+2; 
  const CFuint nbSpPlus3 = nbSpecies+3; 
  
  // apply always the initial condition
  std::vector< Common::SafePtr<TopologicalRegionSet> >& trsList = getTrsList();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< vector<bool> > discardTimeJacob =
    socket_discardTimeJacob.getDataHandle();
  
  for(CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> trs = trsList[iTRS];

    CFLogDebugMax("StrongNoSlipWallIsothermalVTImpl::setup() called for TRS: "
		  << trs->getName() << "\n");

    Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();
        
    for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
      const CFuint localStateID = (*statesIdx)[iState];
      State *const currState = states[localStateID];
      // fix the wall temperature (assume rho_i v T Tv variables)
      
      (*currState)[nbSpecies] = 0.0;
      (*currState)[nbSpPlus1] = 0.0;

      if (dim == DIM_3D) {
	(*currState)[nbSpPlus2] = 0.0;
	(*currState)[nbSpPlus3] = _TWall;
      }
      else {
	(*currState)[nbSpPlus2] = _TWall;
      }
      
      /// @TODO AL: here check that the order in which BCs are setup
      /// is the order of their execution !!!
      discardTimeJacob[localStateID].resize(nbEqs);
      for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
	if (dim == DIM_3D) {
	  discardTimeJacob[localStateID][iVar] = 
	    (iVar == nbSpecies || iVar == nbSpPlus1 || iVar == nbSpPlus2 || iVar == nbSpPlus3) ? 
	    true : false;
	}
	else {
	  discardTimeJacob[localStateID][iVar] = 
	    (iVar == nbSpecies || iVar == nbSpPlus1 || iVar == nbSpPlus2) ? 
	    true : false;
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalVTImpl::executeOnTrs()
{  
  (PhysicalModelStack::getActive()->getDim() == DIM_3D) ? compute3D() : compute2D();
}
//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalVTImpl::compute2D()
{ 
  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "StrongNoSlipWallIsothermalVTImpl::execute() called for TRS: "
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
  const CFuint nbSpecies = 1;
  const CFuint nbSpPlus1 = nbSpecies+1; 
  const CFuint nbSpPlus2 = nbSpecies+2; 
  
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
	  
          const CFuint globalID = idxMapping.getColID(localStateID)*nbEqs;
          // here we have to know how many and which vertices
          // reference the current boundary node to avoid VERY
          // expensive reallocations!!!!
          const CFuint nbEntries = bStatesNeighbors[localStateID].size();
          cf_assert(nbEntries > 0);
	  
          for (CFuint i = 0; i < nbEntries; ++i) {
            const CFuint entryID = 
	      bStatesNeighbors[localStateID][i]->getLocalID();
            const CFuint uID  = globalID + nbSpecies;
            const CFuint vID  = globalID + nbSpPlus1;
	    const CFuint eID  = globalID + nbSpPlus2;
	    
	    if (entryID == localStateID) {
	      
	      // condition for the adiabatic wall
	      CFuint gID = globalID;
	      for (CFuint ib = 0; ib < nbEqs; ++ib, ++gID) {
		if (gID < uID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(eID, gID, 0.0);
		  continue;
		}
		else if (gID == uID) {
		  jacobMatrix->setValue(uID, gID, coeff);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(eID, gID, 0.0);
		  continue;
		}
		else if (gID == vID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, coeff);
		  jacobMatrix->setValue(eID, gID, 0.0);
		  continue;
		}
		else if (gID == eID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(eID, gID, coeff);
		}
	      }
	    }
            else {
	      CFuint globalColID =  idxMapping.getColID(entryID)*nbEqs;
              for (CFuint ib = 0; ib < nbEqs; ++ib, ++globalColID) {
		jacobMatrix->setValue(uID, globalColID, 0.0);
		jacobMatrix->setValue(vID, globalColID, 0.0);
		jacobMatrix->setValue(eID, globalColID, 0.0);
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
      (*currState)[nbSpPlus2] = _TWall;
    }
  }
  
  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalVTImpl::compute3D()
{ 
  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "StrongNoSlipWallIsothermalVTImpl::execute() called for TRS: "
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
  const CFuint nbSpecies = 1;
  const CFuint nbSpPlus1 = nbSpecies+1; 
  const CFuint nbSpPlus2 = nbSpecies+2; 
  const CFuint nbSpPlus3 = nbSpecies+3; 
  
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
	  rhs(localStateID, nbSpPlus2, nbEqs) = 0.0;
	  rhs(localStateID, nbSpPlus3, nbEqs) = coeff*(_TWall - (*currState)[nbSpPlus3]);
	  
          const CFuint globalID = idxMapping.getColID(localStateID)*nbEqs;
          // here we have to know how many and which vertices
          // reference the current boundary node to avoid VERY
          // expensive reallocations!!!!
          const CFuint nbEntries = bStatesNeighbors[localStateID].size();
          cf_assert(nbEntries > 0);
	  
          for (CFuint i = 0; i < nbEntries; ++i) {
            const CFuint entryID = 
	      bStatesNeighbors[localStateID][i]->getLocalID();
            const CFuint uID  = globalID + nbSpecies;
            const CFuint vID  = globalID + nbSpPlus1;
	    const CFuint wID  = globalID + nbSpPlus2;
	    const CFuint eID  = globalID + nbSpPlus3;
	    
	    if (entryID == localStateID) {
	      
	      // condition for the adiabatic wall
	      CFuint gID = globalID;
	      for (CFuint ib = 0; ib < nbEqs; ++ib, ++gID) {
		if (gID < uID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(wID, gID, 0.0);
		  jacobMatrix->setValue(eID, gID, 0.0);
		  continue;
		}
		else if (gID == uID) {
		  jacobMatrix->setValue(uID, gID, coeff);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(wID, gID, 0.0);
		  jacobMatrix->setValue(eID, gID, 0.0);
		  continue;
		}
		else if (gID == vID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, coeff);
		  jacobMatrix->setValue(wID, gID, 0.0);
		  jacobMatrix->setValue(eID, gID, 0.0);
		  continue;
		}
		else if (gID == wID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(wID, gID, coeff);
		  jacobMatrix->setValue(eID, gID, 0.0);
		  continue;
		}
		else if (gID == eID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(wID, gID, 0.0);
		  jacobMatrix->setValue(eID, gID, coeff);
		}
	      }
	    }
            else {
	      CFuint globalColID =  idxMapping.getColID(entryID)*nbEqs;
              for (CFuint ib = 0; ib < nbEqs; ++ib, ++globalColID) {
		jacobMatrix->setValue(uID, globalColID, 0.0);
		jacobMatrix->setValue(vID, globalColID, 0.0);
		jacobMatrix->setValue(wID, globalColID, 0.0);
		jacobMatrix->setValue(eID, globalColID, 0.0);
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
      (*currState)[nbSpPlus2] = 0.0;
      (*currState)[nbSpPlus3] = _TWall;
    }
  }
  
  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
