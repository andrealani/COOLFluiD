#include "FluctSplit/FluctSplit.hh"
#include "StrongSymmetryPlaneImpl.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFL.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongSymmetryPlaneImpl,
		      FluctuationSplitData,
		      FluctSplitModule>
strongSymmetryPlaneImplProvider("StrongSymmetryPlaneImpl");

//////////////////////////////////////////////////////////////////////////////

void StrongSymmetryPlaneImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< vector<CFuint> >
     ("annullVarID","Tells which variables should be annulled.");
}

//////////////////////////////////////////////////////////////////////////////

StrongSymmetryPlaneImpl::StrongSymmetryPlaneImpl(const std::string& name) :
  StrongImplBC(name)
{
  addConfigOptionsTo(this);

  _varIDToAnnull = vector<CFuint>();
  setParameter("annullVarID",&_varIDToAnnull);
}

//////////////////////////////////////////////////////////////////////////////

StrongSymmetryPlaneImpl::~StrongSymmetryPlaneImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void StrongSymmetryPlaneImpl::setup()
{
  StrongImplBC::setup();
  
  // fix for restart
  // apply the initial condition
  vector<Common::SafePtr<TopologicalRegionSet> >& trsList = getTrsList();
  const CFuint nbAnnullVars = _varIDToAnnull.size();
  cf_assert(nbAnnullVars > 0);
  
  DataHandle< vector<bool> > discardTimeJacob =
    socket_discardTimeJacob.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for(CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> trs = trsList[iTRS];

    CFLogDebugMax("StrongSymmetryPlaneImpl::setup() called for TRS: "
		  << trs->getName() << "\n");

    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();
   
    for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
      const CFuint localStateID = (*statesIdx)[iState];
      State *const currState = states[localStateID];
      
      /// @TODO AL: here check that the order in which BCs are setup
      /// is the order of their execution !!!
      discardTimeJacob[localStateID].resize(nbEqs);
      for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
	discardTimeJacob[localStateID][iVar] = false;
      } 
      for (CFuint iv = 0; iv < nbAnnullVars; ++iv) {
	(*currState)[_varIDToAnnull[iv]] = 0.0;
	discardTimeJacob[localStateID][iv] = true;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSymmetryPlaneImpl::executeOnTrs()
{
  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  CFLogDebugMax
    ( "StrongNoSlipWallAdiabaticNS2DImpl::execute() called for TRS: "
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

  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();
  const CFuint nbAnnullVars = _varIDToAnnull.size();

  _useForInitialization = getMethodData().isInitializationPhase();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint localStateID = (*statesIdx)[iState];
    State *const currState = states[localStateID];

    if (!_useForInitialization) {
      if (currState->isParUpdatable()) {
	if (!isUpdated[localStateID]) {
          const CFreal coeff = max(updateCoeff[localStateID]/cfl, 1./cfl);

	  // annull the RHS components
	  for (CFuint iv = 0; iv < nbAnnullVars; ++iv) {
	    rhs(localStateID, _varIDToAnnull[iv], nbEqs) = 0.0;
	  }

          const CFuint globalID = idxMapping.getColID(localStateID)*nbEqs;
          // here we have to know how many and which vertices
          // reference the current boundary node to avoid VERY
          // expensive reallocations!!!!
          const CFuint nbEntries = bStatesNeighbors[localStateID].size();
          cf_assert(nbEntries > 0);

          for (CFuint i = 0; i < nbEntries; ++i) {
            const CFuint entryID =
	      bStatesNeighbors[localStateID][i]->getLocalID();

            if (entryID == localStateID) {
	      CFuint gID = globalID;
	      for (CFuint ib = 0; ib < nbEqs; ++ib, ++gID) {
		for (CFuint iv = 0; iv < nbAnnullVars; ++iv) {
		  const CFuint uID = globalID + _varIDToAnnull[iv];
		  if (gID == uID) {
		    jacobMatrix->setValue(uID, gID, coeff);
		  }
		  else {
		    jacobMatrix->setValue(uID, gID, 0.0);
		  }
		}
	      }
            }
            else {
	      CFuint globalColID =  idxMapping.getColID(entryID)*nbEqs;
              for (CFuint ib = 0; ib < nbEqs; ++ib, ++globalColID) {
		for (CFuint iv = 0; iv < nbAnnullVars; ++iv) {
		  const CFuint uID = globalID + _varIDToAnnull[iv];
		  jacobMatrix->setValue(uID, globalColID, 0.0);
		}
	      }
	    }
          }
          isUpdated[localStateID] = true; // flagging is important!!!!!
        }
      }
    }
    else {
      for (CFuint iv = 0; iv < nbAnnullVars; ++iv) {
	(*currState)[_varIDToAnnull[iv]] = 0.0;
      }
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
