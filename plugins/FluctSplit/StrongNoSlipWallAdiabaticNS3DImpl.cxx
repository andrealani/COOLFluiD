#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/StrongNoSlipWallAdiabaticNS3DImpl.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFL.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongNoSlipWallAdiabaticNS3DImpl, FluctuationSplitData, FluctSplitModule> strongNoSlipWallAdiabaticNS3DImplProvider("StrongNoSlipWallAdiabaticNS3DImpl");

//////////////////////////////////////////////////////////////////////////////
void StrongNoSlipWallAdiabaticNS3DImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("u0","x-coponent of the velocity of the wall");
   options.addConfigOption< CFreal >("v0","y-coponent of the velocity of the wall");
options.addConfigOption< CFreal >("w0","z-coponent of the velocity of the wall");

}


//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallAdiabaticNS3DImpl::StrongNoSlipWallAdiabaticNS3DImpl(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  _useForInitialization(false)
{

 addConfigOptionsTo(this);

  m_u0 = 0.0;
  setParameter("u0",&m_u0);
  m_v0 = 0.0;
  setParameter("v0",&m_v0);
  m_w0 = 0.0;
  setParameter("w0", &m_w0);

}

//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallAdiabaticNS3DImpl::~StrongNoSlipWallAdiabaticNS3DImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongNoSlipWallAdiabaticNS3DImpl::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_bStatesNeighbors);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallAdiabaticNS3DImpl::executeOnTrs()
{
  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "StrongNoSlipWallAdiabaticNS3DImpl::execute() called for TRS: "
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

  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint localStateID = (*statesIdx)[iState];
    State *const currState = states[localStateID];

    if (!_useForInitialization) {
      if (currState->isParUpdatable()) {
	if (!isUpdated[localStateID]) {
          const CFreal coeff = max(updateCoeff[localStateID]/cfl, 1./cfl);

          // velocity has to remain null
          rhs(localStateID, 1, nbEqs) = 0.0;
          rhs(localStateID, 2, nbEqs) = 0.0;
	  rhs(localStateID, 3, nbEqs) = 0.0;

	(*currState)[1] = m_u0;
        (*currState)[2] = m_v0;
	(*currState)[3] = m_w0;



          const CFuint globalID = idxMapping.getColID(localStateID)*nbEqs;
          // here we have to know how many and which vertices
          // reference the current boundary node to avoid VERY
          // expensive reallocations!!!!
          const CFuint nbEntries = bStatesNeighbors[localStateID].size();
          cf_assert(nbEntries > 0);

          for (CFuint i = 0; i < nbEntries; ++i) {
            const CFuint entryID = bStatesNeighbors[localStateID][i]->getLocalID();
            const CFuint uID = globalID + 1;
            const CFuint vID = globalID + 2;
	    const CFuint wID = globalID + 3;

            if (entryID == localStateID) {
	      CFuint gID = globalID;
	      for (CFuint ib = 0; ib < nbEqs; ++ib, ++gID) {
		if (gID == uID) {
		  jacobMatrix->setValue(uID, gID, coeff);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(wID, gID, 0.0);
		  continue;
		}
		else if (gID == vID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, coeff);
		  jacobMatrix->setValue(wID, gID, 0.0);
		  continue;
		}
		else if (gID == wID) {
		  jacobMatrix->setValue(uID, gID, 0.0);
                  jacobMatrix->setValue(vID, gID, 0.0);
                  jacobMatrix->setValue(wID, gID, coeff);
                  continue;
		}
		else {
		  cf_assert((gID != vID) && (gID != uID) && (gID != wID));
		  jacobMatrix->setValue(uID, gID, 0.0);
		  jacobMatrix->setValue(vID, gID, 0.0);
		  jacobMatrix->setValue(wID, gID, 0.0);
		}
	      }
            }
            else {
	      CFuint globalColID =  idxMapping.getColID(entryID)*nbEqs;
              for (CFuint ib = 0; ib < nbEqs; ++ib, ++globalColID) {
		jacobMatrix->setValue(uID, globalColID, 0.0);
		jacobMatrix->setValue(vID, globalColID, 0.0);
		jacobMatrix->setValue(wID, globalColID, 0.0);
	      }
	    }
          }
          isUpdated[localStateID] = true; // flagging is important!!!!!
        }
      }
    }
    else {

      (*currState)[1] = m_u0;
      (*currState)[2] = m_v0;
      (*currState)[3] = m_w0;
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
