#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "StrongNoSlipWallIsothermalNS3DImpl.hh"
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

MethodCommandProvider<StrongNoSlipWallIsothermalNS3DImpl,
		      FluctuationSplitData,
		      FluctSplitNavierStokesModule>
strongNoSlipWallIsothermalNS3DImplProvider("StrongNoSlipWallIsothermalNS3DImpl");

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalNS3DImpl::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >
    ("TWall","Wall temperature.");
}

//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallIsothermalNS3DImpl::StrongNoSlipWallIsothermalNS3DImpl(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  _useForInitialization(false)
{
   addConfigOptionsTo(this);

   _TWall = 0.0;
   setParameter("TWall",&_TWall);
}

//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallIsothermalNS3DImpl::~StrongNoSlipWallIsothermalNS3DImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongNoSlipWallIsothermalNS3DImpl::needsSockets()
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

void StrongNoSlipWallIsothermalNS3DImpl::setup()
{

  FluctuationSplitCom::setup();

  cf_assert(_TWall > 0.);

  SafePtr<EulerTerm> eulerModel = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<EulerTerm>();

  _TWall /= eulerModel->getTempRef();

  // apply the initial condition
  std::vector< Common::SafePtr<TopologicalRegionSet> >& trsList = getTrsList();

  for(CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
    SafePtr<TopologicalRegionSet> trs = trsList[iTRS];

    CFLogDebugMax("StrongNoSlipWallAdiabaticNS3DImpl::setup() called for TRS: "
		  << trs->getName() << "\n");

    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();

    // (total) internal energy at the wall eWall = cv*Tw
    const CFreal gamma = eulerModel->getGamma();
    const CFreal eWall = eulerModel->getR()/(gamma - 1.)*_TWall;

    for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
      const CFuint localStateID = (*statesIdx)[iState];
      State *const currState = states[localStateID];

      (*currState)[1] = 0.0;
      (*currState)[2] = 0.0;
      (*currState)[3] = 0.0;
      (*currState)[4] = (*currState)[0]*eWall;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalNS3DImpl::executeOnTrs()
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

  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();

  SafePtr<EulerTerm> eulerModel = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<EulerTerm>();

  // (total) internal energy at the wall eWall = cv*Tw
  const CFreal gamma = eulerModel->getGamma();
  const CFreal eWall = eulerModel->getR()/(gamma - 1.)*_TWall;

  _useForInitialization = getMethodData().isInitializationPhase();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint localStateID = (*statesIdx)[iState];
    State *const currState = states[localStateID];

    if (!_useForInitialization) {
      if (currState->isParUpdatable()) {
	if (!isUpdated[localStateID]) {
          // const CFreal coeff = max(updateCoeff[localStateID]/cfl, 1./cfl);
	  const CFreal coeff = updateCoeff[localStateID]/cfl;

          // velocity has to remain null
          rhs(localStateID, 1, nbEqs) = 0.0;
          rhs(localStateID, 2, nbEqs) = 0.0;
	  rhs(localStateID, 3, nbEqs) = 0.0;
	  rhs(localStateID, 4, nbEqs) = 0.0;

          const CFuint globalID = idxMapping.getColID(localStateID)*nbEqs;
          // here we have to know how many and which vertices
          // reference the current boundary node to avoid VERY
          // expensive reallocations!!!!
          const CFuint nbEntries = bStatesNeighbors[localStateID].size();
          cf_assert(nbEntries > 0);

          for (CFuint i = 0; i < nbEntries; ++i) {
            const CFuint entryID =
	      bStatesNeighbors[localStateID][i]->getLocalID();
            const CFuint uID = globalID + 1;
            const CFuint vID = globalID + 2;
	    const CFuint wID = globalID + 3;
	    const CFuint eID = globalID + 4;
	    if (entryID == localStateID) {

	      // condition for the adiabatic wall
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

	      // condition for the isothermal wall
	      // reset the ID to the globalID of the first variable
	      gID = globalID;
	      for (CFuint ib = 0; ib < nbEqs; ++ib, ++gID) {
		if (gID == globalID) {
		  // the factor "2" depends on the fact that
		  // the diagonal entry in position eID will receive
		  // an additional contribution (+coeff) from the time
		  // RHS command:
		  // -Ew*2*coeff*dRho + 2*coeff*dRhoE = 0 =>
		  // -Ew*dRho + dRhoE = 0, as we want !!

		  // @TODO this may not work if the time integrator is
		  // != Backward Euler
		  jacobMatrix->setValue(eID, gID, -2.0*coeff*eWall);
		  continue;
		}
		else if (gID == eID) {
		  jacobMatrix->setValue(eID, gID, coeff);
		  continue;
		}
		else {
		  cf_assert((gID != globalID) && (gID != eID));
		  jacobMatrix->setValue(eID, gID, 0.0);
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
      (*currState)[1] = 0.0;
      (*currState)[2] = 0.0;
      (*currState)[3] = 0.0;
      (*currState)[4] = (*currState)[0]*eWall;
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
