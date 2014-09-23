#include "FluctSplit/FluctSplit.hh"
#include "SuperOutletImplMHD2DProjection.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperOutletImplMHD2DProjection,
		      FluctuationSplitData,
		      FluctSplitModule>
superOutletImplMHD2DProjectionProvider("SuperOutletImplMHD2DProjection");

//////////////////////////////////////////////////////////////////////////////

void SuperOutletImplMHD2DProjection::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >
     ("DiagCoeffFactor","Factor to control the diagonal coefficient.");

  options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the supersonic outlet.");
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletImplMHD2DProjection::SuperOutletImplMHD2DProjection(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_isUpdated("isUpdated"),
  socket_states("states"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  socket_updateCoeff("updateCoeff")
{
   addConfigOptionsTo(this);

   _diagCoeffFactor = 1.0;
   setParameter("DiagCoeffFactor",&_diagCoeffFactor);

  _refPhi = 0.;
   setParameter("refPhi",&_refPhi);
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletImplMHD2DProjection::~SuperOutletImplMHD2DProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
SuperOutletImplMHD2DProjection::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_states);
  result.push_back(&socket_bStatesNeighbors);
  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletImplMHD2DProjection::executeOnTrs()
{
  CFAUTOTRACE;
  SafePtr<LinearSystemSolver> lss =
    getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "SuperOutletImplMHD2DProjection::execute() called for TRS: "
  << trs->getName() << "\n");

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  DataHandle< std::valarray<Framework::State*> > bStatesNeighbors =
    socket_bStatesNeighbors.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();
  const CFreal invcfl = 1.0 / cfl;

  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();

  const LSSIdxMapping& idxMapping =
    getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint localStateID = (*statesIdx)[iState];
    State *const currState = states[localStateID];

    if (currState->isParUpdatable()) {

      if (!isUpdated[localStateID]) {
	const CFreal kcoeff = updateCoeff[localStateID] * invcfl;
	const CFreal coeff  = max(kcoeff, _diagCoeffFactor * invcfl);

	// for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	// must include inletState so it will work for unsteady BCs
	// which are updated every iteration
	rhs(localStateID, 8, nbEqs) =
	  (coeff + kcoeff) * (_refPhi - (*currState)[8]);

        // here we have to know how many and which vertices
        // reference the current boundary node to avoid VERY
        // expensive reallocations!!!!
        const CFuint nbEntries = bStatesNeighbors[localStateID].size();
        cf_assert(nbEntries > 0);
	CFuint globalRowID = idxMapping.getColID(localStateID)*nbEqs;
	const CFuint globalID8 = globalRowID + 8;

	for (CFuint i = 0; i < nbEntries; ++i) {
	  const CFuint entryID = bStatesNeighbors[localStateID][i]->getLocalID();
	  CFuint globalID = idxMapping.getColID(entryID)*nbEqs;
	  if (entryID == localStateID) {
	    // reset to 0 all the entries in subblock row 8
	    // except for the diagonal value
	    for (CFuint iEq = 0; iEq < nbEqs; ++iEq, ++globalID) {
	      const CFreal dCoeff = (globalID == globalID8) ? coeff : 0.0;
	      jacobMatrix->setValue(globalID8, globalID, dCoeff);
	    }
	  }
	  else {
	    for (CFuint iEq = 0; iEq < nbEqs; ++iEq, ++globalID) {
	      jacobMatrix->setValue(globalID8, globalID, 0.0);
	    }
	  }
	}

	isUpdated[localStateID] = true; // flagging is important!!!!!
      }
    }
  }

  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletImplMHD2DProjection::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace FluctSplit



} // namespace COOLFluiD
