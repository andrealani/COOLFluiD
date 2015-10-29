#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "TwoLayerSuperInletImpl.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerSuperInletImpl, FluctuationSplitData, FluctSplitSpaceTimeModule> TwoLayerSuperInletImplProvider("TwoLayerSuperInletImpl");

//////////////////////////////////////////////////////////////////////////////

TwoLayerSuperInletImpl::TwoLayerSuperInletImpl(const std::string& name) :
  SuperInletImpl(name),
  socket_interRhs("interRhs"),
  socket_interStates("interStates")
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerSuperInletImpl::~TwoLayerSuperInletImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerSuperInletImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = SuperInletImpl::needsSockets();

  result.push_back(&socket_interRhs);
  result.push_back(&socket_interStates);

  return result;
}


//////////////////////////////////////////////////////////////////////////////

void TwoLayerSuperInletImpl::executeOnTrs()
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
  CFLogDebugMax( "SuperInletImpl::execute() called for TRS: "
  << trs->getName() << "\n");

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> interRhs = socket_interRhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<State*> interStates = socket_interStates.getDataHandle();
  DataHandle< std::valarray<Framework::State*> > bStatesNeighbors = socket_bStatesNeighbors.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  RealVector variables(dim+1);

  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();

  // block accumulator 1*1
  auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(1, 1, nbEqs*2));

  const CFreal invcfl = 1.0 / getMethodData().getCFL()->getCFLValue();

  RealVector inletState(nbEqs);

  bool applyBC = true;
  State dimState;
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    State *const state = states[stateID];
    State *const intState = interStates[stateID];

    if (state->isParUpdatable()) {

      // find if we should apply the BC to this node
      if (m_checkCondition) {
        const CFreal applyBCvalue = m_condition.Eval(&state->getCoordinates()[0]);
        applyBC = ((!isUpdated[stateID]) && (applyBCvalue > 0.0));
      }
      else {
        applyBC = (!isUpdated[stateID]);
      }

      // apply the BC to this node
      if (applyBC){

        //Set the values of the variables xyz + time
        for (CFuint i = 0; i < state->getCoordinates().size();++i){
          variables[i] = state->getCoordinates()[i];
        }
        variables[state->getCoordinates().size()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

        //Evaluate the function
        m_vFunction.evaluate(variables,*_input);

        //Set the state value
        if (_inputAdimensionalValues){
          inletState = *_inputToUpdateVar->transform(_input);
        }
        else{
          dimState = *_inputToUpdateVar->transform(_input);
          _varSet->setAdimensionalValues(dimState, inletState);
        }

        // compute coefficients
	const CFreal kcoeff = updateCoeff[stateID]*invcfl;
	const CFreal coeff  = max(kcoeff, _diagCoeffFactor*invcfl);

        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          // must include inletState so it will work for unsteady BCs
          // which are updated every iteration
          rhs(stateID, iEq, nbEqs) = (inletState[iEq] - (*state)[iEq]);
          interRhs(stateID, iEq, nbEqs) = (inletState[iEq] - (*intState)[iEq]);
        }

        acc->setRowIndex(0, stateID);

        // here we have to know how many and which vertices
        // reference the current boundary node to avoid VERY
        // expensive reallocations!!!!
        const CFuint nbEntries = bStatesNeighbors[stateID].size();
        cf_assert(nbEntries > 0);

        for (CFuint i = 0; i < nbEntries; ++i) {
          const CFuint entryID = bStatesNeighbors[stateID][i]->getLocalID();
          acc->setColIndex(0, entryID);
          acc->setValue(0.0);
          if (entryID == stateID) {
            for (CFuint ib = 0; ib < nbEqs; ++ib) {
              acc->setValue(0,0,ib,ib, coeff);
            }
          }
          jacobMatrix->setValues(*acc);
        }
        isUpdated[stateID] = true; // flagging is important!!!!!
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
