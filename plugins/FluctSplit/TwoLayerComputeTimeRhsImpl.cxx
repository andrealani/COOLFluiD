#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "TwoLayerComputeTimeRhsImpl.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/LSSIdxMapping.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TwoLayerComputeTimeRhsImpl, FluctuationSplitData, FluctSplitSpaceTimeModule> rdsTwoLayerComputeTimeRhsImpl("TwoLayerTimeRhsImpl");

//////////////////////////////////////////////////////////////////////////////

TwoLayerComputeTimeRhsImpl::TwoLayerComputeTimeRhsImpl(const std::string& name) :
  FluctuationSplitCom(name),
  _lss(CFNULL),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_interUpdateCoeff("interUpdateCoeff")
{
}

//////////////////////////////////////////////////////////////////////////////

TwoLayerComputeTimeRhsImpl::~TwoLayerComputeTimeRhsImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TwoLayerComputeTimeRhsImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_interUpdateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerComputeTimeRhsImpl::setup()
{
  FluctuationSplitCom::setup();

  // linear system solver
  _lss = getMethodData().getLinearSystemSolver()[0];

}

//////////////////////////////////////////////////////////////////////////////

void TwoLayerComputeTimeRhsImpl::execute()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle< CFreal> interUpdateCoeff = socket_interUpdateCoeff.getDataHandle();

  SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  const CFreal cfl = getMethodData().getCFL()->getCFLValue();
  const CFreal overCfl = 1.0 / cfl;


  const LSSIdxMapping& idxMapping =
    getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();

  // add the diagonal entries in the jacobian (updateCoeff/CFL)
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    if (states[iState]->isParUpdatable()) {

      CFreal diagValue1 = interUpdateCoeff[iState]*overCfl;
      CFreal diagValue2 = updateCoeff[iState]*overCfl;

      CFuint globalID = idxMapping.getColID(states[iState]->getLocalID())*2*nbEqs;
      for (CFuint iEq = 0; iEq < 2*nbEqs; ++iEq, ++globalID) {
        if (iEq < nbEqs) jacobMatrix->addValue(globalID, globalID, diagValue1);
        else jacobMatrix->addValue(globalID, globalID, diagValue2);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
