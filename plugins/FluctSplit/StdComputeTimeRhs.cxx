#include "FluctSplit/FluctSplit.hh"
#include "StdComputeTimeRhs.hh"
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

MethodCommandProvider<StdComputeTimeRhs, FluctuationSplitData, FluctSplitModule> 
rdsStdComputeTimeRhs("StdTimeRhs");

//////////////////////////////////////////////////////////////////////////////

StdComputeTimeRhs::StdComputeTimeRhs(const std::string& name) :
  FluctuationSplitCom(name),
  _lss(CFNULL),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_discardTimeJacob("discardTimeJacob")
{
}

//////////////////////////////////////////////////////////////////////////////

StdComputeTimeRhs::~StdComputeTimeRhs()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdComputeTimeRhs::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_discardTimeJacob);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdComputeTimeRhs::setup()
{
  FluctuationSplitCom::setup();

  // linear system solver
  _lss = getMethodData().getLinearSystemSolver()[0];

}

//////////////////////////////////////////////////////////////////////////////

void StdComputeTimeRhs::execute()
{
  CFAUTOTRACE;
  
  if (getMethodData().doComputeJacobian()) {
    SafePtr<LSSMatrix> jacobMatrix = _lss->getMatrix();
    
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
    
    const CFuint nbStates = states.size();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    
    const CFreal cfl = getMethodData().getCFL()->getCFLValue();
    const CFreal overCfl = 1.0 / cfl;
    
    const LSSIdxMapping& idxMapping =
      getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();
    
    // add the diagonal entries in the jacobian (updateCoeff/CFL)
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      if (states[iState]->isParUpdatable()) {
	
	const CFreal diagValue = updateCoeff[iState]*overCfl;
	CFuint globalID = idxMapping.getColID
	  (states[iState]->getLocalID())*nbEqs;
	
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq, ++globalID) {
	  jacobMatrix->addValue(globalID, globalID, diagValue);
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
