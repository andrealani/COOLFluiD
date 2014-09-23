#include "BetaSTSplitter.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////
      
MethodStrategyProvider<BetaSTSplitter,
		       FluctuationSplitData,
		       SourceTermSplitter,
                       FluctSplitModule>
betaSTSplitterProvider("Beta");
      
MethodStrategyProvider<BetaSTSplitter,
		       FluctuationSplitData,
		       SourceTermSplitter,
                       FluctSplitModule>
beta1STSplitterProvider("Beta1");
    
MethodStrategyProvider<BetaSTSplitter,
		       FluctuationSplitData,
		       SourceTermSplitter,
                       FluctSplitModule>
beta2STSplitterProvider("Beta2");

//////////////////////////////////////////////////////////////////////////////

BetaSTSplitter::BetaSTSplitter(const std::string& name) :
  SourceTermSplitter(name),
  _tmpSource()
{
  addConfigOptionsTo(this);
  _excludeBStates = false;
  setParameter("ExcludeBStates",&_excludeBStates);
}
    
//////////////////////////////////////////////////////////////////////////////

BetaSTSplitter::~BetaSTSplitter()
{
}

//////////////////////////////////////////////////////////////////////////////

void BetaSTSplitter::distribute(std::vector<RealVector>& residual)
{
  DistributionData& ddata = getMethodData().getDistributionData();
  getMethodData().getSourceTermComputer(ddata.sourceTermID)->setPerturb(true);
  if (_excludeBStates) {
    distributeExceptBStates(residual);
  }
  else {
    const CFuint nbStates = residual.size();
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      residual[iState] -= (*ddata.currBetaMat)[iState]*ddata.phiS;
    }
  }
  getMethodData().getSourceTermComputer(ddata.sourceTermID)->setPerturb(false);
}
    
//////////////////////////////////////////////////////////////////////////////

void BetaSTSplitter::distributeExceptBStates(std::vector<RealVector>& residual)
{
  // AL: this was WRONG !!! only points on the symmetry axis should not 
  //     receive residual contribution
  // DistributionData& ddata = getMethodData().getDistributionData();
  //   DataHandle<bool> isBState = socket_isBState.getDataHandle();
  //   const CFuint nstates = ddata.states->size();
  //   for (CFuint iState = 0; iState < nstates; ++iState) {
  //     if (!isBState[(*ddata.states)[iState]->getLocalID()]){
  //       residual[iState] -= (*ddata.currBetaMat)[iState]*ddata.phiS;
  //     }
  //   }
  
  // only states with R=0 (y=0) should not get the residual
  DistributionData& ddata = getMethodData().getDistributionData();
  DataHandle<bool> isBState = socket_isBState.getDataHandle();
  const CFuint nstates = ddata.states->size();
  for (CFuint iState = 0; iState < nstates; ++iState) {
    const CFreal radius = (*ddata.states)[iState]->getCoordinates()[YY];
    if (radius > 1e-10) {
      cf_assert ((!isBState[(*ddata.states)[iState]->getLocalID()]) || radius > 1e-10);
      residual[iState] -= (*ddata.currBetaMat)[iState]*ddata.phiS;
    }
  }  
}

//////////////////////////////////////////////////////////////////////////////

void BetaSTSplitter::computeSourceTerm(const InwardNormalsData& normalsData)
{
  DistributionData& ddata = getMethodData().getDistributionData();
  getMethodData().getSourceTermComputer(ddata.sourceTermID)->
    computeSourceFSM(ddata.cell,_tmpSource,normalsData);
  
  // transform the residual to distribution variables because the betas
  // are defined in distribution variables
  ddata.phiS = *getMethodData().getSolutionToDistribMatTrans()->
    transformFromRef(&_tmpSource);
}
      
//////////////////////////////////////////////////////////////////////////////

void BetaSTSplitter::setup()
{
  SourceTermSplitter::setup();
  getMethodData().getDistributionData().computeBetas = true;
  _tmpSource.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

void BetaSTSplitter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("ExcludeBStates",
				   "Exclude distribution of the source term to boundary states.");
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
