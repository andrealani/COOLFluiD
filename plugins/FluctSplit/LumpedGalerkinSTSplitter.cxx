#include "LumpedGalerkinSTSplitter.hh"
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
      
MethodStrategyProvider<LumpedGalerkinSTSplitter,
		       FluctuationSplitData,
		       SourceTermSplitter,
                       FluctSplitModule>
lumpedGalerkinSplitterProvider("LumpedGalerkin");

//////////////////////////////////////////////////////////////////////////////

LumpedGalerkinSTSplitter::LumpedGalerkinSTSplitter(const std::string& name) :
  SourceTermSplitter(name)
{
}

//////////////////////////////////////////////////////////////////////////////

LumpedGalerkinSTSplitter::~LumpedGalerkinSTSplitter()
{
}

//////////////////////////////////////////////////////////////////////////////

void LumpedGalerkinSTSplitter::distribute(std::vector<RealVector>& residual)
{
  DistributionData& ddata = getMethodData().getDistributionData();
  SafePtr<vector<State*> > states = getMethodData().getDistributionData().states;
  const CFreal coeff = 1./(1. + PhysicalModelStack::getActive()->getDim());
  for (CFuint iState = 0; iState < states->size(); ++iState) {
    residual[iState] -= ddata.nodalST[iState]*coeff;
 }
}

//////////////////////////////////////////////////////////////////////////////

void LumpedGalerkinSTSplitter::computeSourceTerm(const InwardNormalsData& normalsData)
{
  DistributionData& ddata = getMethodData().getDistributionData();
  SafePtr<ComputeSourceTermFSM> st = 
    getMethodData().getSourceTermComputer(ddata.sourceTermID);
  st->setPerturb(true); // AL: is this needed?
  st->computeSourceFSM(ddata.cell, ddata.nodalST, normalsData);
  st->setPerturb(false); // AL: is this needed?

}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
