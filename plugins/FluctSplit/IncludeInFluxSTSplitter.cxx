#include "IncludeInFluxSTSplitter.hh"
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
      
MethodStrategyProvider<IncludeInFluxSTSplitter,
		       FluctuationSplitData,
		       SourceTermSplitter,
                       FluctSplitModule>
includeInFluxSTSplitterProvider("IncludeInFlux");
      
//////////////////////////////////////////////////////////////////////////////

IncludeInFluxSTSplitter::IncludeInFluxSTSplitter(const std::string& name) :
  SourceTermSplitter(name),
  _tmpSource()
{
}
    
//////////////////////////////////////////////////////////////////////////////

IncludeInFluxSTSplitter::~IncludeInFluxSTSplitter()
{
}

//////////////////////////////////////////////////////////////////////////////

void IncludeInFluxSTSplitter::distribute(std::vector<RealVector>& residual)
{
}
    
//////////////////////////////////////////////////////////////////////////////

void IncludeInFluxSTSplitter::computeSourceTerm(const InwardNormalsData& normalsData)
{
  DistributionData& ddata = getMethodData().getDistributionData();
  getMethodData().getSourceTermComputer(ddata.sourceTermID)->
    computeSourceFSM(ddata.cell,_tmpSource,normalsData);
  
  // transform the residual to distribution variables because the betas
  // are defined in distribution variables
  ddata.phiS = _tmpSource;
}
    
//////////////////////////////////////////////////////////////////////////////

void IncludeInFluxSTSplitter::setup()
{
  SourceTermSplitter::setup();
  getMethodData().getDistributionData().computeBetas = true;
  _tmpSource.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
