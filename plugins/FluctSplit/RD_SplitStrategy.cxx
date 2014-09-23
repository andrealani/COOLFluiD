#include "FluctSplit/RD_SplitStrategy.hh"

#include "Framework/GeometricEntity.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RD_SplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitModule>
rdFluctSplitStrategyProvider("RD");

//////////////////////////////////////////////////////////////////////////////

RD_SplitStrategy::RD_SplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  m_splitter(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

RD_SplitStrategy::~RD_SplitStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void RD_SplitStrategy::computeFluctuation(vector<RealVector>& residual)
{
  setCurrentCell();
  
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DistributionData& ddata = getMethodData().getDistributionData();
  
  // compute the residual and the upwind parameters k in this cell
  m_splitter->computeK(*ddata.states, normals[ddata.cellID]);

  
  const CFuint nbST = getMethodData().getSourceTermSplitter()->size();
  if (nbST == 1) {
    ddata.sourceTermID = 0;
    getMethodData().getSourceTermSplitter(0)->computeSourceTerm(*normals[ddata.cellID]);
  
    m_splitter->distribute(residual);

    getMethodData().getSourceTermSplitter(0)->distribute(residual);

  }
  else {
    m_splitter->distribute(residual);
    
    for (CFuint i = 0; i < nbST; ++i) {
      ddata.sourceTermID = i;
      getMethodData().getSourceTermSplitter(i)->computeSourceTerm(*normals[ddata.cellID]); 
      getMethodData().getSourceTermSplitter(i)->distribute(residual);
    }    
  }

  
}

//////////////////////////////////////////////////////////////////////////////

void RD_SplitStrategy::setCurrentCell()
{
  DistributionData& ddata = getMethodData().getDistributionData();
  getMethodData().getLinearizer()->setUpdateStates(ddata.states);
  ddata.tStates = computeConsistentStates(ddata.states);
}

//////////////////////////////////////////////////////////////////////////////

void RD_SplitStrategy::setup()
{
  FluctuationSplitStrategy::setup();

  cf_assert( !getMethodData().isMultipleSplitter() );
  
  m_splitter = getMethodData().getSplitter();
  getMethodData().getUpdateVar()->setExtraData(true);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
RD_SplitStrategy::needsSockets()
{
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
     result = FluctuationSplitStrategy::needsSockets();

   return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
