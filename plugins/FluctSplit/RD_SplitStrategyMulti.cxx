#include "FluctSplit/RD_SplitStrategyMulti.hh"

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

MethodStrategyProvider<RD_SplitStrategyMulti,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitModule>
                       rdsSplitStrategyMultiProvider("RDMulti");

//////////////////////////////////////////////////////////////////////////////

RD_SplitStrategyMulti::RD_SplitStrategyMulti(const std::string& name) :
  FluctuationSplitStrategy(name),
  m_systemSplitter(CFNULL),
  m_scalarSplitter(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

RD_SplitStrategyMulti::~RD_SplitStrategyMulti()
{
}

//////////////////////////////////////////////////////////////////////////////

void RD_SplitStrategyMulti::
computeFluctuation(vector<RealVector>& residual)
{
  setCurrentCell();

  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DistributionData& ddata = getMethodData().getDistributionData();
  
  // compute the residual and the upwind parameters k in this cell
  m_systemSplitter->computeK(*ddata.states,normals[ddata.cellID]);
  m_scalarSplitter->computeK(*ddata.states,normals[ddata.cellID]);
  
  // distribute the residual
  m_systemSplitter->distributePart(residual);
  m_scalarSplitter->distributePart(residual);
}

//////////////////////////////////////////////////////////////////////////////

void RD_SplitStrategyMulti::setup()
{
  FluctuationSplitStrategy::setup();

  cf_assert(getMethodData().isMultipleSplitter());
  
  m_scalarSplitter = getMethodData().getScalarSplitter();
  m_systemSplitter = getMethodData().getSysSplitter();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
RD_SplitStrategyMulti::needsSockets()
{
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
     result = FluctuationSplitStrategy::needsSockets();

   return result;
}

//////////////////////////////////////////////////////////////////////////////

void RD_SplitStrategyMulti::setCurrentCell()
{
  DistributionData& ddata = getMethodData().getDistributionData();
  getMethodData().getLinearizer()->setUpdateStates(ddata.states);
  ddata.tStates = computeConsistentStates(ddata.states);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
