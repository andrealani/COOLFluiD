#include "FluctSplit/SourceTermSplitStrategy.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<SourceTermSplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitModule>
sourceTermSplitStrategyProvider("SourceTerm");

//////////////////////////////////////////////////////////////////////////////

SourceTermSplitStrategy::SourceTermSplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  m_unitFaceNormals()
{
}

//////////////////////////////////////////////////////////////////////////////

SourceTermSplitStrategy::~SourceTermSplitStrategy()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void SourceTermSplitStrategy::setup()
{
  // first call parent method
  FluctuationSplitStrategy::setup();
  
  // now complement it
  m_unitFaceNormals.resize(MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell());
  for (CFuint  i = 0; i < m_unitFaceNormals.size(); ++i) {
    m_unitFaceNormals[i].resize(PhysicalModelStack::getActive()->getDim());
  }
}

//////////////////////////////////////////////////////////////////////////////

void SourceTermSplitStrategy::computeFluctuation(vector<RealVector>& residual)
{
  setCurrentCell();

  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DistributionData& ddata = getMethodData().getDistributionData();
  
  const CFuint nbST = getMethodData().getSourceTermSplitter()->size();
  for (CFuint i = 0; i < nbST; ++i) {
    ddata.sourceTermID = i;
    getMethodData().getSourceTermSplitter(i)->computeSourceTerm(*normals[ddata.cellID]);
    getMethodData().getSourceTermSplitter(i)->distribute(residual);
  }
}

//////////////////////////////////////////////////////////////////////////////

void SourceTermSplitStrategy::setCurrentCell()
{
  DistributionData& ddata = getMethodData().getDistributionData();
  ddata.tStates = computeConsistentStates(ddata.states);
  cf_assert(ddata.tStates == ddata.states);
}
      
//////////////////////////////////////////////////////////////////////////////

void SourceTermSplitStrategy::unsetup()
{
  // first call parent method
  FluctuationSplitStrategy::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD
