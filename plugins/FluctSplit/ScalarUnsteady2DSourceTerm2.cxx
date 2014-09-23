#include <numeric>

#include "ScalarUnsteady2DSourceTerm2.hh"
#include "InwardNormalsData.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ScalarUnsteady2DSourceTerm2,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitSpaceTimeModule>
su22DSTProvider("ScalarUnsteady2DSource2");

//////////////////////////////////////////////////////////////////////////////

void ScalarUnsteady2DSourceTerm2::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("SourceCoeff","Coefficient of the source");
}

//////////////////////////////////////////////////////////////////////////////

ScalarUnsteady2DSourceTerm2::ScalarUnsteady2DSourceTerm2(const std::string& name) :
  ComputeSourceTermFSM(name),
  m_alpha()
{
  addConfigOptionsTo(this);
   m_alpha = 0.0;
   setParameter("SourceCoeff",&m_alpha);
}

//////////////////////////////////////////////////////////////////////////////

ScalarUnsteady2DSourceTerm2::~ScalarUnsteady2DSourceTerm2()
{
}

//////////////////////////////////////////////////////////////////////////////

void ScalarUnsteady2DSourceTerm2::setup()
{
  
}

//////////////////////////////////////////////////////////////////////////////

void ScalarUnsteady2DSourceTerm2::computeSourceFSM(Framework::GeometricEntity *const cell,
					 RealVector& source,
					 const FluctSplit::InwardNormalsData& normalsData)
{
 
  const RealVector centroid = cell->computeCentroid();
   DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT()/2.0;
  DistributionData& ddata = getMethodData().getDistributionData();
  const vector<State*>& states = *ddata.states;

  const CFuint nbStatesInCell = states.size();
  CFreal alpha = m_alpha/nbStatesInCell;

  source[0] = (*states[0])[0];
  
  for( CFuint iState = 1; iState < nbStatesInCell; ++iState)
    source[0] += (*states[iState])[0] ;

  source[0] *= alpha*volumes[cell->getID()]*dt;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
