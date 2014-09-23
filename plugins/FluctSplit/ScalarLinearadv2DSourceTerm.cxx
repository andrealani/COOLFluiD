#include <numeric>

#include "ScalarLinearadv2DSourceTerm.hh"
#include "InwardNormalsData.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/FluctSplitAdvectionDiffusion.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ScalarLinearadv2DSourceTerm,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitAdvectionDiffusionModule>
sla2DSProvider("ScalarLinearadv2DSource");

//////////////////////////////////////////////////////////////////////////////

void ScalarLinearadv2DSourceTerm::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("SourceCoeff","Coefficient of the source");
}

//////////////////////////////////////////////////////////////////////////////

ScalarLinearadv2DSourceTerm::ScalarLinearadv2DSourceTerm(const std::string& name) :
  ComputeSourceTermFSM(name),
  m_alpha()
{
   addConfigOptionsTo(this);
   m_alpha = 0.0;
   setParameter("SourceCoeff",&m_alpha);
}

//////////////////////////////////////////////////////////////////////////////

ScalarLinearadv2DSourceTerm::~ScalarLinearadv2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void ScalarLinearadv2DSourceTerm::setup()
{
  
}

//////////////////////////////////////////////////////////////////////////////

void ScalarLinearadv2DSourceTerm::computeSourceFSM(Framework::GeometricEntity *const cell,
					 RealVector& source,
					 const FluctSplit::InwardNormalsData& normalsData)
{
//   DistributionData& ddata = getMethodData().getDistributionData();
  
  const vector<State*>& states = *cell->getStates();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  const CFuint nbStatesInCell = states.size();
  CFreal alpha = m_alpha;
 
  source[0] = (*states[0])[0];
  
  for( CFuint iState = 1; iState < nbStatesInCell; ++iState)
    source[0] += (*states[iState])[0] ;

  source[0] *= alpha*(volumes[cell->getID()])/3.0;


    
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
