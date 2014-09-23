#include <numeric>

#include "ScalarRotationadv2DSourceTerm.hh"
#include "InwardNormalsData.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/FluctSplitRotationDiffusion.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ScalarRotationadv2DSourceTerm,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitRotationDiffusionModule>
sra2DSProvider("ScalarRotationadv2DSource");

//////////////////////////////////////////////////////////////////////////////

void ScalarRotationadv2DSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  
}

//////////////////////////////////////////////////////////////////////////////

ScalarRotationadv2DSourceTerm::ScalarRotationadv2DSourceTerm(const std::string& name) :
  ComputeSourceTermFSM(name),
  m_alpha()
{
   addConfigOptionsTo(this);
   m_alpha = 0.0;
  
}

//////////////////////////////////////////////////////////////////////////////

ScalarRotationadv2DSourceTerm::~ScalarRotationadv2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void ScalarRotationadv2DSourceTerm::setup()
{
  
}

//////////////////////////////////////////////////////////////////////////////

void ScalarRotationadv2DSourceTerm::computeSourceFSM(Framework::GeometricEntity *const cell,
					 RealVector& source,
					 const FluctSplit::InwardNormalsData& normalsData)
{
  const vector<State*>& states = *cell->getStates();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  const CFuint nbStatesInCell = states.size();


   const  Node& node = states[0]->getCoordinates();
    CFreal r = sqrt((node[XX] - 0.5)*(node[XX] - 0.5) + (node[YY] - 0.3)*(node[YY] - 0.3));
    if (r < 0.25)
      source[0] = 10.0;
    else
      source[0] = 0.0;
  

  for (CFuint iState = 1; iState < nbStatesInCell; ++iState)
  { 
    const  Node&  node = states[iState]->getCoordinates();
    CFreal r = sqrt((node[XX] - 0.5)*(node[XX] - 0.5) + (node[YY] - 0.3)*(node[YY] - 0.3));
    if (r < 0.25)
      source[0] += 10.0;
    else
      source[0] += 0.0;
  }

  source[0] *= (volumes[cell->getID()])/3.0;



    
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
