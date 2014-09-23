#include <numeric>

#include "ScalarUnsteady2DSourceTermHO.hh"
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

MethodStrategyProvider<ScalarUnsteady2DSourceTermHO,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitSpaceTimeModule>
su2DSThoProvider("ScalarUnsteady2DSourceHO");

//////////////////////////////////////////////////////////////////////////////

void ScalarUnsteady2DSourceTermHO::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("SourceCoeff","Coefficient of the source");
}

//////////////////////////////////////////////////////////////////////////////

ScalarUnsteady2DSourceTermHO::ScalarUnsteady2DSourceTermHO(const std::string& name) :
  ComputeSourceTermFSM(name),
  m_alpha()
{
  addConfigOptionsTo(this);
   m_alpha = 0.0;
   setParameter("SourceCoeff",&m_alpha);
}

//////////////////////////////////////////////////////////////////////////////

ScalarUnsteady2DSourceTermHO::~ScalarUnsteady2DSourceTermHO()
{
}

//////////////////////////////////////////////////////////////////////////////

void ScalarUnsteady2DSourceTermHO::setup()
{
  source_value.resize(6);

  qd0.resize(4); // quadrature points per sub-triangle
  qd1.resize(4); // quadrature points per sub-triangle
  qd2.resize(4); // quadrature points per sub-triangle

  qd0[0] = (1.0/3.0);  qd1[0] = qd0[0]; qd2[0] = qd0[0];
  qd0[1] = 0.6      ;  qd1[1] = 0.2   ; qd2[1] = 0.2   ;
  qd0[2] = 0.2      ;  qd1[2] = 0.6   ; qd2[2] = 0.2   ;
  qd0[3] = 0.2      ;  qd1[3] = 0.2   ; qd2[3] = 0.6   ;
  
  wqd.resize(4); // 3 quadrature points per face
  wqd[0] = -27.0/48.0;
  wqd[1] = 25.0/48.0;
  wqd[2] = 25.0/48.0;
  wqd[3] = 25.0/48.0;

  qdstates.resize(4); // 3 quadrature points per face
  qdstates[0] = new State();
  qdstates[1] = new State();
  qdstates[2] = new State();
  qdstates[3] = new State();
}

//////////////////////////////////////////////////////////////////////////////

void ScalarUnsteady2DSourceTermHO::computeSourceFSM(Framework::GeometricEntity *const cell,
					 RealVector& source,
					 const FluctSplit::InwardNormalsData& normalsData)
{
 
  DistributionData& ddata = getMethodData().getDistributionData();
  const vector<State*>& states = *ddata.states;
  const CFreal dt = SubSystemStatusStack::getActive()->getDT()/2.0;
  std::vector<Framework::State*>& substates = *ddata.subStates;
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  const CFuint nbStatesInCell = states.size();

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState)
  { 
      source_value[iState] = (*states[iState])[0];
 
  }

  const Node& nodeSub0 = substates[0]->getCoordinates();
  const Node& nodeSub1 = substates[1]->getCoordinates();
  const Node& nodeSub2 = substates[2]->getCoordinates();

  const Node& node0 = states[0]->getCoordinates();
  const Node& node1 = states[1]->getCoordinates();
  const Node& node2 = states[2]->getCoordinates();

  const CFreal nx1 = normalsData.getNodalNormComp(0,XX);
  const CFreal nx2 = normalsData.getNodalNormComp(1,XX);
  const CFreal nx3 = normalsData.getNodalNormComp(2,XX);

  const CFreal ny1 = normalsData.getNodalNormComp(0,YY);
  const CFreal ny2 = normalsData.getNodalNormComp(1,YY);
  const CFreal ny3 = normalsData.getNodalNormComp(2,YY);

  const CFreal inv_volume = 1.0 / ddata.cell->computeVolume();         
  const CFreal volume = ddata.cell->computeVolume();  


  CFreal x = qd0[0]*nodeSub0[XX] + qd1[0]*nodeSub1[XX] + qd2[0]*nodeSub2[XX];
  CFreal y = qd0[0]*nodeSub0[YY] + qd1[0]*nodeSub1[YY] + qd2[0]*nodeSub2[YY];
             
  
  CFreal L1 = 1.0 + 0.5*( ( x - node0[XX] )*nx1 + ( y - node0[YY] )*ny1 )*inv_volume;
  CFreal L2 = 1.0 + 0.5*( ( x - node1[XX] )*nx2 + ( y - node1[YY] )*ny2 )*inv_volume;                                   
  CFreal L3 = 1.0 + 0.5*( ( x - node2[XX] )*nx3 + ( y - node2[YY] )*ny3 )*inv_volume;
 
  source[0] = (L1*(2.0*L1 - 1.0)*source_value[0] +  L2*(2.0*L2 - 1.0)*source_value[1] +
     L3*(2.0*L3 - 1.0)*source_value[2] + 4.0*L1*L2*source_value[3] + 4.0*L2*L3*source_value[4] + 
    4.0*L1*L3*source_value[5])*wqd[0]*volume*0.25;

  for (CFuint iQd = 1; iQd < 4; ++iQd)
    {        
  
      x = qd0[iQd]*nodeSub0[XX] + qd1[iQd]*nodeSub1[XX] + qd2[iQd]*nodeSub2[XX];
      y = qd0[iQd]*nodeSub0[YY] + qd1[iQd]*nodeSub1[YY] + qd2[iQd]*nodeSub2[YY];
             
  
      L1 = 1.0 + 0.5*( ( x - node0[XX] )*nx1 + ( y - node0[YY] )*ny1 )*inv_volume;
      L2 = 1.0 + 0.5*( ( x - node1[XX] )*nx2 + ( y - node1[YY] )*ny2 )*inv_volume;                                   
      L3 = 1.0 + 0.5*( ( x - node2[XX] )*nx3 + ( y - node2[YY] )*ny3 )*inv_volume;
  
      source[0] += (L1*(2.0*L1 - 1.0)*source_value[0] +  L2*(2.0*L2 - 1.0)*source_value[1] +
        L3*(2.0*L3 - 1.0)*source_value[2] + 4.0*L1*L2*source_value[3] + 4.0*L2*L3*source_value[4] + 
        4.0*L1*L3*source_value[5])*wqd[iQd]*volume*0.25;

      }  
    source[0] *= m_alpha*dt;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
