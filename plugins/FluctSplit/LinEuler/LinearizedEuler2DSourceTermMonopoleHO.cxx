#include <numeric>

#include "LinearizedEuler2DSourceTermMonopoleHO.hh"
#include "FluctSplit/InwardNormalsData.hh"
//#include "Environment/CFLog.hh"
#include "Framework/State.hh"
#include "Framework/MeshData.hh"
#include "FluctSplitLinEuler.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitLinEuler.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LinearizedEuler2DSourceTermMonopoleHO,
		       FluctuationSplitData,
		       ComputeSourceTermFSM,
		       FluctSplitLinEulerModule>
linEuler2DSTProviderMonopoleHO("LinEuler2DSourceMonopoleHO");

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermMonopoleHO::LinearizedEuler2DSourceTermMonopoleHO(const std::string& name) :
  ComputeSourceTermFSM(name),
  _varSet(CFNULL),
  _temp(),
  m_alpha(-1.0),
  m_eps(-1.0),
  m_freq(-1.0),
  _sourceloc(0),
  qd0(),
  qd1(),
  qd2(),
  wqd(),
  qdstates(),
  source_value()
{
  addConfigOptionsTo(this);

  m_alpha = log(2.)/5.;
  setParameter("alpha", &m_alpha);

  m_eps = 0.01;
  setParameter("eps", &m_eps);

  m_freq = 30.;
  setParameter("freq", &m_freq);

  setParameter("location", &_sourceloc);
}

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermMonopoleHO::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("alpha","Width of the source.");
  options.addConfigOption< CFreal >("eps","Amplitude of the source.");
  options.addConfigOption< CFreal >("freq","Frequency of the source.");
  options.addConfigOption< std::vector<CFreal> >("location","Location of the source.");
}

//////////////////////////////////////////////////////////////////////////////

LinearizedEuler2DSourceTermMonopoleHO::~LinearizedEuler2DSourceTermMonopoleHO()
{
}
/*
//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermMonopoleHO::unsetup()
{
  for (CFuint i = 0; i < qdstates.size(); ++i) {
    deletePtr(qdstates[i]);
    }
}*/

//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermMonopoleHO::setup()
{
  _temp.resize(PhysicalModelStack::getActive()->getNbEq());

//  source_value.resize(6);

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

  
  source_value.resize(6);
  const CFuint nbEq = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iState = 0; iState < 6; iState++ )
    {
      source_value[iState].resize(nbEq);
    }

  _varSet = getMethodData().getUpdateVar().d_castTo<LinEuler2DVarSet>();
}


//////////////////////////////////////////////////////////////////////////////

void LinearizedEuler2DSourceTermMonopoleHO::computeSourceFSM(Framework::GeometricEntity *const cell,
					 RealVector& source,
					 const FluctSplit::InwardNormalsData& normalsData)
{
  
  DistributionData& ddata = getMethodData().getDistributionData();
  const CFreal time = ddata.time;
  const CFuint cellID = ddata.cellID;
  const CFuint nbStatesInCell = ddata.states->size();
  bool First_itP1 = ddata.isfirstP1;
  vector<State*>& states = *ddata.states;
  std::vector<Framework::State*>& substates = *ddata.subStates;
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  CFuint nbEq = PhysicalModelStack::getActive()->getNbEq();
  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];
  const CFreal pi=3.14159265;
  const CFreal alpha = m_alpha;
  const CFreal eps = m_eps;
  const CFreal om = 2.*pi/m_freq;
  const CFreal xs = _sourceloc[0];
  const CFreal ys = _sourceloc[1];
   if (!First_itP1){

  CFreal r;
  CFreal f;
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState)
  { 
    const CFreal x = (states[iState]->getCoordinates())[XX];
    const CFreal y = (states[iState]->getCoordinates())[YY];
    r = alpha*((x-xs)*(x-xs)+(y-ys)*(y-ys));
    if (r >= 75.0)
      f = 0;
    else f=eps*exp(-r);
    
    source_value[iState][0] = f*sin(om*time);
    source_value[iState][1] = 0.0;
    source_value[iState][2] = 0.0;
    source_value[iState][3] = c*c*f*sin(om*time);
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
  const CFreal volume =ddata.cell->computeVolume();  
  CFreal x = qd0[0]*nodeSub0[XX] + qd1[0]*nodeSub1[XX] + qd2[0]*nodeSub2[XX];
  CFreal y = qd0[0]*nodeSub0[YY] + qd1[0]*nodeSub1[YY] + qd2[0]*nodeSub2[YY];
             
  
  CFreal L1 = 1.0 + 0.5*( ( x - node0[XX] )*nx1 + ( y - node0[YY] )*ny1 )*inv_volume;
  CFreal L2 = 1.0 + 0.5*( ( x - node1[XX] )*nx2 + ( y - node1[YY] )*ny2 )*inv_volume;
  CFreal L3 = 1.0 + 0.5*( ( x - node2[XX] )*nx3 + ( y - node2[YY] )*ny3 )*inv_volume;
 
  for (CFuint iEq =0; iEq < nbEq; ++iEq){
    source[iEq] = (L1*(2.0*L1 - 1.0)*source_value[0][iEq] +  L2*(2.0*L2 - 1.0)*source_value[1][iEq] +
       L3*(2.0*L3 - 1.0)*source_value[2][iEq]+ 4.0*L1*L2*source_value[3][iEq] + 4.0*L2*L3*source_value[4][iEq] + 
      4.0*L1*L3*source_value[5][iEq])*wqd[0]*volume*0.25;
  }

  for (CFuint iQd = 1; iQd < 4; ++iQd)
    {        
  
      x = qd0[iQd]*nodeSub0[XX] + qd1[iQd]*nodeSub1[XX] + qd2[iQd]*nodeSub2[XX];
      y = qd0[iQd]*nodeSub0[YY] + qd1[iQd]*nodeSub1[YY] + qd2[iQd]*nodeSub2[YY];
             
  
      L1 = 1.0 + 0.5*( ( x - node0[XX] )*nx1 + ( y - node0[YY] )*ny1 )*inv_volume;
      L2 = 1.0 + 0.5*( ( x - node1[XX] )*nx2 + ( y - node1[YY] )*ny2 )*inv_volume;                                   
      L3 = 1.0 + 0.5*( ( x - node2[XX] )*nx3 + ( y - node2[YY] )*ny3 )*inv_volume;
  
      for (CFuint iEq =0; iEq < nbEq; ++iEq){
        source[iEq] += (L1*(2.0*L1 - 1.0)*source_value[0][iEq] +  L2*(2.0*L2 - 1.0)*source_value[1][iEq] +
          L3*(2.0*L3 - 1.0)*source_value[2][iEq] + 4.0*L1*L2*source_value[3][iEq] + 4.0*L2*L3*source_value[4][iEq] + 
          4.0*L1*L3*source_value[5][iEq])*wqd[iQd]*volume*0.25;
      }
}
}
else{
  source[0] = 0.0;
  source[1] = 0.0;
  source[2] = 0.0;
  source[3] = 0.0;
 
  CFreal r;
  CFreal f;

const CFreal volume =ddata.cell->computeVolume();
 for (CFuint iState = 0; iState < 3; ++iState)
  {
    const CFreal x = (substates[iState]->getCoordinates())[XX];
    const CFreal y = (substates[iState]->getCoordinates())[YY];
    r = alpha*((x-xs)*(x-xs)+(y-ys)*(y-ys));
    if (r >= 80.0){
      f = 0;
    //  CF_DEBUG_POINT;
    }
    else{
     f=eps*exp(-r);
  //   CF_DEBUG_OBJ(f);
    }

    source[0] += f*sin(om*time);
    source[1] += 0.0;
    source[2] += 0.0;
    source[3] += c*c*f*sin(om*time);


  }

  source *= (volume*0.25)/3.0;
 // if (source[0] != 0.0)
   // CF_DEBUG_OBJ(source[0]);
}
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
