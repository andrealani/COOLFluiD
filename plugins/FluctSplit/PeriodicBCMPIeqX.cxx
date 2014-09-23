#include "mpi.h"
#include <iostream>
#include <map>

#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "Common/BadValueException.hh"
#include "Common/PE.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/PeriodicBCMPIeqX.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PeriodicBCMPIeqX,
                      FluctuationSplitData,
                      FluctSplitModule>
periodicBCmpieqXProvider("PeriodicBCMPIeqX");

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCMPIeqX::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Threshold","Threshold of distance that will consider two states matching after the transformation of coordinates");
}

//////////////////////////////////////////////////////////////////////////////

PeriodicBCMPIeqX::PeriodicBCMPIeqX(const std::string& name) :
  FluctuationSplitCom(name),
  socket_states("states"),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated")
{
  addConfigOptionsTo(this);

  m_threshold = 1e-12;
  setParameter("Threshold",&m_threshold);
}

//////////////////////////////////////////////////////////////////////////////

PeriodicBCMPIeqX::~PeriodicBCMPIeqX()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
PeriodicBCMPIeqX::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isUpdated);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCMPIeqX::setup()
{
  FluctuationSplitCom::setup();
  CFuint nP = PE::GetPE().GetProcessorCount();
  CFuint myP = PE::GetPE().GetRank();

  _comm = PE::GetPE().GetCommunicator();
  _my_nP = myP;
  _n_P = nP;
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // get the data handle for the is Updated flags
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  // First we get all the coordinates of all the states of the main TRS
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Common::SafePtr< vector<CFuint> > const applied_trs_statesIdx = trs->getStatesInTrs();
  CFuint nb_st_trs = applied_trs_statesIdx->size();
  m_nb_st_trs = nb_st_trs;
  std::vector<CFint> nFpP(nP,-1);
  MPI_Allgather(&nb_st_trs, 1, MPI_UNSIGNED, &nFpP[0], 1, MPI_UNSIGNED, _comm);
  _countFpP = nFpP;


  const CFuint nEquation = PhysicalModelStack::getActive()->getNbEq();
   _nE = nEquation;
  std::vector<CFreal> NodeCoordinate;
  std::vector<bool> NodeParUp;
  NodeCoordinate.reserve(2*nb_st_trs);
  for (CFuint iState = 0; iState< nb_st_trs; iState++){
   CFLogDebugMed( "iState = " << iState << "\n");	
   const CFuint stateID = (*applied_trs_statesIdx)[iState];
    State *const state = states[stateID];

    NodeCoordinate.push_back((*state).getCoordinates()[XX]);
    NodeCoordinate.push_back((*state).getCoordinates()[YY]);
    if ((*state).isParUpdatable()){
    NodeCoordinate.push_back(1.0); // need to fill a vector that store if the state is updatable
    }
    else NodeCoordinate.push_back(0.0);
     
   }


  CFuint nC = 3; //total number of data of the state
  const CFuint dimq = _n_P*nb_st_trs*nC;
  const CFuint dimq2 = nb_st_trs*nC;
  std::vector<CFreal> SendNodeCoordinate;
  SendNodeCoordinate.reserve(dimq);
  for(CFuint p=0; p<_n_P; p++){
    for(CFuint pp=0; pp<dimq2; pp++){
      SendNodeCoordinate.push_back(NodeCoordinate[pp]);
    }
  }

  vector<int> sendcounts(_n_P);
  vector<int> recvcounts(_n_P);
  vector<int> rdispls(_n_P);
  vector<int> sdispls(_n_P);
  rdispls[0] = 0;
  sdispls[0] = 0;
  CFuint LastDisplacement;

  for(CFuint p=0; p<_n_P; p++){
    sendcounts[p] = nC*_countFpP[_my_nP];
    recvcounts[p] = nC*_countFpP[p];
  }
  if(_n_P>1){
    for(CFuint p=1; p<_n_P; p++){
      rdispls[p] = rdispls[p-1] + recvcounts[p-1];
      sdispls[p] = sdispls[p-1] + sendcounts[p-1];
    }
  }


  _ConnectionStatePeriodic.reserve( nb_st_trs);
  _ConnectionProcessPeriodic.reserve( nb_st_trs);
  vector<CFint> RecvDis(rdispls);
  LastDisplacement = rdispls[_n_P-1] + recvcounts[_n_P-1];
  vector<double> rbuf(LastDisplacement,0.0);
  vector<double> sbuf(SendNodeCoordinate);
  MPI_Alltoallv(&sbuf[0], &sendcounts[0], &sdispls[0], MPI_DOUBLE, &rbuf[0], &recvcounts[0], &rdispls[0], MPI_DOUBLE, _comm);
  std::vector<CFreal> RecvNodeCoordinate(rbuf);


  for(CFuint p=0; p<nP; ++p){
    CFuint count = _countFpP[p];
    for(CFuint iState = 0; iState<nb_st_trs; iState++){
      
            const CFuint StateID = (*applied_trs_statesIdx)[iState];
      State *const state = states[StateID];
      CFuint stateGolbalID = state->getGlobalID();

      CFreal nodeX = state->getCoordinates()[XX];
      CFreal nodeY = state->getCoordinates()[YY];
	for (CFuint iPeriodicState = 0; iPeriodicState<count; ++iPeriodicState){
	  CFreal nodePeriodicXCoordinate = RecvNodeCoordinate[RecvDis[p] + nC*iPeriodicState + 0];
	  CFreal nodePeriodicYCoordinate = RecvNodeCoordinate[RecvDis[p] + nC*iPeriodicState + 1];
	  CFreal isParUpdate = RecvNodeCoordinate[RecvDis[p] + nC*iPeriodicState + 2];

	  if(MathTools::MathChecks::isEqualWithError(nodeX, nodePeriodicXCoordinate, m_threshold)) {
	    if(MathTools::MathChecks::isNotEqualWithError(nodeY,nodePeriodicYCoordinate, m_threshold)){
	      if(isParUpdate == 1.0){

		_ConnectionStatePeriodic.insert(stateGolbalID,iPeriodicState);
		_ConnectionProcessPeriodic.insert(stateGolbalID,p);
		break;
	      
	    }
	  }
	}
      }
  }
  }
  _ConnectionStatePeriodic.sortKeys();
  _ConnectionProcessPeriodic.sortKeys();
 _sendcounts2.resize(_n_P);
 _sendcounts3.resize(_n_P);
  _recvcounts2.resize(_n_P);
  _recvcounts3.resize(_n_P);
  _rdispls2.resize(_n_P);
  _rdispls3.resize(_n_P);
  _sdispls2.resize(_n_P);
  _sdispls3.resize(_n_P);
  _rdispls2[0] = 0;
  _rdispls3[0] = 0;
  _sdispls2[0] = 0;
  _sdispls3[0] = 0;
  _LastDisplacement = 0;
  _LastDisplacement3 = 0;

  for(CFuint p=0; p<_n_P; p++){
    _sendcounts2[p] = _nE*_countFpP[_my_nP];
    _recvcounts2[p] = _nE*_countFpP[p];
    _sendcounts3[p] = _countFpP[_my_nP];
    _recvcounts3[p] = _countFpP[p];

  }
  if(_n_P>1){
    for(CFuint p=1; p<_n_P; p++){
      _rdispls2[p] = _rdispls2[p-1] + _recvcounts2[p-1];
      _sdispls2[p] = _sdispls2[p-1] + _sendcounts2[p-1];
      _rdispls3[p] = _rdispls3[p-1] + _recvcounts3[p-1];
      _sdispls3[p] = _sdispls3[p-1] + _sendcounts3[p-1];
    }
  }
  _LastDisplacement = _rdispls2[_n_P-1] + _recvcounts2[_n_P-1];
  _rbuf2.resize(_LastDisplacement);
  _LastDisplacement3 = _rdispls3[_n_P-1] + _recvcounts3[_n_P-1];
  _rbuf3.resize(_LastDisplacement3);
  _Boundaryrhs.resize(_LastDisplacement);
  _Boundaryupcoef.resize(_LastDisplacement3);
  for(CFuint i=0; i<_LastDisplacement; i++){
    _rbuf2[i] = 0.0;
    _Boundaryrhs[i] = 0.0;
  }

  for(CFuint i=0; i<_LastDisplacement3; i++){
    _rbuf3[i] = 0.0;
    _Boundaryupcoef[i] = 0.0;
  }
  for(CFuint i=0; i<_n_P; i++){
    _RecvDis2.push_back(_rdispls2[i]);
    _RecvDis3.push_back(_rdispls3[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCMPIeqX::executeOnTrs()
{
  DataHandle<CFreal>  rhs         = socket_rhs.getDataHandle();
  DataHandle<CFreal>  updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle< Framework::State*, Framework::GLOBAL > m_states = socket_states.getDataHandle();
  Common::SafePtr< vector<CFuint> > const applied_trs_statesIdx = getCurrentTRS()->getStatesInTrs();
 
  // get the data handle for the is Updated flags
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  
  preProcess();

 for (CFuint is = 0; is < applied_trs_statesIdx->size(); ++is)
  {
    
  const CFuint applied_sid = (*applied_trs_statesIdx)[is];
  State *const state = m_states[applied_sid];
  const CFuint globalID = state->getGlobalID();

  const CFuint ProcessPeriodic = _ConnectionProcessPeriodic.find(globalID);
  const CFuint FacePeriodic = _ConnectionStatePeriodic.find(globalID);

 
  if ((state->isParUpdatable())){
    for(CFuint h=0; h<_nE; h++){
      rhs(applied_sid,h,_nE) += _rbuf2[_RecvDis2[ProcessPeriodic] +_nE*FacePeriodic+h];
      
    }
    
    updateCoeff[applied_sid] += _rbuf3[_RecvDis3[ProcessPeriodic] + FacePeriodic];
    
  }
  
    


  }

}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCMPIeqX::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
 
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCMPIeqX::preProcess()
{
  const CFuint dimq = m_nb_st_trs*_nE*_n_P;
  const CFuint dimq2 = m_nb_st_trs*_nE;// <<-----------------------------------with allgatherv

  DataHandle<CFreal>  rhs         = socket_rhs.getDataHandle();
  DataHandle<CFreal>  updateCoeff = socket_updateCoeff.getDataHandle();

  std::vector<CFreal> Boundaryrhs(dimq2);
  std::vector<CFreal> Boundaryupcoef( m_nb_st_trs);
  std::vector<CFreal> Boundaryupcoef2( m_nb_st_trs);

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  Common::SafePtr< vector<CFuint> > const applied_trs_statesIdx = getCurrentTRS()->getStatesInTrs();

  for(CFuint iState = 0; iState < m_nb_st_trs; iState++){
       const CFuint stateID = (*applied_trs_statesIdx)[iState];

       State *const state = states[stateID];

    for(CFuint e=0; e<_nE; e++){
      Boundaryrhs[_nE*iState + e] = rhs(stateID,e,_nE);
      
    }
 

    Boundaryupcoef[iState] = updateCoeff[stateID];
  
  }
 
  std::vector<CFreal> SendBoundaryrhs;
  std::vector<CFreal> SendBoundaryupcoef;
  SendBoundaryrhs.reserve(dimq);
  SendBoundaryupcoef.reserve(m_nb_st_trs*_n_P);
  for(CFuint p = 0; p <_n_P; p++){
    for(CFuint pp = 0; pp < m_nb_st_trs*_nE; pp++){
      SendBoundaryrhs.push_back(Boundaryrhs[pp]);
    }

    for(CFuint pp = 0; pp < m_nb_st_trs; pp++){
      SendBoundaryupcoef.push_back(Boundaryupcoef[pp]);
    }

  }
  vector<double> sbuf2(SendBoundaryrhs);
  vector<double> sbuf2AG(Boundaryrhs);// <<-----------------------------------with allgatherv
  vector<double> sbuf3(SendBoundaryupcoef);
  vector<double> sbuf3AG(Boundaryupcoef);// <<-----------------------------------with allgatherv

  MPI_Alltoallv(&sbuf2[0], &_sendcounts2[0], &_sdispls2[0], MPI_DOUBLE, &_rbuf2[0], &_recvcounts2[0], &_rdispls2[0], MPI_DOUBLE, _comm);
  MPI_Alltoallv(&sbuf3[0], &_sendcounts3[0], &_sdispls3[0], MPI_DOUBLE, &_rbuf3[0], &_recvcounts3[0], &_rdispls3[0], MPI_DOUBLE, _comm);

  _Boundaryrhs.clear();
  _Boundaryupcoef.clear();

  _Boundaryrhs.reserve(_LastDisplacement);
  _Boundaryupcoef.reserve(_LastDisplacement3);

  for(CFuint i=0; i<_LastDisplacement; i++){
    _Boundaryrhs.push_back(_rbuf2[i]);
  }

for(CFuint i=0; i<_LastDisplacement3; i++){
    _Boundaryupcoef.push_back(_rbuf3[i]);
    
  }
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCMPIeqX::Barrier()
{
  MPI_Barrier(_comm);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
