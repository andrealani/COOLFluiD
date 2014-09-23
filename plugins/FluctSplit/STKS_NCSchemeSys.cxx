#include "STKS_NCSchemeSys.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "MathTools/MatrixInverter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STKS_NCSchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
stksNcSchemeSysProvider("STKS_SysNC");

//////////////////////////////////////////////////////////////////////////////

STKS_NCSchemeSys::STKS_NCSchemeSys(const std::string& name) :
  STKS_SplitterSys(name),
  m_sumKplus(),
  m_sumKplusU(),
  m_uTemp(),
  m_diss(),
  _invK(),
  past_diss(),
  m_time_comp()
{
}

//////////////////////////////////////////////////////////////////////////////

STKS_NCSchemeSys::~STKS_NCSchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKS_NCSchemeSys::setup()
{
   STKS_SplitterSys::setup();

  m_sumKplus.resize(_nbEquations, _nbEquations);
  m_uTemp.resize(_nbEquations);
  m_sumKplusU.resize(_nbEquations);
  m_diss.resize(_nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
  getMethodData().getDistributionData().needDiss = true;
const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

past_diss.resize(maxNbStatesInCell);

  for (CFuint iState = 0; iState < maxNbStatesInCell; ++ iState){
   
   past_diss[iState].resize(_nbEquations);
}
  m_time_comp.resize(maxNbStatesInCell);

  for (CFuint iState = 0; iState < maxNbStatesInCell; ++ iState){
   
   m_time_comp[iState].resize(_nbEquations);
}
}

//////////////////////////////////////////////////////////////////////////////

void STKS_NCSchemeSys::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalar::computePicardJacob()");
}


//////////////////////////////////////////////////////////////////////////////
  void STKS_NCSchemeSys::ComputePastDissipationAndTimeComp( const vector<State*>& tStates )
{
 const CFuint nbStatesInCell = _nbStatesInCell;

 m_sumKplus = *_kPlus[0];
 m_sumKplusU = *_kPlus[0]*(*tStates[0]);

  for (CFuint iState = 1; iState < nbStatesInCell; ++ iState){
    m_sumKplus += *_kPlus[iState];
     m_sumKplusU += *_kPlus[iState]*(*tStates[iState]);
  } 

 _inverter->invert(m_sumKplus, _invK);
  m_uTemp = _invK*m_sumKplusU;

  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
  past_diss[iState]  = *_kPlus[iState]*(*tStates[iState])- (*_kPlus[iState])*m_uTemp;
  m_time_comp[iState] = (-1.)*(*tStates[iState]);
}
}


//////////////////////////////////////////////////////////////////////////////

void STKS_NCSchemeSys::distribute(vector<RealVector>& residual)
{  
  const RealVector& phiT = getMethodData().getDistributionData().phi;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT()/dim;
  const CFreal Area = _cellVolume/(dim+1.);

  m_sumKplus = *_kPlus[0]; 
  m_sumKplusU = *_kPlus[0]*(*tStates[0]);

  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    m_sumKplus  += *_kPlus[iState];
    m_sumKplusU += *_kPlus[iState]*(*tStates[iState]);
  }

  _inverter->invert(m_sumKplus, _invK);
  m_uTemp = _invK*phiT;
  
  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState) { 
    residual[iState] = *_kPlus[iState]*m_uTemp;
  }

  m_uTemp = _invK*m_sumKplusU;
 
  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState) { 
    // And then the dissipation
    m_diss = *_kPlus[iState]*(*tStates[iState]) - *_kPlus[iState]*m_uTemp;

   residual[iState] += timeStep*(m_diss + past_diss[iState]);
  
  residual[iState] += Area*((*tStates[iState]) +  m_time_comp[iState]);
  }

  // AL: compute betas on the fly is more expensive because
  // it involves a matrix*matrix and a matrix*vector
  // while normal scheme requires only two matrix*vector operations
  vector<RealMatrix>& betas = *getMethodData().getDistributionData().currBetaMat;
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    betas[iState] = (*_kPlus[iState])*_invK;
  }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
