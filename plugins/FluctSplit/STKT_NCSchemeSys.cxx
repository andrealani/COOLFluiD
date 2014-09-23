#include "STKT_NCSchemeSys.hh"
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

MethodStrategyProvider<STKT_NCSchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeRiccNcSchemeSysProvider("STKT_SysNC");

//////////////////////////////////////////////////////////////////////////////

STKT_NCSchemeSys::STKT_NCSchemeSys(const std::string& name) :
  STKT_SplitterSys(name),
  m_sumKplus(),
  _invK(),
  m_uTemp(),
  m_diss(),
  m_sumKplusU()
{
}

//////////////////////////////////////////////////////////////////////////////

STKT_NCSchemeSys::~STKT_NCSchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NCSchemeSys::setup()
{
  STKT_SplitterSys::setup();
 
  // Do not need a dissipation contribution from past 
 getMethodData().getDistributionData().needDiss = false;

  m_sumKplus.resize(_nbEquations, _nbEquations);
  m_uTemp.resize(_nbEquations);
  m_sumKplusU.resize(_nbEquations);
  m_diss.resize(_nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NCSchemeSys::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NCSchemeSys::distribute(vector<RealVector>& residual)
{  
  DistributionData& ddata = getMethodData().getDistributionData();
  const vector<State*>& tStates = *ddata.tStates;
  
  m_sumKplus = *_kPlus[0];
  m_sumKplusU = *_kPlus[0]*(*tStates[0]);
  
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    m_sumKplus  += *_kPlus[iState];
    m_sumKplusU += *_kPlus[iState]*(*tStates[iState]);
  }
  
  _inverter->invert(m_sumKplus, _invK);
  
  m_diss  = _invK*m_sumKplusU;
  m_uTemp = _invK*ddata.phi;
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++ iState){ 
    const RealMatrix& kp = *_kPlus[iState];
    residual[iState] = kp*(m_uTemp + (*tStates[iState]) - m_diss);
    
    // AL: compute betas on the fly is more expensive because
    // it involves a matrix*matrix and a matrix*vector
    // while normal scheme requires only two matrix*vector operations
    RealMatrix& betaLDA = (*ddata.currBetaMat)[iState];
    betaLDA = kp*_invK;
  }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
