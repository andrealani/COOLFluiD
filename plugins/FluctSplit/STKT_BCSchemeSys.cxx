#include "Common/NotImplementedException.hh"

#include "MathTools/MatrixInverter.hh"

#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/STKT_BCSchemeSys.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STKT_BCSchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeRiccBcSchemeSysProvider("STKT_SysBC");

//////////////////////////////////////////////////////////////////////////////

STKT_BCSchemeSys::STKT_BCSchemeSys(const std::string& name) :
  STKT_SplitterSys(name),
  m_sumKplus(),
  _invK(),
  m_uTemp(),
  m_diss(),
  m_sumKplusU(),
  m_phiN(),
  m_phitot(),
  m_absphitot(),
  m_phiT()
{
}

//////////////////////////////////////////////////////////////////////////////

STKT_BCSchemeSys::~STKT_BCSchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKT_BCSchemeSys::setup()
{
  STKT_SplitterSys::setup();
   // Do not need a dissipation contribution from past 
 getMethodData().getDistributionData().needDiss = false;
 CFuint maxNbStatesInCell = _kPlus.size();
  m_phiN.resize(maxNbStatesInCell);
  m_diss.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++ iState){
    m_phiN[iState].resize(_nbEquations);
    m_diss[iState].resize(_nbEquations);
 }
 m_sumKplus.resize(_nbEquations,_nbEquations);
 m_uTemp.resize(_nbEquations);
   m_sumKplusU.resize(_nbEquations);
   m_phitot.resize(_nbEquations);
  m_absphitot.resize(_nbEquations);
 _invK.resize(_nbEquations,_nbEquations);
  m_phiT.resize(_nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void STKT_BCSchemeSys::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void STKT_BCSchemeSys::distribute(vector<RealVector>& residual)
{  
    const RealVector& m_phiT = getMethodData().getDistributionData().phi;
  const CFuint nbStatesInCell = _nbStatesInCell;
 const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
CFreal theta;
  m_sumKplus = *_kPlus[0];
  m_sumKplusU = *_kPlus[0]*(*tStates[0]);

  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    m_sumKplus  += *_kPlus[iState];
    m_sumKplusU += *_kPlus[iState]*(*tStates[iState]);
  }
  _inverter->invert(m_sumKplus, _invK);
  m_uTemp = _invK*m_phiT;
  
  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
    residual[iState] = *_kPlus[iState]*m_uTemp;
    m_phiN[iState] = *_kPlus[iState]*m_uTemp;
  }

  m_uTemp = _invK*m_sumKplusU;
 
  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
    // And then the dissipation
    m_diss[iState] = *_kPlus[iState]*(*tStates[iState]) - *_kPlus[iState]*m_uTemp;

    m_phiN[iState] += m_diss[iState];
  } 

  //computation of theta
  m_phitot = m_phiN[0];
  for(CFuint iEq = 0; iEq < _nbEquations; ++iEq)
    m_absphitot[iEq] = std::abs(m_phiN[0][iEq]);

  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
    m_phitot += m_phiN[iState];
    for(CFuint iEq = 0; iEq < _nbEquations; ++iEq)
      m_absphitot[iEq] += std::abs(m_phiN[iState][iEq]);
  }

  for(CFuint iEq = 0; iEq < _nbEquations; ++iEq){
    theta =  std::abs(m_phitot[iEq])/max(MathTools::MathConsts::CFrealEps(), m_absphitot[iEq]);
    for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
      residual[iState][iEq] += theta*m_diss[iState][iEq];
  
    }
  }

    // AL: compute betas on the fly is more expensive because
    // it involves a matrix*matrix and a matrix*vector
    // while normal scheme requires only two matrix*vector operations
    vector<RealMatrix>& betas = *getMethodData().getDistributionData().currBetaMat;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
    {
      betas[iState] = (*_kPlus[iState])*_invK;
    }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
