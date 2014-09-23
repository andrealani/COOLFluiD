#include "Common/NotImplementedException.hh"

#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/STKT_BCSchemeScalar.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STKT_BCSchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeRiccBcSchemeScalarProvider("STKT_ScalarBC");

//////////////////////////////////////////////////////////////////////////////

STKT_BCSchemeScalar::STKT_BCSchemeScalar(const std::string& name) :
  STKT_SplitterScalar(name),
  m_sumKplus(),
  m_uTemp(),
  m_sumKplusU(),
  m_diss(),
  m_phiN(),
  m_phitot(),
  m_absphitot(),
  m_phiT()
{
}

//////////////////////////////////////////////////////////////////////////////

STKT_BCSchemeScalar::~STKT_BCSchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKT_BCSchemeScalar::setup()
{
  STKT_SplitterScalar::setup();
  // Do not need a dissipation contribution from past 
 getMethodData().getDistributionData().needDiss = false;
  CFuint maxNbStatesInCell = _kPlus.size();
  m_phiN.resize(maxNbStatesInCell);
  m_diss.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++ iState){
    m_phiN[iState].resize(_nbEquations);
    m_diss[iState].resize(_nbEquations);
 }
 
  m_sumKplus.resize(_nbEquations);
  m_uTemp.resize(_nbEquations);
  m_sumKplusU.resize(_nbEquations);
  m_phitot.resize(_nbEquations);
  m_absphitot.resize(_nbEquations);
  m_phiT.resize(_nbEquations);

 getMethodData().getDistributionData().needDiss = false;
}

//////////////////////////////////////////////////////////////////////////////

void STKT_BCSchemeScalar::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void STKT_BCSchemeScalar::distribute(vector<RealVector>& residual)
{  
   ///@todo No mesh deformation implemented here!!
  const CFuint nbStatesInCell = _nbStatesInCell;
  DistributionData& ddata = getMethodData().getDistributionData();
  const RealVector& m_phiT = ddata.phi;
  const vector<State*>& tStates = *ddata.tStates;
  CFreal theta;
  
  m_sumKplus = _kPlus[0];
  m_sumKplusU = _kPlus[0]*(*tStates[0]);

  for (CFuint iState = 1; iState < nbStatesInCell; ++ iState){
    m_sumKplus += _kPlus[iState];
    m_sumKplusU += _kPlus[iState]*(*tStates[iState]);
  }

  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
    // We first add the Residual of LDA
    m_uTemp = _kPlus[iState]/m_sumKplus;
    m_phiN[iState] = m_uTemp*m_phiT;
    residual[iState] = m_uTemp*m_phiT;
    // And then the dissipation
    m_diss[iState] = _kPlus[iState]*(*tStates[iState])- (_kPlus[iState]*m_sumKplusU/m_sumKplus);

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
  vector<RealMatrix>& betas = 
      *getMethodData().getDistributionData().currBetaMat;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations
      
      betas[iState] = _kPlus[iState]/m_sumKplus;
    }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
