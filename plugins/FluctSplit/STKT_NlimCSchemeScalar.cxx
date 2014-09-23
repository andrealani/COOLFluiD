#include "STKT_NlimCSchemeScalar.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STKT_NlimCSchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
STKTNlimcSchemeScalarProvider("STKT_ScalarNlimC");

//////////////////////////////////////////////////////////////////////////////

STKT_NlimCSchemeScalar::STKT_NlimCSchemeScalar(const std::string& name) :
  STKT_SplitterScalar(name),
  m_sumKplus(),
  m_uTemp(),
  m_sumKplusU(),
  m_diss(),
  m_phiT()
{
}

//////////////////////////////////////////////////////////////////////////////

STKT_NlimCSchemeScalar::~STKT_NlimCSchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NlimCSchemeScalar::setup()
{
  STKT_SplitterScalar::setup();
   // Do not need a dissipation contribution from past
 getMethodData().getDistributionData().needDiss = false;

  m_sumKplus.resize(_nbEquations);
 m_uTemp.resize(_nbEquations);
   m_sumKplusU.resize(_nbEquations);
 m_diss.resize(_nbEquations);
 m_phiT.resize(_nbEquations);
const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
m_sumBetaNPlus.resize(_nbEquations);
m_phitot.resize(_nbEquations);
  m_phiN.resize(maxNbStatesInCell);
  m_betaN.resize(maxNbStatesInCell);
  m_betaNPlus.resize(maxNbStatesInCell);
  m_beta.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    m_phiN[iState].resize(_nbEquations);
    m_betaN[iState].resize(_nbEquations);
    m_betaNPlus[iState].resize(_nbEquations);
    m_beta[iState].resize(_nbEquations);
  }
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NlimCSchemeScalar::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NlimCSchemeScalar::distribute(vector<RealVector>& residual)
{
   ///@todo No mesh deformation implemented here!!
  const CFuint nbStatesInCell = _nbStatesInCell;

  DistributionData& ddata = getMethodData().getDistributionData();
  const RealVector& m_phiT = ddata.phi;

 const vector<State*>& tStates = *ddata.tStates;

  m_sumKplus = _kPlus[0];
  m_sumKplusU = _kPlus[0]*(*tStates[0]);

  for (CFuint iState = 1; iState < nbStatesInCell; ++ iState){
    m_sumKplus += _kPlus[iState];
     m_sumKplusU += _kPlus[iState]*(*tStates[iState]);
  }

  m_phitot = 0.0;
  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){
    // We first add the Residual of LDA
    m_uTemp = _kPlus[iState]/m_sumKplus;
   m_phiN[iState] = m_uTemp*m_phiT;
  // And then the dissipation
  m_diss = _kPlus[iState]*(*tStates[iState])- (_kPlus[iState]*m_sumKplusU/m_sumKplus);

  m_phiN[iState] += m_diss;

  m_phitot += m_phiN[iState];
}

  //limitation of he N scheme
 m_sumBetaNPlus = 0.0;
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {

    m_betaN[iState] = (m_phiN[iState] + MathTools::MathConsts::CFrealEps()) / (m_phitot + nbStatesInCell*MathTools::MathConsts::CFrealEps());

    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      m_betaNPlus[iState][iEq] = std::max(m_betaN[iState][iEq],MathTools::MathConsts::CFrealEps());
      m_sumBetaNPlus[iEq] += m_betaNPlus[iState][iEq];
    }
  }

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    m_beta[iState] = m_betaNPlus[iState] / m_sumBetaNPlus;
    residual[iState] = m_beta[iState] *  m_phitot;
  }

      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations
    vector<RealMatrix>& betas =  *ddata.currBetaMat;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      betas[iState] = _kPlus[iState]/m_sumKplus;
    }

}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
