#include "FluctSplit/NLimSchemeCScalar.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitScalar.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NLimSchemeCScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
nlimcSchemeScalarProvider("ScalarNLimC");

//////////////////////////////////////////////////////////////////////////////

NLimSchemeCScalar::NLimSchemeCScalar(const std::string& name) :
  RDS_SplitterScalar(name),
  m_sumKplusU(),
  m_sumKplus(),
  m_uTemp(),
  m_uMin(),
  m_temp(),
  m_phiN(),
  m_betaN(),
  m_betaNPlus(),
  m_sumBetaNPlus(),
  m_beta()
{
}

//////////////////////////////////////////////////////////////////////////////

NLimSchemeCScalar::~NLimSchemeCScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void NLimSchemeCScalar::setup()
{
  RDS_SplitterScalar::setup();

  m_sumKplusU.resize(_nbEquations);
  m_sumKplus.resize(_nbEquations);
  m_uTemp.resize(_nbEquations);
  m_uMin.resize(_nbEquations);
  m_temp.resize(_nbEquations);
  m_sumBetaNPlus.resize(_nbEquations);

  const CFuint maxNbStatesInCell = _kPlus.size();

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

void NLimSchemeCScalar::distribute(vector<RealVector>& residual)
{
  const RealVector& phiT = getMethodData().getDistributionData().phi;
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  
  m_sumKplusU = _kPlus[0] * (*tStates[0]);
  m_sumKplus  = _kPlus[0];

  const CFuint nbStatesInCell = _nbStatesInCell;
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    m_sumKplusU += _kPlus[iState] * (*tStates[iState]);
    m_sumKplus  += _kPlus[iState];
  }
  m_temp = ( m_sumKplusU - phiT );
  m_uTemp = m_temp / m_sumKplus;

  m_sumBetaNPlus = 0.0;
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    m_uMin = *tStates[iState] - m_uTemp;
    m_phiN[iState] = _kPlus[iState]*m_uMin;
    m_betaN[iState] = (m_phiN[iState] + MathTools::MathConsts::CFrealEps()) / (phiT + nbStatesInCell*MathTools::MathConsts::CFrealEps());

    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      m_betaNPlus[iState][iEq] = std::max(m_betaN[iState][iEq],MathTools::MathConsts::CFrealEps());
      m_sumBetaNPlus[iEq] += m_betaNPlus[iState][iEq];
    }
  }

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    m_beta[iState] = m_betaNPlus[iState] / m_sumBetaNPlus;
    residual[iState] = m_beta[iState] * phiT;
  }
}

//////////////////////////////////////////////////////////////////////////////

void NLimSchemeCScalar::distributePart(vector<RealVector>& residual)
{
  throw Common::NotImplementedException (FromHere(),"NLimSchemeCScalar::distributePart()");
}

//////////////////////////////////////////////////////////////////////////////

void NLimSchemeCScalar::computePicardJacob(vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"NLimSchemeCScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
