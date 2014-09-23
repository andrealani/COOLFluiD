#include "FluctSplit/NLimDSchemeCScalar.hh"
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

MethodStrategyProvider<NLimDSchemeCScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
nlimdcSchemeScalarProvider("ScalarNLimDC");

//////////////////////////////////////////////////////////////////////////////

NLimDSchemeCScalar::NLimDSchemeCScalar(const std::string& name) :
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
  m_beta(),
  m_theta(),
  m_ubar()
{
}

//////////////////////////////////////////////////////////////////////////////

NLimDSchemeCScalar::~NLimDSchemeCScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void NLimDSchemeCScalar::setup()
{
  RDS_SplitterScalar::setup();

  m_sumKplusU.resize(_nbEquations);
  m_sumKplus.resize(_nbEquations);
  m_uTemp.resize(_nbEquations);
  m_uMin.resize(_nbEquations);
  m_temp.resize(_nbEquations);
  m_sumBetaNPlus.resize(_nbEquations);
  m_theta.resize(_nbEquations);
  m_ubar.resize(_nbEquations);

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

void NLimDSchemeCScalar::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;
  
  m_sumKplusU = _kPlus[0] * (*tStates[0]);
  m_sumKplus  = _kPlus[0];
  m_ubar = *tStates[0];

  CFLogDebugMax( "tStates[0] = " << *tStates[0] << "\n");

  const CFreal h = 0.0004; // this should be got from the geometry
  
  const CFuint nbStatesInCell = _nbStatesInCell;
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    m_sumKplusU += _kPlus[iState] * (*tStates[iState]);
    m_sumKplus  += _kPlus[iState];
    m_ubar  += *tStates[iState];

    CFLogDebugMax( "tStates[" << iState << "] = " << *tStates[iState] << "\n");
  }

  m_ubar /= nbStatesInCell;
  m_temp = ( m_sumKplusU - phiT );
  m_uTemp = m_temp / m_sumKplus;

  CFLogDebugMax( "m_sumKplusU = " << m_sumKplusU << "\n");
  CFLogDebugMax( "m_sumKplus  = " << m_sumKplus << "\n");
  CFLogDebugMax( "m_uTemp     = " << m_uTemp << "\n");
  CFLogDebugMax( "phiT        = " << phiT << "\n\n");

  m_sumBetaNPlus = 0.0;
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    m_uMin = *tStates[iState] - m_uTemp;
    m_phiN[iState] = _kPlus[iState]*m_uMin;

CFLogDebugMax("phiN["<< iState <<"] = " << m_phiN[iState] <<"\n");

    m_betaN[iState] = (m_phiN[iState] + MathTools::MathConsts::CFrealEps()) / (phiT + nbStatesInCell*MathTools::MathConsts::CFrealEps());

CFLogDebugMax("betaN["<< iState <<"] = " << m_betaN[iState] <<"\n");

    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      m_betaNPlus[iState][iEq] = std::max(m_betaN[iState][iEq],MathTools::MathConsts::CFrealEps());
      m_sumBetaNPlus[iEq] += m_betaNPlus[iState][iEq];
    }
CFLogDebugMax("betaNPlus["<< iState <<"] = " << m_betaNPlus[iState] <<"\n");
  }

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    m_beta[iState] = m_betaNPlus[iState] / m_sumBetaNPlus;
    residual[iState] = m_beta[iState] * phiT;
CFLogDebugMax("residual["<< iState <<"] = " << residual[iState] <<"\n");
  }

  // compute the theta function
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    m_theta[iEq] = std::min(1.0,1.0/((std::abs(phiT[iEq])/(m_ubar[iEq]*h*h))+MathTools::MathConsts::CFrealEps()));
  }

  // fix for spurious oscilations by adding a artificial diffusion
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    m_temp = _k[iState] * phiT;
    m_temp /= h;
    m_temp *= m_theta;
    residual[iState] += m_temp;
CFLogDebugMax("fixed_residual["<< iState <<"] = " << residual[iState] <<"\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

      void NLimDSchemeCScalar::distributePart(vector<RealVector>& residual)
{
  throw Common::NotImplementedException (FromHere(),"NLimDSchemeCScalar::distributePart()");
}

//////////////////////////////////////////////////////////////////////////////

void NLimDSchemeCScalar::computePicardJacob(vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"NLimDSchemeCScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
