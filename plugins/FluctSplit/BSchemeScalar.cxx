#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitScalar.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/BSchemeScalar.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////
      
MethodStrategyProvider<BSchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
aScalarBProvider("ScalarB");

//////////////////////////////////////////////////////////////////////////////

void BSchemeScalar::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

BSchemeScalar::BSchemeScalar(const std::string& name) :
  BSchemeBase<RDS_SplitterScalar>(name),
  _sumKminU(),
  _sumKmin(),
  _sumKplus(),
  _uInflow(),
  _uDiff(),
  _phi(),
  _phiLDA(),
  _temp()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

BSchemeScalar::~BSchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeScalar::configure ( Config::ConfigArgs& args )
{
  BSchemeBase<RDS_SplitterScalar>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeScalar::distribute(vector<RealVector>& residual)
{
  DistributionData& distdata = getMethodData().getDistributionData();
  const vector<State*>& tStates = *distdata.tStates;
  const CFuint nbEqs = _nbEquations;
  
  _phi = _k[0]*(*tStates[0]);
  _sumKminU = _kMin[0]*(*tStates[0]);
  _sumKmin = _kMin[0];
  _sumKplus = _kPlus[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _phi += _k[iState]*(*tStates[iState]);
    _sumKminU += _kMin[iState]*(*tStates[iState]);
    _sumKmin  += _kMin[iState];
    _sumKplus += _kPlus[iState];
  }
  _uInflow = _sumKminU / _sumKmin;
  m_uTemp = _phi / _sumKplus;

  m_sumPhiN = 0.0;
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uDiff = *tStates[iState] - _uInflow;
    m_phiN[iState] = _kPlus[iState]*_uDiff;
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      m_sumPhiN[iEq] += std::abs(m_phiN[iState][iEq]);
    }
  }

  computeBlendingCoeff();
  
  if ( m_store_thetas ) storeThetas();

  // blend the residual
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _phiLDA = _kPlus[iState]*m_uTemp;
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      residual[iState][iEq] = m_theta[iEq]*m_phiN[iState][iEq] + (1. - m_theta[iEq])*_phiLDA[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeScalar::distributePart(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  
  _phi.slice(0, _nbEquations) = _k[0].slice(0, _nbEquations) * tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKminU.slice(0, _nbEquations) = _kMin[0].slice(0, _nbEquations) * tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKmin = _kMin[0];
  _sumKplus = _kPlus[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _temp.slice(0,_nbEquations) = _k[iState].slice(0, _nbEquations) * tStates[iState]->slice(_firstVarID, _nbEquations);
    _phi += _temp;
    _temp.slice(0, _nbEquations) = _kMin[iState].slice(0, _nbEquations) * tStates[iState]->slice(_firstVarID, _nbEquations);
    _sumKminU += _temp;
    _sumKmin  += _kMin[iState];
    _sumKplus += _kPlus[iState];
  }

  _uInflow = _sumKminU/_sumKmin;
  m_uTemp = _phi/_sumKplus;

  m_sumPhiN = 0.0;
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uDiff.slice(0, _nbEquations) = tStates[iState]->slice(_firstVarID, _nbEquations) - _uInflow.slice(0, _nbEquations);
    m_phiN[iState] = _kPlus[iState]*_uDiff;
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      m_sumPhiN[iEq] += std::abs(m_phiN[iState][iEq]);
    }
  }

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    m_theta[iEq] = std::abs(_phi[iEq])/max(MathTools::MathConsts::CFrealEps(), m_sumPhiN[iEq]);
  }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _phiLDA = _kPlus[iState]*m_uTemp;
    for (CFuint iEq = _firstVarID, jEq = 0; iEq < _lastVarID; ++iEq, ++jEq) {
      residual[iState][iEq] = m_theta[jEq]*m_phiN[iState][jEq] +
        (1. - m_theta[jEq])*_phiLDA[jEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeScalar::setup()
{
  BSchemeBase<RDS_SplitterScalar>::setup();

  _sumKminU.resize(_nbEquations);
  _sumKmin.resize(_nbEquations);
  _sumKplus.resize(_nbEquations);
  _uInflow.resize(_nbEquations);
  _uDiff.resize(_nbEquations);
  _phi.resize(_nbEquations);
  _phiLDA.resize(_nbEquations);
  _temp.resize(_nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeScalar::computeBlendingCoeff()
{
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq)
  {
    m_theta[iEq] = max ( std::abs(_phi[iEq])/max(MathTools::MathConsts::CFrealEps(), m_sumPhiN[iEq]) , m_min_theta );
  }
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeScalar::addExtraDissipation()
{
  throw Common::NotImplementedException (FromHere(),"BSchemeScalar::addExtraDissipation()");
}

//////////////////////////////////////////////////////////////////////////////
     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
