#include "BSchemeCSysScalar.hh"
#include "MathTools/MatrixInverter.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitSystem.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<BSchemeCSysScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
BSchemeCSysScalarProvider("BSchemeCSysScalar");

//////////////////////////////////////////////////////////////////////////////

BSchemeCSysScalar::BSchemeCSysScalar(const std::string& name) :
  BSchemeBase<RDS_SplitterSysScalar>(name),
  _sumKplusU(),
  _sumKplus(),
  _invK(),
  _uInflow(),
  _uDiff(),
  _phiLDA(),
  _phiN(),
  m_sumKplusU(),
  m_sumKplus(),
  m_invCoeff(),
  m_uInflow(),
  m_uTempSc()
{ 
}
      
//////////////////////////////////////////////////////////////////////////////

BSchemeCSysScalar::~BSchemeCSysScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeCSysScalar::setup()
{
  BSchemeBase<RDS_SplitterSysScalar>::setup();
  
  cf_always_assert(_sysSize > 0 || _scalarEqIDs.size() > 0);
  cf_always_assert(this->m_theta.size() >= _sysSize);
  
  if (_sysSize > 0) {
    // system data
    _sumKplusU.resize(_sysSize);
    _sumKplus.resize(_sysSize, _sysSize);
    _invK.resize(_sysSize, _sysSize);
    _uInflow.resize(_sysSize);
    _uDiff.resize(_sysSize);
    _phiLDA.resize(_sysSize); 
  }
  
  const CFuint scalarSize = _scalarEqIDs.size();
  _phiN.resize(m_maxNbStatesInCell);
  for (CFuint i = 0; i< m_maxNbStatesInCell; ++i) {
    _phiN[i].resize(_sysSize + scalarSize);
  }
  
  if (scalarSize > 0) {
    //scalar data
    m_sumKplusU.resize(scalarSize);
    m_sumKplus.resize(scalarSize);
    m_invCoeff.resize(scalarSize);
    m_uInflow.resize(scalarSize);
    m_uTempSc.resize(scalarSize);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void BSchemeCSysScalar::distribute(vector<RealVector>& residual)
{
  DistributionData& distdata = getMethodData().getDistributionData();
  const vector<State*>& tStates = *distdata.tStates;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  RealVector& phi = distdata.phi;
  
  // full (system + scalar) size
  cf_assert(residual[0].size() == nbEqs);
  cf_assert(phi.size() == nbEqs);
  
  // system size
  cf_assert(_sumKplus.size() ==_sysSize*_sysSize);
  cf_assert(_sumKplusU.size() == _sysSize);
  cf_assert(_kPlus[0]->size() == _sysSize*_sysSize);
  cf_assert(_invK.size() == _sysSize*_sysSize);
  cf_assert(m_uTemp.size() == _sysSize);
  cf_assert(m_theta.size() >= _sysSize);
  
  // scalar size
  cf_assert(m_sumKplus.size() == _scalarEqIDs.size());
  cf_assert(m_sumKplusU.size() == _scalarEqIDs.size());
  cf_assert(_kPlusScalar[0].size() == _scalarEqIDs.size());

  // system part 
  if (_sysSize > 0) {
    _sumKplusU = (*_kPlus[0]) * tStates[0]->slice(_firstVarID, _sysSize);
    _sumKplus  = *_kPlus[0];
  }
  
  // scalar part
  const CFuint scalarSize = _scalarEqIDs.size();
  for (CFuint i = 0; i < scalarSize; ++i) {
    const CFreal kp = _kPlusScalar[0][i];
    m_sumKplusU[i] = kp * (*tStates[0])[_scalarEqIDs[i]];
    m_sumKplus[i] = kp;
  }
  
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    State& ts = *tStates[iState];
    RealMatrix& kpM = (*_kPlus[iState]);

    if (_sysSize > 0) {
      // system part 
      _sumKplusU += kpM * ts.slice(_sysStartID, _sysSize);
      _sumKplus += kpM;
    }
    
    // scalar part
    for (CFuint i = 0; i < scalarSize; ++i) {
      const CFreal kp = _kPlusScalar[iState][i];
      m_sumKplusU[i] += kp * ts[_scalarEqIDs[i]];
      m_sumKplus[i] += kp;
    }
  }
  
  if (_sysSize > 0) {
    _inverter->invert(_sumKplus, _invK);
    _uInflow = _invK * (_sumKplusU - phi.slice(_sysStartID, _sysSize));
    m_uTemp = _invK * phi.slice(_sysStartID, _sysSize);
  }
  
  for (CFuint i = 0; i < scalarSize; ++i) {
    const CFuint eqID = _scalarEqIDs[i];
    m_invCoeff[i] = 1./m_sumKplus[i];
    m_uInflow[i] = m_invCoeff[i]*(m_sumKplusU[i] - phi[eqID]);
    m_uTempSc[i] = m_invCoeff[i]*phi[eqID];
  }
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    RealVector& phiN = _phiN[iState];
    cf_assert(phiN.size() == nbEqs);
    State& ts = *tStates[iState];
    RealMatrix& beta = (*distdata.currBetaMat)[iState];
    
    if (_sysSize > 0) {
      beta = 0.0;
      beta.slice(_sysStartID, _sysStartID, _sysSize, _sysSize) = (*_kPlus[iState])*_invK;
      // system part
      phiN.slice(_sysStartID, _sysSize) = (*_kPlus[iState])*(ts.slice(_sysStartID, _sysSize) - _uInflow);
    }
    
    // scalar part
    for (CFuint i = 0; i < scalarSize; ++i) {
      const CFuint eqID = _scalarEqIDs[i];
      phiN[eqID] = _kPlusScalar[iState][i]*(ts[eqID] - m_uInflow[i]);
      
      // fill the diagonal part of the beta matrix (beta of scalar LDA)
      beta(eqID,eqID) = _kPlusScalar[iState][i]/m_sumKplus[i];
    }
  }
  
  computeBlendingCoeff();
  
  if ( m_store_thetas ) storeThetas();
  
  assert(m_theta.size() == _sysSize + scalarSize);
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    if (_sysSize > 0) {
      // B system part
      _phiLDA = (*_kPlus[iState])*m_uTemp;
      const CFuint endSys = _sysStartID + _sysSize;
      for (CFuint iEq = _sysStartID, jEq = 0; iEq < endSys; ++iEq, ++jEq) {
	const CFreal theta = m_theta[jEq];
	residual[iState][iEq] = theta*_phiN[iState][iEq] +  (1. - theta)*_phiLDA[jEq];
      }
    }
    
    // B scalar part
    for (CFuint i = 0; i < scalarSize; ++i) {
      const CFuint eqID = _scalarEqIDs[i]; 
      const CFreal theta = m_theta[eqID];
      const CFreal phiLDA = _kPlusScalar[iState][i]*m_uTempSc[i];
      residual[iState][eqID] = theta*_phiN[iState][eqID] +  (1. - theta)*phiLDA;
    }
    
    // cout << iState << " => ";
    // cout.precision(16);  cout.setf(ios::scientific,ios::floatfield); cout <<  residual[iState] << endl;
  }  
}
      
//////////////////////////////////////////////////////////////////////////////

void BSchemeCSysScalar::computePicardJacob(vector<RealMatrix*>& jacob)
{
  cout << "BSchemeCSysScalar::computePicardJacob() not implemented!!" << endl;
  abort();
}
      
//////////////////////////////////////////////////////////////////////////////
      
void BSchemeCSysScalar::computeBlendingCoeff()
{
  cout << "BSchemeCSysScalar::computeBlendingCoeff() not implemented!!" << endl;
  abort();
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeCSysScalar::addExtraDissipation()
{
  cout << "BSchemeCSysScalar::addExtraDissipation() not implemented!!" << endl;
  abort();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
