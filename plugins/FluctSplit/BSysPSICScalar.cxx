#include "BSysPSICScalar.hh"
#include "MathTools/MatrixInverter.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitSystem.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/MeshData.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<BSysPSICScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
BSysPSICScalarProvider("BSysPSICScalar");

//////////////////////////////////////////////////////////////////////////////

BSysPSICScalar::BSysPSICScalar(const std::string& name) :
  BSchemeBase<RDS_SplitterSysScalar>(name),
  _sumKplusU(),
  _sumKplus(),
  _invK(),
  _uInflow(),
  _uDiff(),
  _phiLDA(),
  _ut(),
  _phiN(),
  m_sumKplusU(),
  m_sumKplus(),
  m_sumBeta(),
  m_invCoeff()
{ 
}
      
//////////////////////////////////////////////////////////////////////////////

BSysPSICScalar::~BSysPSICScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void BSysPSICScalar::setup()
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
    _ut.resize(_sysSize);
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
    m_sumBeta.resize(scalarSize);  
    m_invCoeff.resize(scalarSize);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void BSysPSICScalar::distribute(vector<RealVector>& residual)
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
  cf_assert(_ut.size() == _sysSize);
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
    _ut = _invK * phi.slice(_sysStartID, _sysSize);
  }
  
  m_sumBeta = 0.0;
    
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
      RealVector& res =  residual[iState];
      const CFuint eqID = _scalarEqIDs[i];
      phiN[eqID] = _kPlusScalar[iState][i]*(ts[eqID] - (m_sumKplusU[i] - phi[eqID])/m_sumKplus[i]);
      res[eqID] = (MathChecks::isNotZero(phi[eqID])) ? phiN[eqID]/phi[eqID] : 0.0;
      
      // fill the diagonal part of the beta matrix (beta of scalar LDA)
      beta(eqID,eqID) = _kPlusScalar[iState][i]/m_sumKplus[i];
      
      m_sumBeta[i] += std::max<CFreal>(0.,res[eqID]);
    }
  }
    
  if (_sysSize > 0) {
    computeBlendingCoeff();
  }
  
  // inverse of useful coefficients
  for (CFuint i = 0; i < scalarSize; ++i) {
    const CFuint eqID = _scalarEqIDs[i];
    m_invCoeff[i] = (MathChecks::isNotZero(m_sumBeta[i])) ? phi[eqID]/m_sumBeta[i] : 0.0;
  }
  
  if (_sysSize > 0) {
    if ( m_store_thetas ) storeThetas();
  }
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    if (_sysSize > 0) {
      // B system part
      _phiLDA = (*_kPlus[iState])*_ut;
      const CFuint endSys = _sysStartID + _sysSize;
      for (CFuint iEq = _sysStartID, jEq = 0; iEq < endSys; ++iEq, ++jEq) {
	const CFreal theta = m_theta[jEq];
	residual[iState][iEq] = theta*_phiN[iState][iEq] +  (1. - theta)*_phiLDA[jEq];
      }
    }
    
    // PSI scalar part
    for (CFuint i = 0; i < scalarSize; ++i) {
      const CFuint eqID = _scalarEqIDs[i]; 
      residual[iState][eqID] = std::max<CFreal>(0., residual[iState][eqID])*m_invCoeff[i];
    }
    
    //  cout << iState << " => ";
    // cout.precision(16);  cout.setf(ios::scientific,ios::floatfield); cout <<  residual[iState] << endl;
  }  
}
      
//////////////////////////////////////////////////////////////////////////////

void BSysPSICScalar::computePicardJacob(vector<RealMatrix*>& jacob)
{
  //  const CFuint nbEqs = _sysSize + _scalarEqIDs.size();
  
  //   // carefull with the signs !!!
  //   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
  //     const CFuint nStart = iState*_nbStatesInCell;
  //     for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
  //       RealMatrix *const block = jacob[nStart + jState];
  
  //       _tempMat = _invK*(*_kMin[jState]);
  //       if (iState == jState) {
  // 	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
  // 	  _tempMat(iEq,iEq) -= 1.0;
  // 	}
  //       }
  
  //       _tempMat *= -1.0;
  
  //       (*block) = (*_kPlus[iState])*_tempMat;
  //     }
  //   }
}

//////////////////////////////////////////////////////////////////////////////

void BSysPSICScalar::computeBlendingCoeff()
{
  cout << "BSysPSICScalar::computeBlendingCoeff() not implemented!!" << endl;
  abort();
 //  for (CFuint iEq = 0; iEq < _sysSize; ++iEq)
//   {
//     m_theta[iEq] = max ( std::abs(_phi[iEq])/max(MathTools::MathConsts::CFrealEps(), m_sumPhiN[iEq]) , m_min_theta );
//   }
}

//////////////////////////////////////////////////////////////////////////////

void BSysPSICScalar::addExtraDissipation()
{
  cout << "BSysPSICScalar::addExtraDissipation() not implemented!!" << endl;
  abort();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
