#include "MathTools/MatrixInverter.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"

#include "FluctSplit/FluctSplitSystem.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/BSchemeSys.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<BSchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
aBSchemeSysProvider("SysB");

//////////////////////////////////////////////////////////////////////////////

BSchemeSys::BSchemeSys(const std::string& name) :
  BSchemeBase<RDS_SplitterSys>(name),
  _sumKminU(),
  _sumKmin(),
  _sumKplus(),
  _k(),
  _invKmin(),
  _invKplus(),
  _uInflow(),
  _phi(),
  _phiLDA(),
  _temp(),
  _betaLDA(),
  _tempMat()
 {
 }

//////////////////////////////////////////////////////////////////////////////

BSchemeSys::~BSchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeSys::setup()
{
  BSchemeBase<RDS_SplitterSys>::setup();

  _sumKminU.resize(_nbEquations);
  _sumKmin.resize(_nbEquations, _nbEquations);
  _sumKplus.resize(_nbEquations, _nbEquations);
  _k.resize(_nbEquations, _nbEquations);
  _invKmin.resize(_nbEquations, _nbEquations);
  _invKplus.resize(_nbEquations, _nbEquations);
  _uInflow.resize(_nbEquations);
  _phi.resize(_nbEquations);
  _phiLDA.resize(_nbEquations);
  _temp.resize(_nbEquations);
  _betaLDA.resize(_nbEquations, _nbEquations);
  _tempMat.resize(_nbEquations, _nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeSys::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  
  _k = (*_kPlus[0]) + (*_kMin[0]);
  _phi = _k * (*tStates[0]);
  _sumKminU = (*_kMin[0]) * (*tStates[0]);
  _sumKmin  = *_kMin[0];
  _sumKplus = *_kPlus[0];

  const CFuint nbStatesInCell = _nbStatesInCell;
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _k = (*_kPlus[iState]) + (*_kMin[iState]);
    _phi += _k*(*tStates[iState]);
    _sumKminU += (*_kMin[iState])*(*tStates[iState]);
    _sumKmin  += *_kMin[iState];
    _sumKplus += *_kPlus[iState];
  }

  _inverter->invert(_sumKmin, _invKmin);
  _uInflow = _invKmin*_sumKminU;

  _inverter->invert(_sumKplus, _invKplus);
  m_uTemp = _invKplus*_phi;

  m_sumPhiN = 0.0;
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    m_phiN[iState] = (*_kPlus[iState])*(*tStates[iState] - _uInflow);
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      m_sumPhiN[iEq] += std::abs(m_phiN[iState][iEq]);
    }
  }
  
  if (m_firstOrderJacob && getMethodData().getDistributionData().isPerturb) {
    m_theta = 1.0;
  }
  else { 
    computeBlendingCoeff();
  }

  if ( m_store_thetas ) storeThetas();
    
  if (!getMethodData().getDistributionData().computeBetas) {
    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
      _phiLDA = (*_kPlus[iState])*m_uTemp;
      for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
	residual[iState][iEq] = m_theta[iEq]*m_phiN[iState][iEq] +
	  (1. - m_theta[iEq])*_phiLDA[iEq];
      }
    }
  }
  else {
    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
      RealMatrix& currBeta = (*getMethodData().getDistributionData().currBetaMat)[iState];
      currBeta = (*_kPlus[iState])*_invKplus;
      _phiLDA = currBeta*_phi;
      for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
	residual[iState][iEq] = m_theta[iEq]*m_phiN[iState][iEq] +
	  (1. - m_theta[iEq])*_phiLDA[iEq];
      }
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

void BSchemeSys::distributePart(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  
  _k = (*_kPlus[0]) + (*_kMin[0]);
  _phi.slice(0, _nbEquations) = _k * tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKminU.slice(0, _nbEquations) = (*_kMin[0]) * tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKmin  = *_kMin[0];
  _sumKplus = *_kPlus[0];

  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _k = (*_kPlus[iState]) + (*_kMin[iState]);
    _phi.slice(0, _nbEquations) += _k * tStates[iState]->slice(_firstVarID, _nbEquations);
    _sumKminU.slice(0, _nbEquations) +=
      (*_kMin[iState]) * tStates[iState]->slice(_firstVarID, _nbEquations);
    _sumKmin  += *_kMin[iState];
    _sumKplus += *_kPlus[iState];
  }

  _inverter->invert(_sumKmin, _invKmin);
  _uInflow = _invKmin * _sumKminU;

  _inverter->invert(_sumKplus, _invKplus);
  m_uTemp = _invKplus * _phi;
  m_sumPhiN = 0.0;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    m_phiN[iState].slice(0, _nbEquations) = (*_kPlus[iState]) *
      (tStates[iState]->slice(_firstVarID, _nbEquations) - _uInflow.slice(0, _nbEquations));
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      m_sumPhiN[iEq] += std::abs(m_phiN[iState][iEq]);
    }
  }
   
  if (m_firstOrderJacob && getMethodData().getDistributionData().isPerturb) {
    // cout << "JACOB" << endl;
    m_theta = 1.0;
  }
  else { 
    // cout << "RHS" << endl;  
    computeBlendingCoeff();
  }
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _phiLDA = (*_kPlus[iState]) * m_uTemp;
    for (CFuint iEq = _firstVarID, jEq = 0; iEq < _lastVarID; ++iEq, ++jEq) {
      residual[iState][iEq] = m_theta[jEq]*m_phiN[iState][jEq] +
        (1. - m_theta[jEq])*_phiLDA[jEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeSys::computePicardJacob(vector<RealMatrix*>& jacob)
{
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    const CFuint nStart = iState*_nbStatesInCell;
    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      RealMatrix *const block = jacob[nStart + jState];

      _tempMat = _invKmin*(*_kMin[jState]);
      if (iState == jState) {
	for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
	  _tempMat(iEq,iEq) -= 1.0;
	}
      }

      _tempMat *= -1.0;

      (*block) = (*_kPlus[iState])*_tempMat;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeSys::computeBlendingCoeff()
{
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq)
  {
    m_theta[iEq] = max ( std::abs(_phi[iEq])/max(MathTools::MathConsts::CFrealEps(), m_sumPhiN[iEq]) , m_min_theta );
  }
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeSys::addExtraDissipation()
{
  throw Common::NotImplementedException (FromHere(),"BSchemeSys::addExtraDissipation()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
