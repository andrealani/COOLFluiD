#include "NSchemeCSysScalar.hh"
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

MethodStrategyProvider<NSchemeCSysScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
ncSchemeSysScalarProvider("SysScalarNC");

//////////////////////////////////////////////////////////////////////////////

NSchemeCSysScalar::NSchemeCSysScalar(const std::string& name) :
  RDS_SplitterSysScalar(name),
  _sumKplusU(),
  _sumKplus(),
  _invK(),
  _tempMat(),
  m_sumKplusU(),
  m_sumKplus()
{ 
}
      
//////////////////////////////////////////////////////////////////////////////

NSchemeCSysScalar::~NSchemeCSysScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCSysScalar::setup()
{
  RDS_SplitterSysScalar::setup();
  
  cf_always_assert(_sysSize > 0);
  cf_always_assert(_scalarEqIDs.size() > 0);
    
  // system data
  _sumKplusU.resize(_sysSize);
  _sumKplus.resize(_sysSize, _sysSize);
  _invK.resize(_sysSize, _sysSize);
  _tempMat.resize(_sysSize, _sysSize);
  
  //scalar data
  const CFuint scalarSize = _scalarEqIDs.size();
  m_sumKplusU.resize(scalarSize);
  m_sumKplus.resize(scalarSize);
}
      
//////////////////////////////////////////////////////////////////////////////

void NSchemeCSysScalar::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  
  // system part 
  _sumKplusU.slice(0, _sysSize) = (*_kPlus[0]) * tStates[0]->slice(_sysStartID, _sysSize);
  _sumKplus  = *_kPlus[0]; 
  
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
    
    // system part 
    _sumKplusU.slice(0, _sysSize) += kpM * ts.slice(_sysStartID, _sysSize);
    _sumKplus += kpM;
    
    // scalar part
    for (CFuint i = 0; i < scalarSize; ++i) {
      const CFreal kp = _kPlusScalar[iState][i];
      m_sumKplusU[i] += kp * ts[_scalarEqIDs[i]];
      m_sumKplus[i] += kp;
    }
  }
  
  _inverter->invert(_sumKplus, _invK);
  RealVector& phi = const_cast<RealVector&>(getMethodData().getDistributionData().phi);
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    RealVector& res = residual[iState];
    State& ts = *tStates[iState];
    
    // system part
    res.slice(_sysStartID, _sysSize) = (*_kPlus[iState])*(ts.slice(_sysStartID, _sysSize) - _invK*(_sumKplusU.slice(0, _sysSize) - phi.slice(_sysStartID, _sysSize)));
    
    //scalar part
    for (CFuint i = 0; i < scalarSize; ++i) {
      const CFuint eqID = _scalarEqIDs[i];
      res[eqID] = _kPlusScalar[iState][i]*(ts[eqID] - (m_sumKplusU[i] - phi[eqID])/m_sumKplus[i]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCSysScalar::computePicardJacob(vector<RealMatrix*>& jacob)
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

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
