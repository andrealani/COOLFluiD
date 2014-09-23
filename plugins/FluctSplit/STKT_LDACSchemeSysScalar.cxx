#include "STKT_LDACSchemeSysScalar.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "MathTools/MatrixInverter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STKT_LDACSchemeSysScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeRiccLDAcSchemeSysScalarProvider("STKT_SysScalarLDAC");

//////////////////////////////////////////////////////////////////////////////

STKT_LDACSchemeSysScalar::STKT_LDACSchemeSysScalar(const std::string& name) :
  STKT_SplitterSysScalar(name),
  m_sumKplus(),
  m_invK(), 
  m_k(),
  m_beta(),
  m_invSumKplusScalar(),
  m_uTemp(),
  m_betaScalar()  
{
}

//////////////////////////////////////////////////////////////////////////////

STKT_LDACSchemeSysScalar::~STKT_LDACSchemeSysScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKT_LDACSchemeSysScalar::setup()
{
  STKT_SplitterSysScalar::setup();
 
  // Do not need a dissipation contribution from past 
  getMethodData().getDistributionData().needDiss = false;
  m_sumKplus.resize(_sysSize, _sysSize);
  m_invK.resize(_sysSize, _sysSize);
  m_k.resize(_sysSize, _sysSize);
  m_beta.resize(_sysSize, _sysSize);
  
  m_invSumKplusScalar.resize(_sysSize);
  m_uTemp.resize(_sysSize);
  m_betaScalar.resize(_scalarEqIDs.size());
}
      
//////////////////////////////////////////////////////////////////////////////

void STKT_LDACSchemeSysScalar::computePicardJacob(std::vector<RealMatrix*>& jacob)
{  
  const CFuint scalarSize = _scalarEqIDs.size();
  RealSliceMatrix::setNbRowsCols(_sysSize, _sysSize);
  
  // carefull with the signs !!!
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    // beta coefficient
    m_beta = (*_kPlus[iState])*m_invK;
    m_betaScalar = _kPlusScalar[iState]*m_invSumKplusScalar;
    
    const CFuint nStart = iState*_nbStatesInCell;
    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      m_k = (*_kPlus[jState]) + (*_kMin[jState]);
      RealMatrix& jacobIJ = (*jacob[nStart + jState]);
      jacobIJ.slice(_sysStartID,_sysStartID) = m_beta*m_k;
      
      for (CFuint i = 0; i < scalarSize; ++i) {
	const CFuint eqID = _scalarEqIDs[i];
	jacobIJ(eqID,eqID) = m_betaScalar[i]*_kScalar[jState][i];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void STKT_LDACSchemeSysScalar::distribute(vector<RealVector>& residual)
{
  RealSliceVector::setSize(_sysSize);
  
  const CFuint scalarSize = _scalarEqIDs.size();
  for (CFuint i = 0; i < scalarSize; ++i) {
    m_invSumKplusScalar[i] = _kPlusScalar[0][i];
  }
  
  m_sumKplus = *_kPlus[0]; 
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    m_sumKplus  += *_kPlus[iState];
    for (CFuint i = 0; i < scalarSize; ++i) {
      m_invSumKplusScalar[i] += _kPlusScalar[iState][i];
    }
  }
  
  _inverter->invert(m_sumKplus, m_invK);
  
  const RealVector& phi = getMethodData().getDistributionData().phi;
  m_uTemp = m_invK*phi.slice(_sysStartID);
  m_invSumKplusScalar = 1./m_invSumKplusScalar;
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) { 
    RealVector& res = residual[iState];
    res.slice(_sysStartID) = (*_kPlus[iState])*m_uTemp;
    
    for (CFuint i = 0; i < scalarSize; ++i) {
      const CFuint eqID = _scalarEqIDs[i];
      res[eqID] = _kPlusScalar[iState][i]*phi[eqID]*m_invSumKplusScalar[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
