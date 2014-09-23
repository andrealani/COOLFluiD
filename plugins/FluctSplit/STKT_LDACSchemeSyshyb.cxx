#include "STKT_LDACSchemeSyshyb.hh"
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

MethodStrategyProvider<STKT_LDACSchemeSyshyb,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeRiccLDAcSchemehybSysProvider("STKT_SysLDAC_hyb");

//////////////////////////////////////////////////////////////////////////////

STKT_LDACSchemeSyshyb::STKT_LDACSchemeSyshyb(const std::string& name) :
  STKT_SplitterSyshyb(name),
  m_sumKplus(),
  _invK(),
  m_uTemp(),
  m_phiT()
{
}

//////////////////////////////////////////////////////////////////////////////

STKT_LDACSchemeSyshyb::~STKT_LDACSchemeSyshyb()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKT_LDACSchemeSyshyb::setup()
{
  STKT_SplitterSyshyb::setup();
 
  // Do not need a dissipation contribution from past 
 getMethodData().getDistributionData().needDiss = false;
  m_uTemp.resize(_nbEquations);
  m_sumKplus.resize(_nbEquations, _nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
  m_phiT.resize(_nbEquations);
  _k.resize(_nbEquations, _nbEquations);
  _beta.resize(_nbEquations, _nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void STKT_LDACSchemeSyshyb::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  // carefull with the signs !!!
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    // beta coefficient
    _beta = (*_kPlus[iState])*_invK;
    const CFuint nStart = iState*_nbStatesInCell;
    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      _k = (*_kPlus[jState]) + (*_kMin[jState]);
      (*jacob[nStart + jState]) = _beta*_k;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void STKT_LDACSchemeSyshyb::distribute(vector<RealVector>& residual)
{  
  m_sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    m_sumKplus  += *_kPlus[iState];
  }
  
  _inverter->invert(m_sumKplus, _invK);
  m_uTemp = _invK*getMethodData().getDistributionData().phi;
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++ iState){ 
    residual[iState] = *_kPlus[iState]*m_uTemp;
  }
  
//     for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
//       // AL: compute betas on the fly is more expensive because
//       // it involves a matrix*matrix and a matrix*vector
//       // while normal scheme requires only two matrix*vector operations
// 
//       RealMatrix& betaLDA = (*getMethodData().getDistributionData().currBetaMat)[iState];
//       betaLDA = (*_kPlus[iState])*_invK;
//     }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
