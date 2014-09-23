#include "STKS_LDACSchemeSys.hh"
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

MethodStrategyProvider<STKS_LDACSchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
stksLDAcSchemeSysProvider("STKS_SysLDAC");

//////////////////////////////////////////////////////////////////////////////

STKS_LDACSchemeSys::STKS_LDACSchemeSys(const std::string& name) :
  STKS_SplitterSys(name),
  m_sumKplus(),
  m_inv_K(),
  m_uTemp(),
  m_phiT()
{
}

//////////////////////////////////////////////////////////////////////////////

STKS_LDACSchemeSys::~STKS_LDACSchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKS_LDACSchemeSys::setup()
{
  STKS_SplitterSys::setup();
 

  // Do not need a dissipation contribution from past 
 getMethodData().getDistributionData().needDiss = false;

  m_uTemp.resize(_nbEquations);
  m_sumKplus.resize(_nbEquations, _nbEquations);
  m_inv_K.resize(_nbEquations, _nbEquations);
 m_phiT.resize(_nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void STKS_LDACSchemeSys::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void STKS_LDACSchemeSys::distribute(vector<RealVector>& residual)
{  
  const CFuint nbStatesInCell = _nbStatesInCell;
  m_phiT = getMethodData().getDistributionData().phi;

  m_sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    m_sumKplus  += *_kPlus[iState];
  }
  _inverter->invert(m_sumKplus, m_inv_K);
  m_uTemp = m_inv_K*m_phiT;
  
  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
    residual[iState] = *_kPlus[iState]*m_uTemp;
 
  }

  vector<RealMatrix>& betas = *getMethodData().getDistributionData().currBetaMat;
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
  {
    betas[iState] = (*_kPlus[iState])*m_inv_K;
  }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
