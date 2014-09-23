#include "STKT_NCSchemeSysquad.hh"
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

MethodStrategyProvider<STKT_NCSchemeSysquad,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeRiccNcquadSchemeSysProvider("STKT_SysNCquad");

//////////////////////////////////////////////////////////////////////////////

STKT_NCSchemeSysquad::STKT_NCSchemeSysquad(const std::string& name) :
  STKT_SplitterSyshyb(name),
  m_sumKplus(),
  _invK(),
  m_uTemp(),
  m_diss(),
  m_sumKplusU(),
  m_phiT()
{
}

//////////////////////////////////////////////////////////////////////////////

STKT_NCSchemeSysquad::~STKT_NCSchemeSysquad()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NCSchemeSysquad::setup()
{
  STKT_SplitterSyshyb::setup();
 
  // Do not need a dissipation contribution from past 
 getMethodData().getDistributionData().needDiss = false;

  m_sumKplus.resize(_nbEquations, _nbEquations);
  m_uTemp.resize(_nbEquations);
  m_sumKplusU.resize(_nbEquations);
  m_diss.resize(_nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
   m_phiT.resize(_nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NCSchemeSysquad::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NCSchemeSysquad::distribute(vector<RealVector>& residual)
{  
   DistributionData& ddata = getMethodData().getDistributionData();
  const RealVector& m_phiT = ddata.phi;
  const CFuint nbStatesInCell = _nbStatesInCell;
 const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;

  m_sumKplus = *_kPlus[0];
  m_sumKplusU = *_kPlus[0]*(*tStates[0]);

  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    m_sumKplus  += *_kPlus[iState];
    m_sumKplusU += *_kPlus[iState]*(*tStates[iState]);
  }
  _inverter->invert(m_sumKplus, _invK);
  m_uTemp = _invK*m_phiT;
  
  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
    residual[iState] = *_kPlus[iState]*m_uTemp;
  }

  m_uTemp = _invK*m_sumKplusU;
 
  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
    // And then the dissipation
    m_diss = *_kPlus[iState]*(*tStates[iState]) - *_kPlus[iState]*m_uTemp;

    residual[iState] += m_diss;
  }


    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations

      RealMatrix& betaLDA = (*getMethodData().getDistributionData().currBetaMat)[iState];
      betaLDA = (*_kPlus[iState])*_invK;
    }

}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
