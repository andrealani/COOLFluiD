#include "STKT_NCSchemeScalar.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STKT_NCSchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeRiccNcSchemeScalarProvider("STKT_ScalarNC");

//////////////////////////////////////////////////////////////////////////////

STKT_NCSchemeScalar::STKT_NCSchemeScalar(const std::string& name) :
  STKT_SplitterScalar(name),
  m_sumKplus(),
  m_uTemp(),
m_diss(),
  m_sumKplusU(),
  m_phiT()
{
}

//////////////////////////////////////////////////////////////////////////////

STKT_NCSchemeScalar::~STKT_NCSchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NCSchemeScalar::setup()
{
  STKT_SplitterScalar::setup();
   // Do not need a dissipation contribution from past 
 getMethodData().getDistributionData().needDiss = false;

  m_sumKplus.resize(_nbEquations);
 m_uTemp.resize(_nbEquations);
   m_sumKplusU.resize(_nbEquations);
 m_diss.resize(_nbEquations);
 m_phiT.resize(_nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NCSchemeScalar::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NCSchemeScalar::distribute(vector<RealVector>& residual)
{  
   ///@todo No mesh deformation implemented here!!
  const CFuint nbStatesInCell = _nbStatesInCell;
 DistributionData& ddata = getMethodData().getDistributionData();
  const RealVector& m_phiT = ddata.phi;

 const vector<State*>& tStates = *ddata.tStates;

  
 m_sumKplus = _kPlus[0];
 m_sumKplusU = _kPlus[0]*(*tStates[0]);

  for (CFuint iState = 1; iState < nbStatesInCell; ++ iState){
    m_sumKplus += _kPlus[iState];
     m_sumKplusU += _kPlus[iState]*(*tStates[iState]);
  }

  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
    // We first add the Residual of LDA
    m_uTemp = _kPlus[iState]/m_sumKplus;
   residual[iState] = m_uTemp*m_phiT; 
  // And then the dissipation
  m_diss = _kPlus[iState]*(*tStates[iState])- (_kPlus[iState]*m_sumKplusU/m_sumKplus);

  residual[iState] += m_diss;

}

  vector<RealMatrix>& betas = 
      *ddata.currBetaMat;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations
      
      betas[iState] = _kPlus[iState]/m_sumKplus;
    }

}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
