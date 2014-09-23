#include "Common/NotImplementedException.hh"

#include "MathTools/MatrixInverter.hh"

#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/STKT_BDNSSchemeSys.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STKT_BDNSSchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeRiccBdnsSchemeSysProvider("STKT_SysBDNS");

//////////////////////////////////////////////////////////////////////////////

STKT_BDNSSchemeSys::STKT_BDNSSchemeSys(const std::string& name) :
  STKT_SplitterSys(name),
  m_sumKplus(),
  _invK(),
  m_uTemp(),
  m_diss(),
  m_sumKplusU(),
  m_phiN(),
  m_phitot(),
  m_absphitot(),
  m_phiT(),
   socket_dampingCoeff("dampingCoeff")
{

  // addConfigOptionsTo(this);

   // m_r0 = 0.0;
   // setParameter("r0",&m_r0);

   // m_beta = 0.0;
   // setParameter("beta",&m_beta);
}

//////////////////////////////////////////////////////////////////////////////

STKT_BDNSSchemeSys::~STKT_BDNSSchemeSys()
{
}
//////////////////////////////////////////////////////////////////////////////
// void STKT_BDNSSchemeSys::defineConfigOptions(Config::OptionList& options)
// {
//   options.addConfigOption< CFreal > ("r0","Maximum limit of X where LDA is applied");
//   options.addConfigOption< CFreal >("beta","Maximum limit of Y where LDA is applied");
// }
//////////////////////////////////////////////////////////////////////////////

 std::vector<Common::SafePtr<BaseDataSocketSink> >
 STKT_BDNSSchemeSys::needsSockets()
 {
   std::vector<Common::SafePtr<BaseDataSocketSink> > result = STKT_SplitterSys::needsSockets();
   result.push_back(&socket_dampingCoeff);
         return result;
       }
//////////////////////////////////////////////////////////////////////////////
void STKT_BDNSSchemeSys::setup()
{
  STKT_SplitterSys::setup();
   // Do not need a dissipation contribution from past 
 getMethodData().getDistributionData().needDiss = false;
 CFuint maxNbStatesInCell = _kPlus.size();
  m_phiN.resize(maxNbStatesInCell);
  m_diss.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++ iState){
    m_phiN[iState].resize(_nbEquations);
    m_diss[iState].resize(_nbEquations);
 }
 m_sumKplus.resize(_nbEquations,_nbEquations);
 m_uTemp.resize(_nbEquations);
   m_sumKplusU.resize(_nbEquations);
   m_phitot.resize(_nbEquations);
  m_absphitot.resize(_nbEquations);
 _invK.resize(_nbEquations,_nbEquations);
  m_phiT.resize(_nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void STKT_BDNSSchemeSys::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void STKT_BDNSSchemeSys::distribute(vector<RealVector>& residual)
{  
    const RealVector& m_phiT = getMethodData().getDistributionData().phi;
  const CFuint nbStatesInCell = _nbStatesInCell;
 const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
CFreal theta;

 DataHandle<CFreal> dampingCoeff = socket_dampingCoeff.getDataHandle();
   
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
    m_phiN[iState] = *_kPlus[iState]*m_uTemp;
  }

   m_uTemp = _invK*m_sumKplusU;
 
   for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
     // And then the dissipation
     m_diss[iState] = *_kPlus[iState]*(*tStates[iState]) - *_kPlus[iState]*m_uTemp;

   } 

   CFuint IDstate = (tStates[0])->getLocalID();	
   theta = dampingCoeff[IDstate];	
  // //computation of theta
  for (CFuint iState = 1; iState < nbStatesInCell; ++ iState){
       IDstate = (tStates[iState])->getLocalID();
       theta += dampingCoeff[IDstate];
       }
   theta /= nbStatesInCell;

   for(CFuint iEq = 0; iEq < _nbEquations; ++iEq){
    for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
       CFuint IDstate = (tStates[iState])->getLocalID();
  //      // residual[iState][iEq] += theta*m_diss[iState][iEq];
       residual[iState][iEq] += theta*m_diss[iState][iEq];
  
     }
   }

    // AL: compute betas on the fly is more expensive because
    // it involves a matrix*matrix and a matrix*vector
    // while normal scheme requires only two matrix*vector operations
    vector<RealMatrix>& betas = *getMethodData().getDistributionData().currBetaMat;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
    {
      betas[iState] = (*_kPlus[iState])*_invK;
    }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
