#include "Common/NotImplementedException.hh"

#include "MathTools/MatrixInverter.hh"

#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FluctSplit/FluctSplitSystem.hh"
#include "FluctSplit/BDNSSchemeSys.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<BDNSSchemeSys,
		       FluctuationSplitData,
		       Splitter,
           FluctSplitSystemModule>

RiccBdnsSchemeSysProvider("SysBDNS");

//////////////////////////////////////////////////////////////////////////////

BDNSSchemeSys::BDNSSchemeSys(const std::string& name) :
  RDS_SplitterSys(name),
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
}

//////////////////////////////////////////////////////////////////////////////

BDNSSchemeSys::~BDNSSchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void BDNSSchemeSys::setup()
{
  RDS_SplitterSys::setup();
   // Do not need a dissipation contribution from past 
 getMethodData().getDistributionData().needDiss = false;
 CFuint maxNbStatesInCell = _kPlus.size();
 m_phiN.resize(maxNbStatesInCell);
 m_diss.resize(maxNbStatesInCell);
 for (CFuint iState = 0; iState < maxNbStatesInCell; ++ iState)
 {
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

std::vector<Common::SafePtr<BaseDataSocketSink> >
BDNSSchemeSys::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = RDS_SplitterSys::needsSockets();
  result.push_back(&socket_dampingCoeff);
        return result;
}

//////////////////////////////////////////////////////////////////////////////

void BDNSSchemeSys::distribute(vector<RealVector>& residual)
{  
  const RealVector& m_phiT = getMethodData().getDistributionData().phi;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  CFreal theta;

  DataHandle<CFreal> dampingCoeff = socket_dampingCoeff.getDataHandle();

  m_sumKplus  = *_kPlus[0];
  m_sumKplusU = *_kPlus[0]*(*tStates[0]);

  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState){
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
       //CFuint IDstate = (tStates[iState])->getLocalID();
       // residual[iState][iEq] += theta*m_diss[iState][iEq];
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

void BDNSSchemeSys::distributePart(vector<RealVector>& residual)
{

// THIS NEEDS TO BE IMPLEMENTED PROPERLY
//  _sumKplus = *_kPlus[0];
//  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
//    _sumKplus  += *_kPlus[iState];
//  }

//  _inverter->invert(_sumKplus, _invK);

//  RealVector& phi = getMethodData().getDistributionData().phi;
//  _uTemp.slice(0, _nbEquations) = _invK * phi.slice(_firstVarID, _nbEquations);

//  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
//    residual[iState].slice(_firstVarID, _nbEquations) =
//      (*_kPlus[iState]) * _uTemp.slice(0, _nbEquations);
//  }
}


//////////////////////////////////////////////////////////////////////////////

void BDNSSchemeSys::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
// THIS NEEDS TO BE IMPLEMENTED PROPERLY
//  throw Common::NotImplementedException (FromHere(),"SysBDNS::computePicardJacob()");
  // carefull with the signs !!!
//  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
//    // beta coefficient
//    _beta = (*_kPlus[iState])*_invK;
//    const CFuint nStart = iState*_nbStatesInCell;

//    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
//      RealMatrix *const block = jacob[nStart + jState];

//      _k =  (*_kPlus[jState]) + (*_kMin[jState]);
//      (*block) = _beta*_k;
//    }
//  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
