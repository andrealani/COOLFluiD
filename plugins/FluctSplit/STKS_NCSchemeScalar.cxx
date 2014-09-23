#include "STKS_NCSchemeScalar.hh"
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

MethodStrategyProvider<STKS_NCSchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
stksNcSchemeScalarProvider("STKS_ScalarNC");

//////////////////////////////////////////////////////////////////////////////

STKS_NCSchemeScalar::STKS_NCSchemeScalar(const std::string& name) :
  STKS_SplitterScalar(name),
  m_sumKplus(),
  m_uTemp(),
  m_diss(),
  m_sumKplusU(),
  past_diss(),
  m_time_comp()
{
 
}

//////////////////////////////////////////////////////////////////////////////

STKS_NCSchemeScalar::~STKS_NCSchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKS_NCSchemeScalar::setup()
{
  STKS_SplitterScalar::setup();

  m_sumKplus.resize(_nbEquations);
  m_uTemp.resize(_nbEquations);
  m_sumKplusU.resize(_nbEquations);
  m_diss.resize(_nbEquations);

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  past_diss.resize(maxNbStatesInCell);
  m_time_comp.resize(maxNbStatesInCell);

  for (CFuint iState = 0; iState < maxNbStatesInCell; ++ iState){ 
    past_diss[iState].resize(_nbEquations);
    m_time_comp[iState].resize(_nbEquations);
  } 
  

  //Need a dissipation contribution from the past
  getMethodData().getDistributionData().needDiss = true;
}

//////////////////////////////////////////////////////////////////////////////

void STKS_NCSchemeScalar::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeCaraNSchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////
  void STKS_NCSchemeScalar::ComputePastDissipationAndTimeComp( const vector<State*>& tStates )
{
  const CFuint nbStatesInCell = _nbStatesInCell;
 

  m_sumKplus = _kPlus[0];
  m_sumKplusU = _kPlus[0]*(*tStates[0]);

  for (CFuint iState = 1; iState < nbStatesInCell; ++ iState){
    m_sumKplus += _kPlus[iState];
    m_sumKplusU += _kPlus[iState]*(*tStates[iState]);
  } 

  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){ 
  past_diss[iState]  = _kPlus[iState]*(*tStates[iState])- (_kPlus[iState]*m_sumKplusU/m_sumKplus);
  m_time_comp[iState]= (-1.)*(*tStates[iState]);
  }
}


//////////////////////////////////////////////////////////////////////////////

void STKS_NCSchemeScalar::distribute(vector<RealVector>& residual)
{  
   ///@todo No mesh deformation implemented here!!
  const CFuint nbStatesInCell = _nbStatesInCell;
  DistributionData& ddata = getMethodData().getDistributionData();
  const RealVector& m_phiT = ddata.phi;
  const vector<State*>& tStates = *ddata.tStates;
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT()/dim;
  const CFreal Area = _cellVolume/(dim+1.);

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
    residual[iState] += timeStep*(m_diss + past_diss[iState]);
    
    //Finaly component from the time
    residual[iState] += Area*((*tStates[iState]) +  m_time_comp[iState]);
  }

  vector<RealMatrix>& betas = 
      *getMethodData().getDistributionData().currBetaMat;
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
