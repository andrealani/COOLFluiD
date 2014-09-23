#include "STKS_NlimCSchemeScalar.hh"
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

MethodStrategyProvider<STKS_NlimCSchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
stksNlimcSchemeScalarProvider("STKS_ScalarNlimC");

//////////////////////////////////////////////////////////////////////////////

STKS_NlimCSchemeScalar::STKS_NlimCSchemeScalar(const std::string& name) :
  STKS_SplitterScalar(name),
  m_sumKplus(),
  m_uTemp(),
  m_diss(),
  m_sumKplusU(),
  past_diss(),
  m_time_comp(),
  m_phiN(),
  m_betaN(),
  m_betaNPlus(),
  m_phitot(),
  m_sumBetaNPlus(),
  m_beta()
{

}

//////////////////////////////////////////////////////////////////////////////

STKS_NlimCSchemeScalar::~STKS_NlimCSchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKS_NlimCSchemeScalar::setup()
{
  STKS_SplitterScalar::setup();

  m_sumKplus.resize(_nbEquations);
 m_uTemp.resize(_nbEquations);
   m_sumKplusU.resize(_nbEquations);
 m_diss.resize(_nbEquations);

 const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

past_diss.resize(maxNbStatesInCell);

  for (CFuint iState = 0; iState < maxNbStatesInCell; ++ iState){

   past_diss[iState].resize(_nbEquations);
}
  m_time_comp.resize(maxNbStatesInCell);

  for (CFuint iState = 0; iState < maxNbStatesInCell; ++ iState){

   m_time_comp[iState].resize(_nbEquations);
}
getMethodData().getDistributionData().needDiss = true;

m_sumBetaNPlus.resize(_nbEquations);
m_phitot.resize(_nbEquations);
  m_phiN.resize(maxNbStatesInCell);
  m_betaN.resize(maxNbStatesInCell);
  m_betaNPlus.resize(maxNbStatesInCell);
  m_beta.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    m_phiN[iState].resize(_nbEquations);
    m_betaN[iState].resize(_nbEquations);
    m_betaNPlus[iState].resize(_nbEquations);
    m_beta[iState].resize(_nbEquations);
  }

}

//////////////////////////////////////////////////////////////////////////////

void STKS_NlimCSchemeScalar::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////
  void STKS_NlimCSchemeScalar::ComputePastDissipationAndTimeComp( const vector<State*>& tStates )
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

void STKS_NlimCSchemeScalar::distribute(vector<RealVector>& residual)
{
  /// @todo No mesh deformation implemented here!!
  const CFuint nbStatesInCell = _nbStatesInCell;
  DistributionData& ddata = getMethodData().getDistributionData();
  const RealVector& m_phiT = ddata.phi;
//   const RealVector& phiT_time = ddata.phi_time;
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

  m_phitot = 0.0;
  for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){
    // We first add the Residual of LDA
    m_uTemp = _kPlus[iState]/m_sumKplus;
   m_phiN[iState] = m_uTemp*m_phiT;
  // And then the dissipation
  m_diss = _kPlus[iState]*(*tStates[iState])- (_kPlus[iState]*m_sumKplusU/m_sumKplus);

  m_phiN[iState] += timeStep*(m_diss + past_diss[iState]);

  m_phiN[iState] += Area*((*tStates[iState]) +  m_time_comp[iState]);

  m_phitot += m_phiN[iState];
}

  //limitation of he N scheme
 m_sumBetaNPlus = 0.0;
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {

    m_betaN[iState] = (m_phiN[iState] + MathTools::MathConsts::CFrealEps()) / (m_phitot + nbStatesInCell*MathTools::MathConsts::CFrealEps());

    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      m_betaNPlus[iState][iEq] = std::max(m_betaN[iState][iEq],MathTools::MathConsts::CFrealEps());
      m_sumBetaNPlus[iEq] += m_betaNPlus[iState][iEq];
    }
  }

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    m_beta[iState] = m_betaNPlus[iState] / m_sumBetaNPlus;
    residual[iState] = m_beta[iState] *  m_phitot;
  }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
