#include "STKT_NlimCSchemeS_nad.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "MathTools/MatrixInverter.hh"

#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STKT_NlimCSchemeS_nad,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeRiccNlimcSchemeSysNadProvider("STKT_SysNlimC_nad");

////////////////////////////////////////////////////////////////////////////

void STKT_NlimCSchemeS_nad::defineConfigOptions(Config::OptionList& options)
{
}
//////////////////////////////////////////////////////////////////////////////

STKT_NlimCSchemeS_nad::STKT_NlimCSchemeS_nad(const std::string& name) :
  STKT_SplitterSys(name),
  m_sumKplus(),
_invK(),
  m_uTemp(),
m_diss(),
  m_sumKplusU(),
      _phy(),
  _phyChar(),
  _residualChar(0),
  _rightEigenVector(),
  _leftEigenVector(),
  _sumBeta(),
  _beta(0),
  _betaLim(0),
  m_phiT()
{
  addConfigOptionsTo(this);


}

//////////////////////////////////////////////////////////////////////////////

STKT_NlimCSchemeS_nad::~STKT_NlimCSchemeS_nad()
{
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NlimCSchemeS_nad::setup()
{
  STKT_SplitterSys::setup();

  // Do not need a dissipation contribution from past 
 getMethodData().getDistributionData().needDiss = false;

  m_sumKplus.resize(_nbEquations, _nbEquations);
  m_uTemp.resize(_nbEquations);
  m_sumKplusU.resize(_nbEquations);
  m_diss.resize(_nbEquations);
  _invK.resize(_nbEquations, _nbEquations); 
  _phy.resize(_nbEquations);
  _phyChar.resize(_nbEquations);
  _rightEigenVector.resize(_nbEquations, _nbEquations);
  _leftEigenVector.resize(_nbEquations, _nbEquations);
  _sumBeta.resize(_nbEquations);
  m_phiT.resize(_nbEquations);

CFuint maxNbStatesInCell = _kPlus.size();

  _residualChar.resize(maxNbStatesInCell);
for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _residualChar[iState].resize(_nbEquations);
  }

  _beta.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _beta[iState].resize(_nbEquations);
  }

  _betaLim.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _betaLim[iState].resize(_nbEquations);
  }

    
   m_cterm = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<EulerTerm>();

   m_normal.resize(2);
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NlimCSchemeS_nad::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SpaceTimeLDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void STKT_NlimCSchemeS_nad::distribute(vector<RealVector>& residual)
{  
  
  const RealVector& m_phiT = getMethodData().getDistributionData().phi;
  const CFuint nbStatesInCell = _nbStatesInCell;
 const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;

SafePtr<EulerTerm> m_cterm = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<EulerTerm>();

  const RealVector& lData =  m_cterm->getPhysicalData();

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

   CFreal ux = lData[EulerTerm::VX];
   CFreal uy = lData[EulerTerm::VY];
   CFreal norm_u = sqrt(ux*ux+uy*uy);
   CFreal inv_norm_u = 1.0/norm_u;
   CFreal cos = 1.0/sqrt(2.0);

   if(norm_u > 1.e-2){
    m_normal[XX] = ux*inv_norm_u;
    m_normal[YY] = uy*inv_norm_u;
  }
  else {
    m_normal[XX] = cos;
    m_normal[YY] = cos;
  } 

 /// Compute the cell residual...
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _phy[iEq] = 0.0;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      _phy[iEq] += residual[iState][iEq];
    }
  }

  /// Transform into characteristic variables the residual, phy
  // Compute the EigenVectors
  getMethodData().getDistribVar()->computeEigenValuesVectors(_rightEigenVector,
							 _leftEigenVector,
							 _eValues,
							 m_normal);

  //Transform the residual into characteristic variables
  _phyChar = _leftEigenVector *_phy;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

      _residualChar[iState] = _leftEigenVector * residual[iState];

  }

  // Do the limiting just as in scalar case
  

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    if (std::abs(_phyChar[iEq])>MathTools::MathConsts::CFrealEps()){
      _sumBeta[iEq] = 0.;
      for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

        _beta[iState][iEq] = _residualChar[iState][iEq]/_phyChar[iEq];

        _sumBeta[iEq] += max(0.,_beta[iState][iEq]);

      }
    }
  }

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    if (std::abs(_phyChar[iEq])>MathTools::MathConsts::CFrealEps()){
      for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
        _betaLim[iState][iEq] = max(0.,_beta[iState][iEq])/_sumBeta[iEq];
        _residualChar[iState][iEq] = _betaLim[iState][iEq] * _phyChar[iEq];
        }
    }
    else{
      for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
        _residualChar[iState][iEq] = 0.;
      }
    }
  }

  //Transform the residual back into conservative variables
  _phy = _rightEigenVector * _phyChar;
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
     residual[iState] = _rightEigenVector * _residualChar[iState];
  }
  
  // AL: compute betas on the fly is more expensive because
  // it involves a matrix*matrix and a matrix*vector
  // while normal scheme requires only two matrix*vector operations
  vector<RealMatrix>& betas = *getMethodData().getDistributionData().currBetaMat;
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    betas[iState] = (*_kPlus[iState])*_invK;
  }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
