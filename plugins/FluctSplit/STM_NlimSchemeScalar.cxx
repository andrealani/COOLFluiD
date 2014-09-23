#include "STM_NlimSchemeScalar.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STM_NlimSchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeNlimSchemeScalarProvider("STM_ScalarNlim");

//////////////////////////////////////////////////////////////////////////////

STM_NlimSchemeScalar::STM_NlimSchemeScalar(const std::string& name) :
  STM_SplitterScalar(name),
  _sumKmin(),
  _sumKminU(),
  _sumKplus(),
  _uMin(),
  _uTemp(),
  _uInflow(),
  _phy(),
  _sumBeta(),
  _sumBetaLim(),
  _beta(0),
  _betaLim(0)
{
}

//////////////////////////////////////////////////////////////////////////////

STM_NlimSchemeScalar::~STM_NlimSchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void STM_NlimSchemeScalar::setup()
{
  STM_SplitterScalar::setup();

  _sumKmin.resize(_nbEquations);
  _sumKminU.resize(_nbEquations);
  _sumKplus.resize(_nbEquations);
  _uMin.resize(_nbEquations);
  _uTemp.resize(_nbEquations);
  _uInflow.resize(_nbEquations);
  _phy.resize(_nbEquations);
  _sumBeta.resize(_nbEquations);
  _sumBetaLim.resize(_nbEquations);

  CFuint maxNbStatesInCell = _kPlus.size();

  _beta.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _beta[iState].resize(_nbEquations);
  }

  _betaLim.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _betaLim[iState].resize(_nbEquations);
  }

}

//////////////////////////////////////////////////////////////////////////////

void STM_NlimSchemeScalar::distributePast(const vector<State*>& tStates)
{
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint nbEqs = _nbEquations;
  const CFreal Area = _cellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+1.);

  RealVector& past_residuals = *getMethodData().getDistributionData().past_residuals;


  CFLogDebugMax( "============== TEMPORAL PART ============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]);
    _uMin *= -pastArea;
    for (CFuint j=0; j< nbEqs; ++j){
        past_residuals[(iState*nbEqs)+j] = _uMin[j];
    }
  }

  _sumKmin = _kMin[0];
  _sumKplus = _kPlus[0];
  _sumKminU = _kMin[0]*(*tStates[0]);
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += _kMin[iState];
    _sumKplus  += _kPlus[iState];
    _sumKminU += _kMin[iState] * (*tStates[iState]);
  }
  _uInflow = _sumKminU/_sumKmin;

  CFLogDebugMax( "============== SPATIAL PART: PAST STATES ============ " << "\n");
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]) - _uInflow;
    _uMin *= _kPlus[iState];
    _uMin *= (timeStep/2.);
    past_residuals.slice(iState*nbEqs, nbEqs) += _uMin.slice(0, nbEqs);
    CFLogDebugMax( "residual[" << iState << "] = " << past_residuals[iState] << "\n");
  }

  CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]);
    _uMin /= 2.;
    _uMin *= (Area - pastArea);
    past_residuals.slice(iState*nbEqs, nbEqs) -= _uMin.slice(0, nbEqs);
  }

 
  
}

//////////////////////////////////////////////////////////////////////////////

void STM_NlimSchemeScalar::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal Area = _cellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+1.);

  const CFuint nbEqs = _nbEquations;

  const RealVector& past_residuals = *getMethodData().getDistributionData().past_residuals;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
      for (CFuint jEq = 0; jEq <  nbEqs; ++jEq)
        residual[iState][jEq] = past_residuals[iState*nbEqs + jEq];

  CFLogDebugMax( "============== TEMPORAL PART ============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = *tStates[iState];
    _uMin *= Area;
    residual[iState] += _uMin;
  }

  CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]);
    _uMin /= 2.;
    _uMin *= (Area - pastArea);
    residual[iState] -= _uMin;
  }

 _sumKmin = _kMin[0];
  _sumKplus = _kPlus[0];
  _sumKminU = _kMin[0]*(*tStates[0]);
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += _kMin[iState];
      _sumKplus  += _kPlus[iState];
    _sumKminU += _kMin[iState] * (*tStates[iState]);
  }
  CFLogDebugMax( "_sumKmin = " << _sumKmin << "\n");
  _uInflow = _sumKminU/_sumKmin;

  CFLogDebugMax( "============== SPATIAL PART: CURRENT STATES ============ " << "\n");
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]) - _uInflow;
    _uMin *= _kPlus[iState];
    _uMin *= (timeStep/2.);
    residual[iState] += _uMin;
    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
  }


  // Compute the cell residual...
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _phy[iEq] = 0.0;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      _phy[iEq] += residual[iState][iEq];
    }
  }


  CFreal Eps = 0.;
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _sumBeta = 0.0;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      _beta[iState][iEq] = 0.;
      if (std::abs(_phy[iEq])>Eps){
        _beta[iState][iEq] = residual[iState][iEq]/_phy[iEq];
        _sumBeta[iEq] += max(0.,_beta[iState][iEq]);
        }
    }
  }

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      _betaLim[iState][iEq] = 0.;
      if ((std::abs(_phy[iEq])>Eps)){
        _betaLim[iState][iEq] = max(0.,_beta[iState][iEq])/_sumBeta[iEq];
      }
      residual[iState][iEq] = _betaLim[iState][iEq]*_phy[iEq];
    }
  }
 vector<RealMatrix>& betas = 
      *getMethodData().getDistributionData().currBetaMat;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations
      
      betas[iState] = _kPlus[iState]/_sumKplus;
    }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
