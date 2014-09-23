#include "SpaceTime2LayerNlimSchemeScalar.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "Framework/MethodStrategyProvider.hh"
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

MethodStrategyProvider<SpaceTime2LayerNlimSchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
SpaceTime2LayerNlimSchemeScalarProvider("ST2ScalarNlim");

//////////////////////////////////////////////////////////////////////////////

SpaceTime2LayerNlimSchemeScalar::SpaceTime2LayerNlimSchemeScalar(const std::string& name) :
  STM_SplitterScalar(name),
  _sumKmin(),
  _sumKminU(),
  _uInflow(),
  _uMin(),
  _uTemp(),
  _phy(),
  _sumBeta(),
  _sumBetaLim(),
  _beta(0),
  _betaLim(0)
{
}

//////////////////////////////////////////////////////////////////////////////

SpaceTime2LayerNlimSchemeScalar::~SpaceTime2LayerNlimSchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNlimSchemeScalar::setup()
{
  STM_SplitterScalar::setup();

  _sumKmin.resize(_nbEquations);
  _sumKminU.resize(_nbEquations);
  _uInflow.resize(_nbEquations);
  _uMin.resize(_nbEquations);
  _uTemp.resize(_nbEquations);
  _phy.resize(_nbEquations);
  _sumBeta.resize(_nbEquations);
  _sumBetaLim.resize(_nbEquations);

  CFuint maxNbStatesInCell = _kPlus.size();

  _beta.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _beta[iState].resize(_nbEquations);
  }

  _interBeta.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _interBeta[iState].resize(_nbEquations);
  }

  _resTemp.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _resTemp[iState].resize(_nbEquations);
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNlimSchemeScalar::distributePast(const vector<State*>& tStates)
{
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
  _sumKminU = _kMin[0]*(*tStates[0]);
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += _kMin[iState];
    _sumKminU += _kMin[iState] * (*tStates[iState]);
  }
  _uInflow = _sumKminU/_sumKmin;

  CFLogDebugMax( "============== SPATIAL PART: PAST STATES ============ " << "\n");
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]) - _uInflow;
    _uMin *= _kPlus[iState];
    _uMin *= (_timeStep/2.);
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

void SpaceTime2LayerNlimSchemeScalar::distributeInterK1(const vector<State*>& tStates,
                                     vector<RealVector>& residual)
{
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint nbEqs = _nbEquations;
  const CFreal Area = _cellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+1.);


  CFLogDebugMax( "============== TEMPORAL PART ============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = (Area*(*tStates[iState]));
    for (CFuint j = 0; j<nbEqs; ++j){
      residual[iState][j] += _uMin[j];
    }
  }

  CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]);
    _uMin /= 2.;
    _uMin *= (Area - pastArea);
    for (CFuint j = 0; j<nbEqs; ++j){
      residual[iState][j] -= _uMin[j];
      }
  }

  _sumKmin = _kMin[0];
  _sumKminU = _kMin[0]*(*tStates[0]);
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += _kMin[iState];
    _sumKminU += _kMin[iState] * (*tStates[iState]);
  }
  CFLogDebugMax( "_sumKmin = " << _sumKmin << "\n");
  _uInflow = _sumKminU/_sumKmin;

  CFLogDebugMax( "============== SPATIAL PART: CURRENT STATES ============ " << "\n");
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]) - _uInflow;
    _uMin *= _kPlus[iState];
    _uMin *= (_timeStep/2.);
    for (CFuint iEq=0; iEq<nbEqs; ++iEq){
      residual[iState][iEq] += _uMin[iEq];
    }
    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
  }

  /// Limiting
  // Compute the residual of cell K1
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
      _betaLim = 0.;
      if ((std::abs(_phy[iEq])>Eps)){
        _betaLim = max(0.,_beta[iState][iEq])/_sumBeta[iEq];
      }
      residual[iState][iEq] = _betaLim*_phy[iEq];
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNlimSchemeScalar::distributeInterK2(const vector<State*>& tStates,
                                     vector<RealVector>& residual)
{
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint nbEqs = _nbEquations;

  _sumKmin = _kMin[0];
  _sumKminU = _kMin[0]*(*tStates[0]);
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += _kMin[iState];
    _sumKminU += _kMin[iState] * (*tStates[iState]);
  }
  _uInflow = _sumKminU/_sumKmin;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]) - _uInflow;
    _uMin *= _kPlus[iState];
    _uMin *= (_timeStep/2.);
    for (CFuint j = 0; j<nbEqs; ++j){
      _resTemp[iState][j] = _uMin[j];
    }
    CFLogDebugMax( "residual[" << iState << "] = " << _resTemp[iState] << "\n");
  }

}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNlimSchemeScalar::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal Area = _cellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFuint nbEqs = _nbEquations;
  //const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  CFLogDebugMax( "============== TEMPORAL PART ============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = *tStates[iState];
    _uMin *= Area;
    _uTemp = *_interStates[iState];
    _uTemp *= pastArea;
    _uMin -= _uTemp;
    for (CFuint j = 0; j < nbEqs; ++j){
      residual[iState][nbEqs+j] = _uMin[j];
    }
  }

  CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = *tStates[iState] + *_interStates[iState];
    _uMin /= 2.;
    _uMin *= (Area - pastArea);
    for (CFuint j = 0; j < nbEqs; ++j){
      residual[iState][nbEqs+j] -= _uMin[j];
    }
  }

  _sumKmin = _kMin[0];
  _sumKminU = _kMin[0]*(*tStates[0]);
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += _kMin[iState];
    _sumKminU += _kMin[iState] * (*tStates[iState]);
  }
  CFLogDebugMax( "_sumKmin = " << _sumKmin << "\n");
  _uInflow = _sumKminU/_sumKmin;

  CFLogDebugMax( "============== SPATIAL PART: CURRENT STATES ============ " << "\n");
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]) - _uInflow;
    _uMin *= _kPlus[iState];
    _uMin *= (_timeStep/2.);
    for (CFuint j = 0; j < nbEqs; ++j){
      residual[iState][nbEqs+j] += _uMin[j];
    }
    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
  }

  // Compute the residual of cell K2
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _phy[iEq] = 0.0;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // Contribution from inter states
      _phy[iEq] += _resTemp[iState][iEq];
      // Contribution from current states
      _phy[iEq] += residual[iState][nbEqs+iEq];
    }
  }

  CFreal Eps = 0.;
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _sumBeta = 0.0;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      _beta[iState][iEq] = 0.;
      _interBeta[iState][iEq] = 0.;
      if (std::abs(_phy[iEq])>Eps){
        _beta[iState][iEq] = residual[iState][nbEqs+iEq]/_phy[iEq];
        _interBeta[iState][iEq] = _resTemp[iState][iEq]/_phy[iEq];
        _sumBeta[iEq] += max(0.,_beta[iState][iEq]);
        _sumBeta[iEq] += max(0.,_interBeta[iState][iEq]);
        }
    }
  }

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      if ((std::abs(_phy[iEq])>Eps)){
        _betaLim = max(0.,_beta[iState][iEq])/_sumBeta[iEq];
        _interBetaLim = max(0.,_interBeta[iState][iEq])/_sumBeta[iEq];
        residual[iState][nbEqs+iEq] = _betaLim * _phy[iEq];
        residual[iState][iEq] += _interBetaLim * _phy[iEq];
      }
      else{
        residual[iState][nbEqs+iEq] = 0.;
        residual[iState][iEq] = 0.;
      }
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
