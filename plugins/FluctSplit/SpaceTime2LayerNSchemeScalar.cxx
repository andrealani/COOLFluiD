#include "SpaceTime2LayerNSchemeScalar.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<SpaceTime2LayerNSchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTime2LayerNSchemeScalarProvider("ST2ScalarN");

//////////////////////////////////////////////////////////////////////////////

SpaceTime2LayerNSchemeScalar::SpaceTime2LayerNSchemeScalar(const std::string& name) :
  STM_SplitterScalar(name),
  _sumKmin(),
  _sumKminU(),
  _uInflow(),
  _uMin()
{
}

//////////////////////////////////////////////////////////////////////////////

SpaceTime2LayerNSchemeScalar::~SpaceTime2LayerNSchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNSchemeScalar::setup()
{
  STM_SplitterScalar::setup();

  _sumKmin.resize(_nbEquations);
  _sumKminU.resize(_nbEquations);
  _uInflow.resize(_nbEquations);
  _uMin.resize(_nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNSchemeScalar::distributePast(const vector<State*>& tStates)
{
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint nbEqs = _nbEquations;
  const CFreal Area = _cellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+1.);

  RealVector& past_residuals = *getMethodData().getDistributionData().past_residuals;

//   CFLogDebugMax( "============== TEMPORAL PART ============ " << "\n");
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

//   CFLogDebugMax( "============== SPATIAL PART: PAST STATES ============ " << "\n");
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]) - _uInflow;
    _uMin *= _kPlus[iState];
    _uMin *= (_timeStep/2.);
    past_residuals.slice(iState*nbEqs, nbEqs) += _uMin.slice(0, nbEqs);
//     CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
  }

//   CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]);
    _uMin /= 2.;
    _uMin *= (Area - pastArea);
    past_residuals.slice(iState*nbEqs, nbEqs) -= _uMin.slice(0, nbEqs);
  }

  
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNSchemeScalar::distributeInterK1(const vector<State*>& tStates,
						     vector<RealVector>& residual)
{
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint nbEqs = _nbEquations;

  const CFreal Area = _cellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+1.);

//  CFLogDebugMax( "============== TEMPORAL PART ============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = (Area*(*tStates[iState]));
    for (CFuint j = 0; j<nbEqs; ++j){
      residual[iState][j] += _uMin[j];
    }
  }

//  CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]);
    _uMin /= 2.;
    _uMin *= (Area - pastArea);
    for (CFuint j=0; j<nbEqs; ++j){
      residual[iState][j] -= _uMin[j];
      }
  }

  _sumKmin = _kMin[0];
  _sumKminU = _kMin[0]*(*tStates[0]);
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += _kMin[iState];
    _sumKminU += _kMin[iState] * (*tStates[iState]);
  }
//  CFLogDebugMax( "_sumKmin = " << _sumKmin << "\n");
  _uInflow = _sumKminU/_sumKmin;

//  CFLogDebugMax( "============== SPATIAL PART: CURRENT STATES ============ " << "\n");
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]) - _uInflow;
    _uMin *= _kPlus[iState];
    _uMin *= (_timeStep/2.);
    for (CFuint j=0; j<nbEqs; ++j){
      residual[iState][j] += _uMin[j];
    }
//    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNSchemeScalar::distributeInterK2(const vector<State*>& tStates,
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
      residual[iState][j] += _uMin[j];
    }
//     CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
  }

}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNSchemeScalar::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal Area = _cellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFuint nbEqs = _nbEquations;
  //const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

//   CFLogDebugMax( "============== TEMPORAL PART ============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = *tStates[iState] - *_interStates[iState];
    _uMin *= Area;
    for (CFuint j=0; j < nbEqs; ++j){
      residual[iState][nbEqs+j] = _uMin[j];
    }
  }

//   CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = *tStates[iState] + *_interStates[iState];
    _uMin /= 2.;
    _uMin *= (Area - pastArea);
    for (CFuint j=0; j < nbEqs; ++j){
      residual[iState][nbEqs+j] -= _uMin[j];
    }
  }

  _sumKmin = _kMin[0];
  _sumKminU = _kMin[0]*(*tStates[0]);
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += _kMin[iState];
    _sumKminU += _kMin[iState] * (*tStates[iState]);
  }
//   CFLogDebugMax( "_sumKmin = " << _sumKmin << "\n");
  _uInflow = _sumKminU/_sumKmin;

//   CFLogDebugMax( "============== SPATIAL PART: CURRENT STATES ============ " << "\n");
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]) - _uInflow;
    _uMin *= _kPlus[iState];
    _uMin *= (_timeStep/2.);
    for (CFuint j=0; j < nbEqs; ++j){
      residual[iState][nbEqs+j] += _uMin[j];
    }
//     CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
  }

}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
