#include "SpaceTime2LayerNlimSchemeSys.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "MathTools/MatrixInverter.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Common/BadValueException.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<SpaceTime2LayerNlimSchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
SpaceTime2LayerNlimSchemeSysProvider("ST2SysNlim");

//////////////////////////////////////////////////////////////////////////////

SpaceTime2LayerNlimSchemeSys::SpaceTime2LayerNlimSchemeSys(const std::string& name) :
  STM_SplitterSys(name),
  _sumKmin(),
  _invK(),
  _sumKminU(),
  _uInflow(),
  _tempMat(),
  _uMin(),
  _uTemp(),
  _phy(),
  _phyChar(),
  _residualChar(0),
  _rightEigenVector(),
  _leftEigenVector(),
  _sumBeta(),
  _sumBetaLim(),
  _resTemp(0),
  _beta(0),
  _interBeta(0),
  _betaLim(0),
  _interBetaLim(0)
{
}

//////////////////////////////////////////////////////////////////////////////

SpaceTime2LayerNlimSchemeSys::~SpaceTime2LayerNlimSchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNlimSchemeSys::setup()
{
  STM_SplitterSys::setup();

  _sumKmin.resize(_nbEquations,_nbEquations);
  _invK.resize(_nbEquations,_nbEquations);
  _sumKminU.resize(_nbEquations);
  _uInflow.resize(_nbEquations);
  _tempMat.resize(_nbEquations, _nbEquations);
  _uMin.resize(_nbEquations);
  _uTemp.resize(_nbEquations);
  _phy.resize(_nbEquations);
  _phyChar.resize(_nbEquations);
  _rightEigenVector.resize(_nbEquations, _nbEquations);
  _leftEigenVector.resize(_nbEquations, _nbEquations);
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

  _residualChar.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _residualChar[iState].resize(_nbEquations);
  }

  _interResidualChar.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _interResidualChar[iState].resize(_nbEquations);
  }

}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNlimSchemeSys::distributePast(const vector<State*>& tStates)
{
  const CFreal Area = _cellVolume;
  const CFreal pastArea = _pastCellVolume;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint nbEqs = _nbEquations;
  CFreal coef1;
  CFreal coef2;
  
  RealVector& past_residuals = *getMethodData().getDistributionData().past_residuals;

  switch (PhysicalModelStack::getActive()->getDim()) {
  case DIM_1D:
    coef1 = 0.25;
    coef2 = 0.;
    break;
  case DIM_2D:
    coef1 = 1./12.;
    coef2 = 0.125;
    break;
  case DIM_3D:
    coef1 = 2./36.;
    coef2 = 1./72.;
    break;
  default:
    std::string msg = std::string("Bad dimension. Cannot be larger than 3D !!");
    throw BadValueException (FromHere(),msg);
  }


  CFLogDebugMax( "============== SPATIAL PART - PAST ============ " << "\n");

  _sumKmin  = *_kMin[0];
  _sumKminU = *_kMin[0]* *tStates[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += *_kMin[iState];
    _sumKminU += *_kMin[iState]* *tStates[iState];
  }

  _inverter->invert(_sumKmin, _invK);
  _uTemp = _invK * _sumKminU;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin  = *tStates[iState];
    _uMin -= _uTemp;
    _uInflow = *_kPlus[iState]*_uMin;
    _uMin  = _uInflow * (_timeStep/2.);
    for (CFuint iEq=0; iEq < nbEqs; ++iEq){
        past_residuals[(iState*nbEqs)+iEq] = _uMin[iEq];
    }
  }

  CFLogDebugMax( "============== TIME PART - PAST ============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

  // ********************************************
  // Conservative Variables !!!!
  // ********************************************
    _uMin = *_consStates[iState];
    _uMin *= coef1;

  // Transformed Linearized variables
    _uTemp = *tStates[iState];
    _uTemp *= coef2;

    _uMin += _uTemp;
    _uMin *= pastArea;

    for (CFuint iEq=0; iEq < nbEqs; ++iEq){
        past_residuals[(iState*nbEqs)+iEq] -= _uMin[iEq];
    }
  }

  CFLogDebugMax( "============== DEFORMING MESH PART - PAST ============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = *_consStates[iState];
    _uMin *= coef1;

  // Transformed Linearized variables
    _uTemp = *tStates[iState];
    _uTemp *= coef2;

    _uMin += _uTemp;

    _uMin /= 2.;
    _uMin *= (Area-pastArea);
    for (CFuint iEq=0; iEq < nbEqs; ++iEq){
        past_residuals[(iState*nbEqs)+iEq] -= _uMin[iEq];
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNlimSchemeSys::distributeInterK1(const vector<State*>& tStates,
                                     vector<RealVector>& residual)
{
  const CFreal Area = _cellVolume;
  const CFreal pastArea = _pastCellVolume;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint nbEqs = _nbEquations;
  CFreal coef1;
  CFreal coef2;

  switch (PhysicalModelStack::getActive()->getDim()) {
  case DIM_1D:
    coef1 = 0.25;
    coef2 = 0.;
    break;
  case DIM_2D:
    coef1 = 1./12.;
    coef2 = 0.125;
    break;
  case DIM_3D:
    coef1 = 2./36.;
    coef2 = 1./72.;
    break;
  default:
    std::string msg = std::string("Bad dimension. Cannot be larger than 3D !!");
    throw BadValueException (FromHere(),msg);
  }

  CFLogDebugMax( "============== SPATIAL PART: CURRENT STATES ============ " << "\n");

  // for efficiency reason iState = 0 is treated out of the loop
  // this allows to avoid resetting _sumKminU and _sumKmin to 0.0
  // after every loop (operation that would be done
  // nbTimeSteps*nbStatesPerCell*nbCells times!!)
  _sumKmin  = *_kMin[0];
  _sumKminU = *_kMin[0]* *tStates[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += *_kMin[iState];
    _sumKminU += *_kMin[iState]* *tStates[iState];
  }

  _inverter->invert(_sumKmin, _invK);
  _uTemp = _invK * _sumKminU;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin  = *tStates[iState];
    _uMin -= _uTemp;
    _uInflow = *_kPlus[iState]*_uMin;
    _uMin  = _uInflow * (_timeStep/2.);
    for (CFuint iEq=0; iEq < nbEqs; ++iEq){
      residual[iState][iEq] += _uMin[iEq];
    }

  }

  CFLogDebugMax( "============== TIME PART - CURRENT ============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
  // ********************************************
  // Conservative Variables !!!!
  // ********************************************
    _uMin = *_consStates[iState];
    _uMin *= coef1;

  // Transformed Linearized variables
    _uTemp = *tStates[iState];
    _uTemp *= coef2;

    _uMin += _uTemp;
    _uMin *= Area;

    for (CFuint iEq=0; iEq < nbEqs; ++iEq){
      residual[iState][iEq] += _uMin[iEq];
    }
  }

  CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
  // Conservative Variables
    _uMin = *_consStates[iState];
    _uMin *= coef1;

  // Transformed Linearized variables
    _uTemp = *tStates[iState];
    _uTemp *= coef2;

    _uMin += _uTemp;

    _uMin /= 2.;
    _uMin *= (Area-pastArea);
    for (CFuint iEq=0; iEq < nbEqs; ++iEq){
      residual[iState][iEq] -= _uMin[iEq];
      }
  }

  /// Limiting
    
  /// Compute the cell residual...(cell K1)
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
							 _adimNormal);

  //Transform the residual into characteristic variables
  _phyChar = _leftEigenVector * _phy;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _residualChar[iState].slice(0, _nbEquations) = _leftEigenVector * residual[iState].slice(0, _nbEquations);
  }

  // Do the limiting just as in scalar case
  CFreal Eps = 0.;

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _sumBeta[iEq] = 0.;
    if (std::abs(_phyChar[iEq])>Eps){
      for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
        _beta[iState][iEq] = _residualChar[iState][iEq]/_phyChar[iEq];
        _sumBeta[iEq] += max(0.,_beta[iState][iEq]);
      }
    }
  }

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    if (std::abs(_phyChar[iEq])>Eps){
      for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
        _betaLim = max(0.,_beta[iState][iEq])/_sumBeta[iEq];
        _residualChar[iState][iEq] = _betaLim * _phyChar[iEq];
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
    residual[iState].slice(0, _nbEquations) += _rightEigenVector * _residualChar[iState].slice(0, _nbEquations);
  }

}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNlimSchemeSys::distributeInterK2(const vector<State*>& tStates,
                                     vector<RealVector>& residual)
{
  // unused // const CFreal Area = _cellVolume;
  // unused // const CFreal pastArea = _pastCellVolume;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint nbEqs = _nbEquations;
  CFreal coef1;
  CFreal coef2;

  switch (PhysicalModelStack::getActive()->getDim()) {
  case DIM_1D:
    coef1 = 0.25;
    coef2 = 0.;
    break;
  case DIM_2D:
    coef1 = 1./12.;
    coef2 = 0.125;
    break;
  case DIM_3D:
    coef1 = 2./36.;
    coef2 = 1./72.;
    break;
  default:
    std::string msg = std::string("Bad dimension. Cannot be larger than 3D !!");
    throw BadValueException (FromHere(),msg);
  }

  // for efficiency reason iState = 0 is treated out of the loop
  // this allows to avoid resetting _sumKminU and _sumKmin to 0.0
  // after every loop (operation that would be done
  // nbTimeSteps*nbStatesPerCell*nbCells times!!)
  _sumKmin  = *_kMin[0];
  _sumKminU = *_kMin[0]* *tStates[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += *_kMin[iState];
    _sumKminU += *_kMin[iState]* *tStates[iState];
  }

  _inverter->invert(_sumKmin, _invK);
  _uTemp = _invK * _sumKminU;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin  = *tStates[iState];
    _uMin -= _uTemp;
    _uInflow = *_kPlus[iState]*_uMin;
    _uMin  = _uInflow * (_timeStep/2.);
    for (CFuint iEq=0; iEq < nbEqs; ++iEq){
      _resTemp[iState][iEq] = _uMin[iEq];
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNlimSchemeSys::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const CFreal Area = _cellVolume;
  const CFreal pastArea = _pastCellVolume;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint nbEqs = _nbEquations;
  CFreal coef1;
  CFreal coef2;

  switch (PhysicalModelStack::getActive()->getDim()) {
  case DIM_1D:
    coef1 = 0.25;
    coef2 = 0.;
    break;
  case DIM_2D:
    coef1 = 1./12.;
    coef2 = 0.125;
    break;
  case DIM_3D:
    coef1 = 2./36.;
    coef2 = 1./72.;
    break;
  default:
    std::string msg = std::string("Bad dimension. Cannot be larger than 3D !!");
    throw BadValueException (FromHere(),msg);
  }

  CFLogDebugMax( "============== SPATIAL PART: CURRENT STATES ============ " << "\n");

  // for efficiency reason iState = 0 is treated out of the loop
  // this allows to avoid resetting _sumKminU and _sumKmin to 0.0
  // after every loop (operation that would be done
  // nbTimeSteps*nbStatesPerCell*nbCells times!!)
  _sumKmin  = *_kMin[0];
  _sumKminU = *_kMin[0]* *tStates[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += *_kMin[iState];
    _sumKminU += *_kMin[iState]* *tStates[iState];
  }

  _inverter->invert(_sumKmin, _invK);
  _uTemp = _invK * _sumKminU;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    _uMin  = *tStates[iState];
    _uMin -= _uTemp;
    _uInflow = *_kPlus[iState]*_uMin;
    _uMin  = _uInflow * (_timeStep/2.);
    for (CFuint j=0; j<nbEqs; ++j){
      residual[iState][nbEqs+j] = _uMin[j];
    }
  }

  CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
  // Conservative Variables
    _uMin = *_consStates[iState] - *_interConsStates[iState];
    _uMin *= coef1;

  // Transformed Linearized variables
    _uTemp = *tStates[iState] - *_interStates[iState];
    _uTemp *= coef2;

    _uMin += _uTemp;

    _uMin /= 2.;
    _uMin *= (Area-pastArea);
    for (CFuint iEq=0; iEq<nbEqs; ++iEq){
      residual[iState][nbEqs+iEq] -= _uMin[iEq];
      }
  }

  CFLogDebugMax( "============== TIME PART - CURRENT ============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
  // ********************************************
  // Conservative Variables !!!!
  // ********************************************
    _uMin = *_consStates[iState] - *_interConsStates[iState];
    _uMin *= coef1;

  // Transformed Linearized variables
    _uTemp = *tStates[iState] - *_interStates[iState];
    _uTemp *= coef2;

    _uMin += _uTemp;
    _uMin *= Area;

    for (CFuint j=0; j<nbEqs; ++j){
      residual[iState][nbEqs+j] += _uMin[j];
    }
  }

  /// Limiting for K2
  /// Compute the cell residual...
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _phy[iEq] = 0.0;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // Contribution from inter states
      _phy[iEq] += _resTemp[iState][iEq];
      // Contribution from current states
      _phy[iEq] += residual[iState][nbEqs+iEq];
    }
  }

  /// Transform into characteristic variables the residuals
  // Compute the EigenVectors
  getMethodData().getDistribVar()->computeEigenValuesVectors(_rightEigenVector,
							 _leftEigenVector,
							 _eValues,
							 _adimNormal);

  //Transform the residual into characteristic variables
  _phyChar = _leftEigenVector * _phy;
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _interResidualChar[iState] = _leftEigenVector * _resTemp[iState];
    _residualChar[iState].slice(0, _nbEquations) = _leftEigenVector * residual[iState].slice(_nbEquations, _nbEquations);
  }

  // Do the limiting just as in scalar case
  CFreal Eps = MathTools::MathConsts::CFrealEps();
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    if (std::abs(_phyChar[iEq])>Eps){
    _sumBeta[iEq] = 0.;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
        _interBeta[iState][iEq] = _interResidualChar[iState][iEq]/_phyChar[iEq];
        _beta[iState][iEq] = _residualChar[iState][iEq]/_phyChar[iEq];
        _sumBeta[iEq] += max(0.,_interBeta[iState][iEq]);
        _sumBeta[iEq] += max(0.,_beta[iState][iEq]);
      }
    }
  }

    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      if (std::abs(_phyChar[iEq])>Eps){
      for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
          _betaLim = max(0.,_beta[iState][iEq])/_sumBeta[iEq];
          _interBetaLim = max(0.,_interBeta[iState][iEq])/_sumBeta[iEq];
          _residualChar[iState][iEq] = _betaLim * _phyChar[iEq];
          _interResidualChar[iState][iEq] = _interBetaLim * _phyChar[iEq];
        }
      }
      else{
        for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
          _residualChar[iState][iEq] = 0.;
          _interResidualChar[iState][iEq] = 0.;
        }
      }
  }

  //Transform the residual back into conservative variables
  _phy = _rightEigenVector * _phyChar;
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _resTemp[iState] = _rightEigenVector * _interResidualChar[iState];
  }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      residual[iState][iEq] += _resTemp[iState][iEq];
    }
  }

 for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
   _resTemp[iState] = _rightEigenVector * _residualChar[iState];
 }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      residual[iState][nbEqs+iEq] = _resTemp[iState][iEq];
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
