#include "SpaceTime2LayerNSchemeSys.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "MathTools/MatrixInverter.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<SpaceTime2LayerNSchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTime2LayerNSchemeSysProvider("ST2SysN");

//////////////////////////////////////////////////////////////////////////////

SpaceTime2LayerNSchemeSys::SpaceTime2LayerNSchemeSys(const std::string& name) :
  STM_SplitterSys(name),
  _sumKmin(),
  _invK(),
  _sumKminU(),
  _uInflow(),
  _tempMat(),
  _uMin(),
  _uTemp()
{
}

//////////////////////////////////////////////////////////////////////////////

SpaceTime2LayerNSchemeSys::~SpaceTime2LayerNSchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNSchemeSys::setup()
{
  STM_SplitterSys::setup();

  _sumKmin.resize(_nbEquations,_nbEquations);
  _invK.resize(_nbEquations,_nbEquations);
  _sumKminU.resize(_nbEquations);
  _uInflow.resize(_nbEquations);
  _tempMat.resize(_nbEquations, _nbEquations);
  _uMin.resize(_nbEquations);
  _uTemp.resize(_nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNSchemeSys::distributePast(const vector<State*>& tStates)
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
    for (CFuint j=0; j< nbEqs; ++j){
        past_residuals[(iState*nbEqs)+j] = _uMin[j];
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

    for (CFuint j=0; j< nbEqs; ++j){
        past_residuals[(iState*nbEqs)+j] -= _uMin[j];
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
    for (CFuint j=0; j< nbEqs; ++j){
        past_residuals[(iState*nbEqs)+j] -= _uMin[j];
    }
  }


}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNSchemeSys::distributeInterK1(const vector<State*>& tStates,
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

}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNSchemeSys::distributeInterK2(const vector<State*>& tStates,
                                     vector<RealVector>& residual)
{
  // unused //  const CFreal Area = _cellVolume;
  // unused //  const CFreal pastArea = _pastCellVolume;
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
      residual[iState][iEq] += _uMin[iEq];
    }

  }

}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime2LayerNSchemeSys::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const CFreal Area = _cellVolume;
  const CFreal pastArea = _pastCellVolume;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint nbEqs = _nbEquations;
  CFreal coef1 = 0.0;
  CFreal coef2 = 0.0;

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
    for (CFuint iEq=0; iEq<nbEqs; ++iEq){
      residual[iState][nbEqs+iEq] = _uMin[iEq];
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

    for (CFuint iEq=0; iEq < nbEqs; ++iEq){
      residual[iState][nbEqs+iEq] += _uMin[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
