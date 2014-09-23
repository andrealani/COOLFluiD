#include "STM_BSchemeSys.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "MathTools/MatrixInverter.hh"
#include "Common/BadValueException.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STM_BSchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeBSchemeSysProvider("STM_SysB");

//////////////////////////////////////////////////////////////////////////////

STM_BSchemeSys::STM_BSchemeSys(const std::string& name) :
  STM_SplitterSys(name),
  _sumKmin(),
  _sumKplus(),
  _invK(),
  _uMin(),
  _tempMat(),
  _phyt(),
  _phys(),
  _identity(),
  _sumKminU(),
  _uInflow(),
  _uTemp()
{
}

//////////////////////////////////////////////////////////////////////////////

STM_BSchemeSys::~STM_BSchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void STM_BSchemeSys::setup()
{
  STM_SplitterSys::setup();

  _sumKmin.resize(_nbEquations,_nbEquations);
  _invK.resize(_nbEquations,_nbEquations);
  _uMin.resize(_nbEquations);
  _tempMat.resize(_nbEquations,_nbEquations);
  _phys.resize(_nbEquations);
  _identity.resize(_nbEquations);
  _identity = 1.0;
   _sumKminU.resize(_nbEquations);
   _uInflow.resize(_nbEquations);
   _uTemp.resize(_nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void STM_BSchemeSys::distributePast(const vector<State*>& tStates)
{
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint nbEqs = _nbEquations;
  CFreal timeStep = SubSystemStatusStack::getActive()->getDT();
  const CFreal Area = _cellVolume;
  const CFreal pastArea = _pastCellVolume;
  CFreal coef1;
  CFreal coef2;
  residual_lda.resize(nbStatesInCell*nbEqs);
  residual_n.resize(nbStatesInCell*nbEqs);

  //Compute the residual of N scheme

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
    _uInflow *= (timeStep/2.);
    for (CFuint j=0; j< nbEqs; ++j){
        residual_n[(iState*nbEqs)+j] = _uInflow[j];
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
        residual_n[(iState*nbEqs)+j] -= _uMin[j];
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
        residual_n[(iState*nbEqs)+j] -= _uMin[j];
    }
  }

  _phyt.resize(_nbStatesInCell);
  for(CFuint i = 0; i < _nbStatesInCell; ++i){
    _phyt[i].resize(_nbEquations);
 }

  //Compute the residual of LDA scheme

  CFreal coeff_pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+2.);

  _phys = 0.0;
  coef1 = 1./3.;
  coef2 = 1./2.;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    // ********************************************
    // Conservative Variables !!!!
    // ********************************************

    _phyt[iState] = *_consStates[iState];
    _phyt[iState] *= -coef1;
    // Transformed Linearized variables
    _uMin = (*tStates[iState]);
    _uMin *= -coef2;
    _phyt[iState] += _uMin;
    _phyt[iState] *= coeff_pastArea;
    // Compute K matrix
    _tempMat = *_kPlus[iState];
    _tempMat += *_kMin[iState];
    // Compute Phy_s
    _phys += _tempMat * *tStates[iState];
    }

  _phys *= timeStep/2.;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _tempMat = *_kPlus[iState] * _invK;
    _tempMat *= -1.;
    _identity = 1./6.;
    _tempMat += _identity;
    _uMin = _tempMat * _phyt[iState];
    for (CFuint j=0; j< nbEqs; ++j){
      residual_lda[(iState*nbEqs)+j] = _uMin[j];
    }

    _identity = 1./12.;
    _tempMat = *_kPlus[iState] * _invK;
    _tempMat += _identity;

    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      if (jState != iState){
        for (CFuint j=0; j< nbEqs; ++j){

          _uMin = _tempMat * _phyt[jState];

          residual_lda[(iState*nbEqs)+j] -= _uMin[j];
        }
      }
    }
  }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _tempMat = *_kPlus[iState] * _invK;
    for (CFuint j=0; j< nbEqs; ++j){

      _uMin = _tempMat * _phys;

      residual_lda[(iState*nbEqs)+j] -= _uMin[j];
  }
   CFLogDebugMax( "residual[" << iState << "] = " << residual_lda[iState] << "\n");
 }


  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
    for (CFuint jEq = 0; jEq < nbEqs ; ++jEq){
      (*getMethodData().getDistributionData().past_residuals)[iState*nbEqs + jEq] = residual_lda[iState*nbEqs + jEq];
      (*getMethodData().getDistributionData().past_residuals_order1)[iState*nbEqs + jEq] = residual_n[iState*nbEqs + jEq];
    }

}

//////////////////////////////////////////////////////////////////////////////

void STM_BSchemeSys::distribute(vector<RealVector>& residual)
{
  // for efficiency reason iState = 0 is treated out of the loop
  // this allows to avoid resetting _sumKminU and _sumKmin to 0.0
  // after every loop (operation that would be done
  // nbTimeSteps*nbStatesPerCell*nbCells times!!)
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;

  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();
  const CFreal Area = _cellVolume;
  const CFreal pastArea = _pastCellVolume;
  const CFuint nbEqs = _nbEquations;

  residual_lda.resize(nbStatesInCell*nbEqs);

  residual_n.resize(nbStatesInCell*nbEqs);
  CFreal res;
  CFreal res_abs;
  CFreal theta;


   _phyt.resize(_nbStatesInCell);
   for(CFuint i = 0; i < _nbStatesInCell; ++i){
     _phyt[i].resize(_nbEquations);
   }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
       for (CFuint jEq = 0; jEq <  nbEqs; ++jEq){
         residual_lda[iState*nbEqs + jEq] = (*getMethodData().getDistributionData().past_residuals)[iState*nbEqs + jEq];
         residual_n[iState*nbEqs + jEq] = (*getMethodData().getDistributionData().past_residuals_order1)[iState*nbEqs + jEq];
  }

   //Compute the residual of N scheme
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
    _uInflow *= (timeStep/2.);
  for (CFuint jEq = 0; jEq <  nbEqs; ++jEq)
    residual_n[iState*nbEqs+jEq] += _uInflow[jEq];
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
    for (CFuint jEq = 0; jEq <  nbEqs; ++jEq)
      residual_n[iState*nbEqs+jEq] += _uMin[jEq];
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
    for (CFuint jEq = 0; jEq <  nbEqs; ++jEq)
      residual_n[iState*nbEqs+jEq] -= _uMin[jEq];
  }

  //Compute the residual of LDA scheme
  CFreal coeff_Area = _cellVolume*(1./(PhysicalModelStack::getActive()->getDim()+2.)) ;

   _phys = 0.0;
   coef1 = 1./3.;
   coef2 = 1./2.;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
   // ********************************************
   // Conservative Variables !!!!
   // ********************************************
     _phyt[iState] = *_consStates[iState];
     _phyt[iState] *= coef1;

   // Transformed Linearized variables
     _uMin = (*tStates[iState]);
     _uMin *= coef2;

     _phyt[iState] += _uMin;
     _phyt[iState] *= coeff_Area;

    // Compute K matrix
    _tempMat = *_kPlus[iState];
    _tempMat += *_kMin[iState];
    // Compute Phy_s
    _phys += _tempMat * *tStates[iState];
  }
  _phys *= timeStep/2.;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _tempMat = *_kPlus[iState] * _invK;
    _tempMat *= -1.;
    _identity = 1./6.;
    _tempMat += _identity;
    _uTemp = _tempMat * _phyt[iState];
  for (CFuint jEq = 0; jEq <  nbEqs; ++jEq)
    residual_lda[iState*nbEqs+jEq] += _uTemp[jEq];

    _identity = 1./12.;
    _tempMat = *_kPlus[iState] * _invK;
    _tempMat += _identity;
    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      if (jState != iState){
        _uTemp = _tempMat * _phyt[jState];
  for (CFuint jEq = 0; jEq <  nbEqs; ++jEq)
    residual_lda[iState*nbEqs+jEq] -= _uTemp[jEq];

      }
    }
  }


  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _tempMat = *_kPlus[iState] * _invK;
    _uTemp = _tempMat * _phys;
  for (CFuint jEq = 0; jEq <  nbEqs; ++jEq)
    residual_lda[iState*nbEqs+jEq] -= _uTemp[jEq];

    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
  }
theta = -MathTools::MathConsts::CFrealMax();
  //compute the blending coefficient theta and blend the two schemes
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
   {
    res = residual_n[iEq];
    res_abs = std::fabs(residual_n[iEq]);

    for (CFuint iState = 1; iState < _nbStatesInCell; ++iState){
      res += residual_n[iState*nbEqs+iEq];
      res_abs += std::fabs(residual_n[iState*nbEqs+iEq]);
   }

 theta = max(theta,std::fabs(res)/max(MathTools::MathConsts::CFrealEps(), res_abs));

 }
   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState){
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
   {
    residual[iState][iEq] = theta*residual_n[iState*nbEqs+iEq] + (1.-theta)*residual_lda[iState*nbEqs+iEq];
}
  }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
