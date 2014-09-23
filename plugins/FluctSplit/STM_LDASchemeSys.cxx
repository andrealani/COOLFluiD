#include "STM_LDASchemeSys.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "MathTools/MatrixInverter.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STM_LDASchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeLDASchemeSysProvider("STM_SysLDA");

//////////////////////////////////////////////////////////////////////////////

STM_LDASchemeSys::STM_LDASchemeSys(const std::string& name) :
  STM_SplitterSys(name),
  _sumKmin(),
  _sumKplus(),
  _invK(),
  _uMin(),
  _tempMat(),
  _phyt(),
  _phys(),
  _identity()
{
}

//////////////////////////////////////////////////////////////////////////////

STM_LDASchemeSys::~STM_LDASchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void STM_LDASchemeSys::setup()
{
  STM_SplitterSys::setup();

  _sumKmin.resize(_nbEquations,_nbEquations);
  _sumKplus.resize(_nbEquations,_nbEquations);
  _invK.resize(_nbEquations,_nbEquations);
  _uMin.resize(_nbEquations);
  _tempMat.resize(_nbEquations,_nbEquations);
  _phys.resize(_nbEquations);
  _identity.resize(_nbEquations);
  _identity = 1.0;
}


//////////////////////////////////////////////////////////////////////////////

void STM_LDASchemeSys::distributePast(const vector<State*>& tStates)
{
  // unused //  const CFreal Area = _cellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+2.);
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint nbEqs = _nbEquations;
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();
 
  RealVector& past_residuals = *getMethodData().getDistributionData().past_residuals;
 
  _phyt.resize(_nbStatesInCell);
  for(CFuint i = 0; i < _nbStatesInCell; ++i){
    _phyt[i].resize(_nbEquations);
 }

  _sumKmin  = *_kMin[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += *_kMin[iState];
  }
  _inverter->invert(_sumKmin, _invK);

  _phys = 0.0;
  CFreal coef1 = 1./3.;
  CFreal coef2 = 1./2.;


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
    _phyt[iState] *= pastArea;
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
      past_residuals[(iState*nbEqs)+j] = _uMin[j];
    }

    _identity = 1./12.;
    _tempMat = *_kPlus[iState] * _invK;
    _tempMat += _identity;
    
    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      if (jState != iState){
        for (CFuint j=0; j< nbEqs; ++j){
          _uMin = _tempMat * _phyt[jState];
          past_residuals[(iState*nbEqs)+j] -= _uMin[j];
        }
      }
    }
  }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _tempMat = *_kPlus[iState] * _invK;
    for (CFuint j=0; j< nbEqs; ++j){
      _uMin = _tempMat * _phys;
      past_residuals[(iState*nbEqs)+j] -= _uMin[j];
    }
    CFLogDebugMax( "residual[" << iState << "] = " << past_residuals[iState] << "\n");
  }

}

//////////////////////////////////////////////////////////////////////////////

void STM_LDASchemeSys::distribute(vector<RealVector>& residual)
{
  // for efficiency reason iState = 0 is treated out of the loop
  // this allows to avoid resetting _sumKminU and _sumKmin to 0.0
  // after every loop (operation that would be done
  // nbTimeSteps*nbStatesPerCell*nbCells times!!)
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;

  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();
  const CFreal Area = _cellVolume*(1./(PhysicalModelStack::getActive()->getDim()+2.)) ;
  // unused //  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFuint nbEqs = _nbEquations;

  const RealVector& past_residuals = *getMethodData().getDistributionData().past_residuals;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
      for (CFuint jEq = 0; jEq <  nbEqs; ++jEq)
        residual[iState][jEq] = past_residuals[iState*nbEqs + jEq];

   _phyt.resize(_nbStatesInCell);
   for(CFuint i = 0; i < _nbStatesInCell; ++i){
     _phyt[i].resize(_nbEquations);
   }

   _sumKmin  = *_kMin[0];
   _sumKplus  = *_kPlus[0];
   for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
     _sumKmin  += *_kMin[iState];
     _sumKplus  += *_kPlus[iState];
   }

   _inverter->invert(_sumKmin, _invK);

   _phys = 0.0;
   CFreal coef1 = 1./3.;
   CFreal coef2 = 1./2.;
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
     _phyt[iState] *= Area;

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
    residual[iState] += _tempMat * _phyt[iState];

    _identity = 1./12.;
    _tempMat = *_kPlus[iState] * _invK;
    _tempMat += _identity;
    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      if (jState != iState){
        residual[iState] -= _tempMat * _phyt[jState];
      }
    }
  }


  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _tempMat = *_kPlus[iState] * _invK;
    residual[iState] -= (_tempMat * _phys);
    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
  }

    vector<RealMatrix>& betas =
      *getMethodData().getDistributionData().currBetaMat;

    _inverter->invert(_sumKplus, _invK);

    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations

      RealMatrix& betaLDA = betas[iState];
      betaLDA = (*_kPlus[iState])*_invK;
    }


}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
