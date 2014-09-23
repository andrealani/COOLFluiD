#include "STM_LDASchemeScalar.hh"
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

MethodStrategyProvider<STM_LDASchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeLDASchemeScalarProvider("STM_ScalarLDA");

//////////////////////////////////////////////////////////////////////////////

STM_LDASchemeScalar::STM_LDASchemeScalar(const std::string& name) :
  STM_SplitterScalar(name),
  _sumKmin(),
  _sumKplus(),
  _uMin(),
  _uTemp()
{
}

//////////////////////////////////////////////////////////////////////////////

STM_LDASchemeScalar::~STM_LDASchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void STM_LDASchemeScalar::setup()
{
  STM_SplitterScalar::setup();

  _sumKmin.resize(_nbEquations);
  _sumKplus.resize(_nbEquations);
  _uMin.resize(_nbEquations);
  _uTemp.resize(_nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void STM_LDASchemeScalar::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"STM_LDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void STM_LDASchemeScalar::distributePast(const vector<State*>& tStates)
{
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();
  const CFreal Area = _cellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFuint nbEqs = _nbEquations;
 
  RealVector& past_residuals = *getMethodData().getDistributionData().past_residuals;
  
  _sumKmin = _kMin[0];
  _sumKplus = _kPlus[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += _kMin[iState];
    _sumKplus  += _kPlus[iState];
  }
  CFLogDebugMax( "_sumKmin = " << _sumKmin << "\n");

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin  = (*tStates[iState]);
    _uMin *= pastArea ;
    _uTemp = _kPlus[iState]/_sumKmin;
    _uMin *= (1./6. - _uTemp);
    _uMin *= -1.;
    for (CFuint j=0; j< nbEqs; ++j){
        past_residuals[(iState*nbEqs)+j] = _uMin[j];
    }

    _uTemp += 1./12.;
    _uMin   = 0.0;
    for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
    if (jState != iState){
      _uMin  += (*tStates[jState]);
      }
    }
    _uMin *= _uTemp;
    _uMin *= pastArea;
    past_residuals.slice(iState*nbEqs, nbEqs) += _uMin.slice(0, nbEqs);

    _uMin = 0.0;
    for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
      _uMin += _k[jState]*(*tStates[jState]);
     }
    _uTemp = _kPlus[iState]/_sumKmin;
    _uMin *= _uTemp;
    _uMin *= timeStep;
    _uMin /= 2.;
    past_residuals.slice(iState*nbEqs, nbEqs) -= _uMin.slice(0, nbEqs);
    }

  CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]);
    _uMin /= 2.;
    _uMin *= (Area-pastArea);
    past_residuals.slice(iState*nbEqs, nbEqs) += _uMin.slice(0, nbEqs);
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

void STM_LDASchemeScalar::distribute(vector<RealVector>& residual)
{
  // Compute sumKmin
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();
  const CFreal Area = _cellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
   const CFuint nbEqs = _nbEquations;
   const RealVector& past_residuals = *getMethodData().getDistributionData().past_residuals;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
      for (CFuint jEq = 0; jEq <  nbEqs; ++jEq)
        residual[iState][jEq] = past_residuals[iState*nbEqs + jEq];

  _sumKmin = _kMin[0];
  _sumKplus = _kPlus[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += _kMin[iState];
    _sumKplus  += _kPlus[iState];
  }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = 0.0 ;
    _uTemp = 0.0 ;

    // Compute (u_k - u_0)
    _uMin  = (*tStates[iState]);
    _uMin *= Area ;
    _uTemp = _kPlus[iState]/_sumKmin;
    _uMin *= (1./6. - _uTemp);
    residual[iState] += _uMin;
    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");

    _uTemp += 1./12.;
    _uMin  = 0.0;
    for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
    if (jState != iState){
      _uMin  += (*tStates[jState]);
      }
    }
    _uMin *= _uTemp;
    _uMin *= Area;
    residual[iState] -= _uMin;
    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");

    _uMin = 0.0;
    for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
      _uMin  += (_k[jState] * (*tStates[jState]));
     }
    _uTemp = _kPlus[iState]/_sumKmin;
    _uMin *= _uTemp;
    _uMin *= timeStep;
    _uMin /= 2.;
    residual[iState] -= _uMin;
    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
  }

  CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = *tStates[iState];
    _uMin /= 2.;
    _uMin *= (Area-pastArea);
    residual[iState] -= _uMin;
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
