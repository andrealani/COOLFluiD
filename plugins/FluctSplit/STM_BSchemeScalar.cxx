#include "STM_BSchemeScalar.hh"
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

MethodStrategyProvider<STM_BSchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeBSchemeScalarProvider("STM_ScalarB");

//////////////////////////////////////////////////////////////////////////////

STM_BSchemeScalar::STM_BSchemeScalar(const std::string& name) :
  STM_SplitterScalar(name),
  _sumKmin(),
  _sumKplus(),
  _sumKminU(),
  _uInflow(),
  _uMin(),
  _uTemp()
{
}

//////////////////////////////////////////////////////////////////////////////

STM_BSchemeScalar::~STM_BSchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void STM_BSchemeScalar::setup()
{
  STM_SplitterScalar::setup();

  _sumKmin.resize(_nbEquations);
  _sumKplus.resize(_nbEquations);
  _sumKminU.resize(_nbEquations);
  _uInflow.resize(_nbEquations);
  _uMin.resize(_nbEquations);
  _uTemp.resize(_nbEquations);
 

}

//////////////////////////////////////////////////////////////////////////////

void STM_BSchemeScalar::distributePast(const vector<State*>& tStates)
{
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();
  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFuint nbEqs = _nbEquations;

  RealVector residual_lda;
  residual_lda.resize(nbStatesInCell*nbEqs);

  RealVector residual_n;
  residual_n.resize(nbStatesInCell*nbEqs);

  CFLogDebugMax( "============== TEMPORAL PART ============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]);
    _uMin *= -pastArea;
    for (CFuint j=0; j< nbEqs; ++j){
        residual_n[iState*nbEqs + j] = _uMin[j];
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
    for (CFuint j=0; j< nbEqs; ++j){
    residual_n[iState*nbEqs + j] += _uMin[j];
  
    CFLogDebugMax( "residual_order1[" << iState*nbEqs+j << "] = " << residual_n[iState*nbEqs+j] << "\n");
 } }
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin  = (*tStates[iState]);
    _uMin *= pastArea ;
    _uTemp = _kPlus[iState]/_sumKmin;
    _uMin *= (1./6. - _uTemp);
    _uMin *= -1.;
    for (CFuint j=0; j< nbEqs; ++j){
        residual_lda[(iState*nbEqs)+j] = _uMin[j];
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
    residual_lda.slice(iState*nbEqs, nbEqs) += _uMin.slice(0, nbEqs);
  
    _uMin = 0.0;
    for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
      _uMin += _k[jState]*(*tStates[jState]);
     }  
    _uTemp = _kPlus[iState]/_sumKmin;
    _uMin *= _uTemp;
    _uMin *= timeStep;
    _uMin /= 2.;
    residual_lda.slice(iState*nbEqs, nbEqs) -= _uMin.slice(0, nbEqs);
    }
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState){
    for (CFuint jEq = 0; jEq < nbEqs ; ++jEq){
      
        
      (*getMethodData().getDistributionData().past_residuals)[iState*nbEqs + jEq] = residual_lda[iState*nbEqs + jEq];
     (*getMethodData().getDistributionData().past_residuals_order1)[iState*nbEqs + jEq] = residual_n[iState*nbEqs + jEq];
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

void STM_BSchemeScalar::distribute(vector<RealVector>& residual)
{
  // Compute sumKmin
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFreal timeStep = SubSystemStatusStack::getActive()->getDT();
  const CFreal Area = _cellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
  const CFreal pastArea = _pastCellVolume/(PhysicalModelStack::getActive()->getDim()+1.);
   const CFuint nbEqs = _nbEquations;
 RealVector residual_lda;
  residual_lda.resize(nbStatesInCell*nbEqs);

  RealVector residual_n;
  residual_n.resize(nbStatesInCell*nbEqs);
   CFreal res;
  CFreal res_abs;
  CFreal theta;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
       for (CFuint jEq = 0; jEq <  nbEqs; ++jEq){
         residual_lda[iState*nbEqs + jEq] = (*getMethodData().getDistributionData().past_residuals)[iState*nbEqs + jEq];
          residual_n[iState*nbEqs + jEq] = (*getMethodData().getDistributionData().past_residuals_order1)[iState*nbEqs + jEq];
  }

 CFLogDebugMax( "============== TEMPORAL PART ============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = (Area*(*tStates[iState]));
    for (CFuint j=0; j< nbEqs; ++j)
    residual_n[iState*nbEqs+j] += _uMin[j];
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
    for (CFuint j=0; j< nbEqs; ++j){
    residual_n[iState*nbEqs+j] += _uMin[j];
    CFLogDebugMax( "residualorder1[" << iState << "] = " << residual_n[iState*nbEqs+j] << "\n");
}
  }


  CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = (*tStates[iState]);
    _uMin /= 2.;
    _uMin *= (Area - pastArea);
  for (CFuint j=0; j< nbEqs; ++j)
    residual_n[iState*nbEqs + j] -= _uMin[j];
  }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = 0.0 ;
    _uTemp = 0.0 ;

    // Compute (u_k - u_0)
    _uMin  = (*tStates[iState]);
    _uMin *= Area ;
    _uTemp = _kPlus[iState]/_sumKmin;
    _uMin *= (1./6. - _uTemp);
  for (CFuint j=0; j< nbEqs; ++j)
    residual_lda[iState*nbEqs + j] += _uMin[j];
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

  for (CFuint j=0; j< nbEqs; ++j)
    residual_lda[iState*nbEqs + j] -= _uMin[j];
    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");

    _uMin = 0.0;
    for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
      _uMin  += (_k[jState] * (*tStates[jState]));
     }
    _uTemp = _kPlus[iState]/_sumKmin;
    _uMin *= _uTemp;
    _uMin *= timeStep;
    _uMin /= 2.;
  for (CFuint j=0; j< nbEqs; ++j)
    residual_lda[iState*nbEqs + j] -= _uMin[j];
    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
  }

  CFLogDebugMax( "============== DEFORMING MESH PART============ " << "\n");
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = *tStates[iState];
    _uMin /= 2.;
    _uMin *= (Area-pastArea);
   for (CFuint j=0; j< nbEqs; ++j)
    residual_lda[iState*nbEqs + j] -= _uMin[j];
  }
for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
   { 
    res = residual_n[iEq];
    res_abs = std::fabs(residual_n[iEq]);
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState){
    res += residual_n[iState*nbEqs+iEq];
    res_abs += std::fabs(residual_n[iState*nbEqs+iEq]);
  }
theta = std::fabs(res)/max(MathTools::MathConsts::CFrealEps(), res_abs);
for (CFuint iState = 0; iState < _nbStatesInCell; ++iState){
  residual[iState][iEq] = theta*residual_n[iState*nbEqs+iEq] + (1.-theta)*residual_lda[iState*nbEqs+iEq];
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
