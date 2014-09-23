#include "STM_LDACSchemeScalar.hh"
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

MethodStrategyProvider<STM_LDACSchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeLDACSchemeScalarProvider("STM_ScalarLDAC");

//////////////////////////////////////////////////////////////////////////////

STM_LDACSchemeScalar::STM_LDACSchemeScalar(const std::string& name) :
  STM_SplitterScalar(name),
  _sumKmin(),
  _sumKplus(),
  _uMin(),
  _uTemp()
{
}

//////////////////////////////////////////////////////////////////////////////

STM_LDACSchemeScalar::~STM_LDACSchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void STM_LDACSchemeScalar::setup()
{
  STM_SplitterScalar::setup();

  _sumKmin.resize(_nbEquations);
  _sumKplus.resize(_nbEquations);
  _uMin.resize(_nbEquations);
  _uTemp.resize(_nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void STM_LDACSchemeScalar::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"STM_LDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void STM_LDACSchemeScalar::distributePast(const vector<State*>& tStates)
{
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint dim  = PhysicalModelStack::getActive()->getDim();
  const CFreal dt = (SubSystemStatusStack::getActive()->getDT()/dim);
  const CFreal Area = _cellVolume/(dim+1.);
  const CFuint nbEqs = _nbEquations;
   DistributionData& ddata = getMethodData().getDistributionData();

  RealVector m_phi = ddata.phi;
  RealVector& past_residuals = *ddata.past_residuals;
  
  _sumKmin = _kMin[0];
  _sumKplus = _kPlus[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += _kMin[iState];
    _sumKplus  += _kPlus[iState];
  }
  CFLogDebugMax( "_sumKmin = " << _sumKmin << "\n");

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin  = (*tStates[iState]);
    _uMin *= Area ;
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
    _uMin *= Area;
    past_residuals.slice(iState*nbEqs, nbEqs) += _uMin.slice(0, nbEqs);
    

    _uMin = ((_kPlus[iState]/_sumKmin)*ddata.phi)*dt;
    past_residuals.slice(iState*nbEqs, nbEqs) -= _uMin.slice(0, nbEqs);
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

void STM_LDACSchemeScalar::distribute(vector<RealVector>& residual)
{
  // Compute sumKmin
  DistributionData& ddata = getMethodData().getDistributionData();
  const vector<State*>& tStates = *ddata.states;
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint dim  = PhysicalModelStack::getActive()->getDim();
  const CFreal timeStep = (SubSystemStatusStack::getActive()->getDT())/dim;
  const CFreal Area = _cellVolume/(dim+1.);
   const CFuint nbEqs = _nbEquations;
   const RealVector& past_residuals = *ddata.past_residuals;
  RealVector m_phi = ddata.phi;

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

    _uTemp = _kPlus[iState]/_sumKmin;
// if ((ddata.cellID == 1309) || (ddata.cellID == 1310)){
// 
//       CF_DEBUG_OBJ(ddata.cellID);
// CF_DEBUG_OBJ(_uTemp);
// 
// }
    _uMin = _uTemp*m_phi;
// if ((ddata.cellID == 1309) || (ddata.cellID == 1310)){
// 
//       CF_DEBUG_OBJ(ddata.cellID);
// CF_DEBUG_OBJ(_uMin);
// 
// }
    _uMin *= timeStep;
    residual[iState] -= _uMin;/*
if ((ddata.cellID == 1309) || (ddata.cellID == 1310)){

      CF_DEBUG_OBJ(ddata.cellID);
CF_DEBUG_OBJ(_uMin);

}*/
    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
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
