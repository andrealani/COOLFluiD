#include "STM_LDACSchemeSys.hh"
#include "Common/NotImplementedException.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "MathTools/MatrixInverter.hh"
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

MethodStrategyProvider<STM_LDACSchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSpaceTimeModule>
spaceTimeLDACSchemeSysProvider("STM_SysLDAC");

//////////////////////////////////////////////////////////////////////////////

STM_LDACSchemeSys::STM_LDACSchemeSys(const std::string& name) :
  STM_SplitterSys(name),
  _sumKmin(),
  _sumKplus(),
  _uMin(),
  _temp(),
  _uTemp(),
  _temp_mat(),
  _invK(),
  _identity()
{
}

//////////////////////////////////////////////////////////////////////////////

STM_LDACSchemeSys::~STM_LDACSchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void STM_LDACSchemeSys::setup()
{
  STM_SplitterSys::setup();
 _invK.resize(_nbEquations,_nbEquations);
  _sumKmin.resize(_nbEquations,_nbEquations);
  _sumKplus.resize(_nbEquations,_nbEquations);
  _uMin.resize(_nbEquations);
 _uTemp.resize(_nbEquations,_nbEquations);
_temp_mat.resize(_nbEquations,_nbEquations);
 _identity.resize(_nbEquations);
_temp.resize(_nbEquations);
  _identity = 1.0;
}

//////////////////////////////////////////////////////////////////////////////

void STM_LDACSchemeSys::computePicardJacob(std::vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"STM_LDASchemeScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void STM_LDACSchemeSys::distributePast(const vector<State*>& tStates)
{
  const CFuint nbStatesInCell = _nbStatesInCell;
  const CFuint dim  = PhysicalModelStack::getActive()->getDim();
  const CFreal dt = (SubSystemStatusStack::getActive()->getDT()/dim);
  const CFreal Area = _cellVolume/(dim+1.);
  const CFuint nbEqs = _nbEquations;
   DistributionData& ddata = getMethodData().getDistributionData();

  RealVector m_phi = ddata.phi;
  RealVector& past_residuals = *ddata.past_residuals;

  _sumKmin = *_kMin[0];
  _sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += *_kMin[iState];
    _sumKplus  += *_kPlus[iState];
  }

  _inverter->invert(_sumKplus, _invK);
  CFLogDebugMax( "_sumKmin = " << _sumKmin << "\n");

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _temp  = (*tStates[iState]);
    _temp *= Area ;
    _uTemp = *_kPlus[iState]*_invK;     
//if ((ddata.cellID == 1) ){

  //   CF_DEBUG_OBJ(ddata.cellID);
// CF_DEBUG_OBJ(_uTemp);

// }

    _identity = 1./6.;
    _temp_mat = _uTemp;
//if ((ddata.cellID == 1) ){

  //   CF_DEBUG_OBJ(ddata.cellID);
//CF_DEBUG_OBJ(_temp_mat);

// }

_temp_mat += _identity;
// if ((ddata.cellID == 1) ){

  //   CF_DEBUG_OBJ(ddata.cellID);
// CF_DEBUG_OBJ(_temp_mat);

// }

    _uMin = _temp_mat*_temp;
    _uMin *= -1.;

//if ((ddata.cellID == 1) ){

  //   CF_DEBUG_OBJ(ddata.cellID);
 //CF_DEBUG_OBJ(_temp_mat);

// }



    for (CFuint j=0; j< nbEqs; ++j){


        past_residuals[(iState*nbEqs)+j] = _uMin[j];

    } 

//if ((ddata.cellID == 1) ){

  //   CF_DEBUG_OBJ(ddata.cellID);
 //CF_DEBUG_OBJ(_uMin);

// }
   
    _identity = -1./12.;
    _temp_mat = _uTemp;
_temp_mat += _identity; 
    _temp   = 0.0;
    for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
    if (jState != iState){
      _temp  += (*tStates[jState]);
     }
    }
    _uMin = _temp_mat*_temp;
    _uMin *= Area;
    past_residuals.slice(iState*nbEqs, nbEqs) -= _uMin.slice(0, nbEqs);
  //  if ((ddata.cellID == 1) ){

    // CF_DEBUG_OBJ(ddata.cellID);
// CF_DEBUG_OBJ(_uMin);

// }

  _uMin = _uTemp*m_phi;
// if ((ddata.cellID == 1309) || (ddata.cellID == 1310)){
//
//       CF_DEBUG_OBJ(ddata.cellID);
// CF_DEBUG_OBJ(_uMin);
//
// }
    _uMin *= dt;
    past_residuals.slice(iState*nbEqs, nbEqs) += _uMin.slice(0, nbEqs);

//if ((ddata.cellID == 1) ){

  //   CF_DEBUG_OBJ(ddata.cellID);
// CF_DEBUG_OBJ(_uMin);

//CF_DEBUG_OBJ(m_phi);
//}

    }
  vector<RealMatrix>& betas = 
      *getMethodData().getDistributionData().currBetaMat;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations
      
      betas[iState] = *_kPlus[iState]*_invK;
    }
}


//////////////////////////////////////////////////////////////////////////////

void STM_LDACSchemeSys::distribute(vector<RealVector>& residual)
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

  _sumKmin = *_kMin[0];
  _sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKmin  += *_kMin[iState];
    _sumKplus  += *_kPlus[iState];
  }
_inverter->invert(_sumKplus, _invK);
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin = 0.0 ;
    _uTemp = 0.0 ;

    // Compute (u_k - u_0)
    _temp  = (*tStates[iState]);
    _temp *= Area ;
    _uTemp = *_kPlus[iState]*_invK;
 
    _identity = 1./6.;
_temp_mat =_uTemp;
//if ((ddata.cellID == 1) ){

  //   CF_DEBUG_OBJ(ddata.cellID);
// CF_DEBUG_OBJ(_uTemp);

// }


    _temp_mat += _identity;

//if ((ddata.cellID == 1) ){

  //   CF_DEBUG_OBJ(ddata.cellID);
// CF_DEBUG_OBJ(_temp_mat);

// }
    
_uMin = _temp_mat*_temp;
    residual[iState] += _uMin;

    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");

    _identity = -1./12.;
_temp_mat =_uTemp;
//if ((ddata.cellID == 1) ){

  //   CF_DEBUG_OBJ(ddata.cellID);
// CF_DEBUG_OBJ(_temp_mat);

// }

    _temp_mat += _identity;
//if ((ddata.cellID == 1) ){

  //   CF_DEBUG_OBJ(ddata.cellID);
 //CF_DEBUG_OBJ(_temp_mat);

// }

    _temp  = 0.0;
    for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
   if (jState != iState){
    _temp += (*tStates[jState]);
     }
    }
    _uMin = _temp_mat*_temp;
   _uMin *= Area;
   residual[iState] += _uMin;

    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");

    _uMin = 0.0;
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
    residual[iState] += _uMin;
//if ((ddata.cellID == 1) ){

  //   CF_DEBUG_OBJ(ddata.cellID);
// CF_DEBUG_OBJ(_uMin);

//CF_DEBUG_OBJ(m_phi);
//}




/*if ((ddata.cellID == 1309) || (ddata.cellID == 1310)){

      CF_DEBUG_OBJ(ddata.cellID);
CF_DEBUG_OBJ(_uMin);

}*/
    CFLogDebugMax( "residual[" << iState << "] = " << residual[iState] << "\n");
  }

  
_inverter->invert(_sumKplus, _invK);
    vector<RealMatrix>& betas = 
      *getMethodData().getDistributionData().currBetaMat;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations
      
      betas[iState] = *_kPlus[iState]*_invK;
    }

}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
