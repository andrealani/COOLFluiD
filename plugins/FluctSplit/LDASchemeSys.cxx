#include "LDASchemeSys.hh"
#include "MathTools/MatrixInverter.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitSystem.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LDASchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
ldaSchemeSysProvider("SysLDA");

//////////////////////////////////////////////////////////////////////////////

LDASchemeSys::LDASchemeSys(const std::string& name) :
  RDS_SplitterSys(name),
  _phi(),
  _uTemp(),
  _temp(),
  _sumKplus(),
  _k(),
  _invK(),
  _beta()
{
}

//////////////////////////////////////////////////////////////////////////////

LDASchemeSys::~LDASchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeSys::setup()
{
  RDS_SplitterSys::setup();

  _phi.resize(_nbEquations);
  _uTemp.resize(_nbEquations);
  _temp.resize(_nbEquations);
  _sumKplus.resize(_nbEquations, _nbEquations);
  _k.resize(_nbEquations, _nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
  _beta.resize(_nbEquations, _nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeSys::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;

  _k = (*_kPlus[0]) + (*_kMin[0]);
  _phi = _k * (*tStates[0]);
  _sumKplus = *_kPlus[0];

  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _k = (*_kPlus[iState]) + (*_kMin[iState]);
    _phi += _k*(*tStates[iState]);
    _sumKplus += *_kPlus[iState];
  }

//   cout << _sumKplus << "\n";

  _inverter->invert(_sumKplus, _invK);
  
  _uTemp = _invK*_phi;
  

  for (CFuint iState = 0; iState < _nbStatesInCell; iState++) {
    residual[iState] = (*_kPlus[iState])*_uTemp;
    
    if (getMethodData().getDistributionData().computeBetas) {
      (*getMethodData().getDistributionData().currBetaMat)[iState] = 
	(*_kPlus[iState])*_invK;
    }
  }
  
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeSys::distributePart(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  
  _k = (*_kPlus[0]) + (*_kMin[0]);
  _phi.slice(0, _nbEquations) = _k * tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _k = (*_kPlus[iState]) + (*_kMin[iState]);
    _phi.slice(0, _nbEquations) += _k * tStates[iState]->slice(_firstVarID, _nbEquations);
    _sumKplus += *_kPlus[iState];
  }

  _inverter->invert(_sumKplus, _invK);
  _uTemp = _invK * _phi;

  for (CFuint iState = 0; iState < _nbStatesInCell; iState++) {
    residual[iState].slice(_firstVarID, _nbEquations) =
      (*_kPlus[iState])*_uTemp.slice(0, _nbEquations);
  }
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeSys::computePicardJacob(vector<RealMatrix*>& jacob)
{
  // carefull with the signs !!!
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    // beta coefficient
    _beta = (*_kPlus[iState])*_invK;
    const CFuint nStart = iState*_nbStatesInCell;

    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      RealMatrix *const block = jacob[nStart + jState];

      _k =  (*_kPlus[jState]) + (*_kMin[jState]);
      (*block) = _beta*_k;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
