#include "NSchemeScalar.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitScalar.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NSchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
nSchemeScalarProvider("ScalarN");

//////////////////////////////////////////////////////////////////////////////

NSchemeScalar::NSchemeScalar(const std::string& name) :
  RDS_SplitterScalar(name),
  _sumKminU(),
  _sumKmin(),
  _uTemp(),
  _uMin(),
  _temp()
{
}

//////////////////////////////////////////////////////////////////////////////

NSchemeScalar::~NSchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeScalar::setup()
{
  RDS_SplitterScalar::setup();

  _sumKminU.resize(_nbEquations);
  _sumKmin.resize(_nbEquations);
  _uTemp.resize(_nbEquations);
  _uMin.resize(_nbEquations);
  _temp.resize(_nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void NSchemeScalar::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
    
  _sumKminU = _kMin[0] * (*tStates[0]);
  _sumKmin  = _kMin[0];

  const CFuint nbStatesInCell = _nbStatesInCell;
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKminU += _kMin[iState] * (*tStates[iState]);
    _sumKmin  += _kMin[iState];
  }

  _uTemp = _sumKminU / _sumKmin;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    residual[iState] = _kPlus[iState] * (*tStates[iState] - _uTemp);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeScalar::distributePart(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;

  _sumKminU.slice(0, _nbEquations) = _kMin[0].slice(0, _nbEquations) * tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKmin = _kMin[0];

  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _sumKminU.slice(0, _nbEquations) += _kMin[iState].slice(0, _nbEquations) * tStates[iState]->slice(_firstVarID, _nbEquations);
    _sumKmin  += _kMin[iState];
  }
  _uTemp = _sumKminU/_sumKmin;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uMin.slice(0, _nbEquations) =
      tStates[iState]->slice(_firstVarID, _nbEquations) - _uTemp.slice(0, _nbEquations);

    residual[iState].slice(_firstVarID, _nbEquations) =
      _kPlus[iState].slice(0, _nbEquations) * _uMin.slice(0, _nbEquations);
    
    // residual[iState](_firstVarID) =
    // _kPlus[iState](0) * (tStates[iState](_firstVarID) - _uTemp(0));
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeScalar::computePicardJacob(vector<RealMatrix*>& jacob)
{
  // carefull with the signs !!!
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    const CFuint nStart = iState*_nbStatesInCell;
    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      RealMatrix *const block = jacob[nStart + jState];

      _temp = _kMin[jState]/_sumKmin;
      if (iState == jState) {
        _temp -= 1.0;
      }
      _uTemp = _kPlus[iState]*_temp;

      _uTemp *= -1.0;

      (*block) = _uTemp;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
