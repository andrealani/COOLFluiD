#include "PSISchemeScalar.hh"
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

MethodStrategyProvider<PSISchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
psiSchemeScalarProvider("ScalarPSI");

//////////////////////////////////////////////////////////////////////////////

PSISchemeScalar::PSISchemeScalar(const std::string& name) :
  RDS_SplitterScalar(name),
  _sumKminU(),
  _sumKmin(),
  _uTemp(),
  _uDiff(),
  _phy(),
  _sumBeta(),
  _invCoeff(),
  _temp()
{
}

//////////////////////////////////////////////////////////////////////////////

PSISchemeScalar::~PSISchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void PSISchemeScalar::setup()
{
  RDS_SplitterScalar::setup();

  _sumKminU.resize(_nbEquations);
  _sumKmin.resize(_nbEquations);
  _uTemp.resize(_nbEquations);
  _uDiff.resize(_nbEquations);
  _phy.resize(_nbEquations);
  _sumBeta.resize(_nbEquations);
  _invCoeff.resize(_nbEquations);
  _temp.resize(_nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void PSISchemeScalar::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;

  _phy = _k[0]*(*tStates[0]);
  _sumKminU = _kMin[0]*(*tStates[0]);
  _sumKmin = _kMin[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _phy += _k[iState]*(*tStates[iState]);
    _sumKminU += _kMin[iState]*(*tStates[iState]);
    _sumKmin  += _kMin[iState];
  }
  _uTemp = _sumKminU / _sumKmin;
  _sumBeta = 0.0;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    residual[iState] = _kPlus[iState]*(*tStates[iState] - _uTemp);
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      residual[iState][iEq] = (MathChecks::isNotZero(_phy[iEq]) ? residual[iState][iEq]/_phy[iEq] : 0.0);
      _sumBeta[iEq] += std::max<CFreal>(0., residual[iState][iEq]);
    }
  }

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _invCoeff[iEq] = (MathChecks::isNotZero(_sumBeta[iEq]) ? _phy[iEq]/_sumBeta[iEq] : 0.0);
  }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      residual[iState][iEq] = std::max<CFreal>(0., residual[iState][iEq]);
    }
    residual[iState] *= _invCoeff;
  }
}

//////////////////////////////////////////////////////////////////////////////

void PSISchemeScalar::distributePart(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  
  _phy.slice(0, _nbEquations) = _k[0].slice(0, _nbEquations) * tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKminU.slice(0, _nbEquations) = _kMin[0].slice(0, _nbEquations)* tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKmin = _kMin[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _temp.slice(0, _nbEquations) = _k[iState].slice(0, _nbEquations) * tStates[iState]->slice(_firstVarID, _nbEquations);
    _phy += _temp;
    _temp.slice(0, _nbEquations) = _kMin[iState].slice(0, _nbEquations) * tStates[iState]->slice(_firstVarID, _nbEquations);
    _sumKminU += _temp;
    _sumKmin  += _kMin[iState];
  }

  _uTemp = _sumKminU/_sumKmin;
  _sumBeta = 0.0;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _uDiff.slice(0, _nbEquations) = tStates[iState]->slice(_firstVarID, _nbEquations) - _uTemp.slice(0, _nbEquations);
    residual[iState].slice(_firstVarID, _nbEquations) = _kPlus[iState].slice(0, _nbEquations) * _uDiff.slice(0, _nbEquations);

    for (CFuint iEq = _firstVarID, jEq = 0; iEq < _lastVarID; ++iEq, ++jEq) {
      residual[iState][iEq] = (std::abs(_phy[jEq]) > MathTools::MathConsts::CFrealEps()) ?
	residual[iState][iEq]/_phy[jEq] : 0.0;
      _sumBeta[jEq] += std::max<CFreal>(0., residual[iState][iEq]);
    }
  }

  for (CFuint jEq = 0; jEq < _nbEquations; ++jEq) {
    _invCoeff[jEq] = (std::abs(_sumBeta[jEq]) > MathTools::MathConsts::CFrealEps()) ?
      _phy[jEq]/_sumBeta[jEq] : 0.0;
  }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    for (CFuint iEq = _firstVarID, jEq = 0; iEq < _lastVarID; ++iEq, ++jEq) {
      residual[iState][iEq] = std::max<CFreal>(0., residual[iState][iEq]);
      residual[iState][iEq] *= _invCoeff[jEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
