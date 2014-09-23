#include "PSISchemeSys.hh"
#include "MathTools/MatrixInverter.hh"
#include "FluctSplit/FluctSplitSystem.hh"
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

MethodStrategyProvider<PSISchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
psiSchemeSysProvider("SysPSI");

//////////////////////////////////////////////////////////////////////////////

PSISchemeSys::PSISchemeSys(const std::string& name) :
  RDS_SplitterSys(name),
   _k(),
  _sumKmin(),
  _invK(),
  _sumKminU(),
  _uTemp(),
  _uDiff(),
  _phi(),
  _sumBeta(),
  _invCoeff(),
  _temp()
{
}

//////////////////////////////////////////////////////////////////////////////

PSISchemeSys::~PSISchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void PSISchemeSys::setup()
{
  RDS_SplitterSys::setup();

   _k.resize(_nbEquations, _nbEquations);
  _sumKmin.resize(_nbEquations, _nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
  _sumKminU.resize(_nbEquations);
  _uTemp.resize(_nbEquations);
  _uDiff.resize(_nbEquations);
  _phi.resize(_nbEquations);
  _sumBeta.resize(_nbEquations);
  _invCoeff.resize(_nbEquations);
  _temp.resize(_nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void PSISchemeSys::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  _k = (*_kPlus[0]) + (*_kMin[0]);
  _phi = _k * (*tStates[0]);
  _sumKminU = (*_kMin[0]) * (*tStates[0]);
  _sumKmin  = *_kMin[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _k = (*_kPlus[iState]) + (*_kMin[iState]);
    _phi += _k*(*tStates[iState]);
    _sumKminU += (*_kMin[iState])*(*tStates[iState]);
    _sumKmin  += *_kMin[iState];
  }

  _inverter->invert(_sumKmin, _invK);
  _uTemp = _invK * _sumKminU;
  _sumBeta = 0.0;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    residual[iState] =
      (*_kPlus[iState]) * (*tStates[iState] - _uTemp);

    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      residual[iState][iEq] = (std::abs(_phi[iEq]) > MathTools::MathConsts::CFrealEps()) ?
  residual[iState][iEq]/_phi[iEq] : 0.0;
      _sumBeta[iEq] += std::max<CFreal>(0., residual[iState][iEq]);
    }
  }

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _invCoeff[iEq] = (std::abs(_sumBeta[iEq]) > MathTools::MathConsts::CFrealEps()) ?
      _phi[iEq]/_sumBeta[iEq] : 0.0;
  }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      residual[iState][iEq] = std::max<CFreal>(0., residual[iState][iEq]);
    }
    residual[iState] *= _invCoeff;
  }
}

//////////////////////////////////////////////////////////////////////////////

void PSISchemeSys::distributePart(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;

  _k = (*_kPlus[0]) + (*_kMin[0]);
  _phi.slice(0, _nbEquations) = _k * tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKminU.slice(0, _nbEquations) = (*_kMin[0]) * tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKmin  = *_kMin[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _k = (*_kPlus[iState]) + (*_kMin[iState]);
    _phi.slice(0, _nbEquations) += _k * tStates[iState]->slice(_firstVarID, _nbEquations);
    _sumKminU.slice(0, _nbEquations) += (*_kMin[iState]) * tStates[iState]->slice(_firstVarID, _nbEquations);
    _sumKmin  += *_kMin[iState];
  }
  _inverter->invert(_sumKmin, _invK);
  _uTemp = _invK * _sumKminU;
  _sumBeta = 0.0;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    residual[iState].slice(_firstVarID, _nbEquations) =
      (*_kPlus[iState])*(tStates[iState]->slice(_firstVarID, _nbEquations) - _uTemp.slice(0, _nbEquations));

    for (CFuint iEq = _firstVarID, jEq = 0; iEq < _lastVarID; ++iEq, ++jEq) {
      residual[iState][iEq] = (std::abs(_phi[jEq]) > MathTools::MathConsts::CFrealEps()) ? residual[iState][iEq]/_phi[jEq] : 0.0;
      _sumBeta[jEq] += std::max<CFreal>(0., residual[iState][iEq]);
    }
  }

  for (CFuint jEq = 0; jEq < _nbEquations; ++jEq) {
    _invCoeff[jEq] = (std::abs(_sumBeta[jEq]) > MathTools::MathConsts::CFrealEps()) ? _phi[jEq]/_sumBeta[jEq] : 0.0;
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
