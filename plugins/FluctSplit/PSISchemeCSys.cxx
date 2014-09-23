#include "PSISchemeCSys.hh"
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

MethodStrategyProvider<PSISchemeCSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
psicSchemeSysProvider("SysPSIC");

//////////////////////////////////////////////////////////////////////////////

PSISchemeCSys::PSISchemeCSys(const std::string& name) :
  RDS_SplitterSys(name),
  _sumKplus(),
  _invK(),
  _sumKplusU(),
  _uInflow(),
  _uDiff(),
  _sumBeta(),
  _invCoeff(),
  _temp()
{
}

//////////////////////////////////////////////////////////////////////////////

PSISchemeCSys::~PSISchemeCSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void PSISchemeCSys::setup()
{
  RDS_SplitterSys::setup();

  _sumKplus.resize(_nbEquations, _nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
  _sumKplusU.resize(_nbEquations);
  _uInflow.resize(_nbEquations);
  _uDiff.resize(_nbEquations);
  _sumBeta.resize(_nbEquations);
  _invCoeff.resize(_nbEquations);
  _temp.resize(_nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void PSISchemeCSys::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;

  _sumKplusU = (*_kPlus[0])*(*tStates[0]);
  _sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _sumKplusU += (*_kPlus[iState])*(*tStates[iState]);
    _sumKplus  += *_kPlus[iState];
  }

  _inverter->invert(_sumKplus, _invK);
  _uInflow = _invK * (_sumKplusU - phiT);
  _sumBeta = 0.0;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    residual[iState] =
      (*_kPlus[iState]) * (*tStates[iState] - _uInflow);

    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      residual[iState][iEq] = (std::abs(phiT[iEq]) > MathTools::MathConsts::CFrealEps()) ?
	residual[iState][iEq]/phiT[iEq] : 0.0;
      _sumBeta[iEq] += std::max<CFreal>(0., residual[iState][iEq]);
    }
  }

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _invCoeff[iEq] = (std::abs(_sumBeta[iEq]) > MathTools::MathConsts::CFrealEps()) ?
      phiT[iEq]/_sumBeta[iEq] : 0.0;
  }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      residual[iState][iEq] = std::max<CFreal>(0., residual[iState][iEq]);
    }
    residual[iState] *= _invCoeff;
  }
}

//////////////////////////////////////////////////////////////////////////////

void PSISchemeCSys::distributePart(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;
  
  _sumKplusU.slice(0, _nbEquations) = (*_kPlus[0]) * tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKplus  = *_kPlus[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _sumKplusU.slice(0, _nbEquations) +=
      (*_kPlus[iState]) * tStates[iState]->slice(_firstVarID, _nbEquations);
    _sumKplus  += *_kPlus[iState];
  }

  _inverter->invert(_sumKplus, _invK);

  RealVector& phi = const_cast<RealVector&>(phiT);

  _sumKplusU.slice(0, _nbEquations) -= phi.slice(_firstVarID, _nbEquations);
  _uInflow = _invK * _sumKplusU;
  _sumBeta = 0.0;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    residual[iState].slice(_firstVarID, _nbEquations) = (*_kPlus[iState]) *
      (tStates[iState]->slice(_firstVarID, _nbEquations) -  _uInflow.slice(0, _nbEquations));

    for (CFuint iEq = _firstVarID, jEq = 0; iEq < _lastVarID; ++iEq, ++jEq) {
      residual[iState][iEq] = (std::abs(phiT[jEq]) > MathTools::MathConsts::CFrealEps()) ?
	residual[iState][iEq]/phiT[jEq] : 0.0;
      _sumBeta[jEq] += std::max<CFreal>(0., residual[iState][iEq]);
    }
  }

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _invCoeff[iEq] = (std::abs(_sumBeta[iEq]) > MathTools::MathConsts::CFrealEps()) ?
      phiT[iEq]/_sumBeta[iEq] : 0.0;
  }

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    for (CFuint iEq = _firstVarID, jEq = 0; iEq < _lastVarID; iEq++, jEq++) {
      residual[iState][iEq]  = std::max<CFreal>(0., residual[iState][iEq]);
      residual[iState][iEq] *= _invCoeff[jEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
