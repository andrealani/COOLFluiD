#include "LDASchemeScalar.hh"
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

MethodStrategyProvider<LDASchemeScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
ldaSchemeScalarProvider("ScalarLDA");

//////////////////////////////////////////////////////////////////////////////

LDASchemeScalar::LDASchemeScalar(const std::string& name) :
  RDS_SplitterScalar(name),
  _phi(),
  _sumKplus(),
  _uTemp(),
  _temp()
{
}

//////////////////////////////////////////////////////////////////////////////

LDASchemeScalar::~LDASchemeScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeScalar::setup()
{
  RDS_SplitterScalar::setup();

  _phi.resize(_nbEquations);
  _sumKplus.resize(_nbEquations);
  _uTemp.resize(_nbEquations);
  _temp.resize(_nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeScalar::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  _phi = _k[0]*(*tStates[0]);
  _sumKplus = _kPlus[0];

  const CFuint nbStatesInCell = _nbStatesInCell;
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _phi += _k[iState]*(*tStates[iState]);
    _sumKplus += _kPlus[iState];
  }
  _uTemp = _phi / _sumKplus;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    residual[iState] = _kPlus[iState]*_uTemp;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeScalar::distributePart(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  _phi.slice(0, _nbEquations) = _k[0].slice(0, _nbEquations) * tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKplus = _kPlus[0];

  const CFuint nbStatesInCell = _nbStatesInCell;
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _phi.slice(0, _nbEquations) += _k[iState].slice(0, _nbEquations) * tStates[iState]->slice(_firstVarID, _nbEquations);
    _sumKplus += _kPlus[iState];
  }
  _uTemp = _phi/_sumKplus;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    residual[iState].slice(_firstVarID, _nbEquations) =
      _kPlus[iState].slice(0, _nbEquations)*_uTemp.slice(0, _nbEquations);
  }
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeScalar::computePicardJacob(vector<RealMatrix*>& jacob)
{
  // carefull with the signs !!!
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    // beta coefficient
    _uTemp = _kPlus[iState]/_sumKplus;
    const CFuint nStart = iState*_nbStatesInCell;

    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      RealMatrix *const block = jacob[nStart + jState];

      _temp = _uTemp*_k[jState];
      (*block) = _temp;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
