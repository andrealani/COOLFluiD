#include "FluctSplit/NSchemeSys.hh"
#include "MathTools/MatrixInverter.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitSystem.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NSchemeSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
nSchemeSysProvider("SysN");

//////////////////////////////////////////////////////////////////////////////

NSchemeSys::NSchemeSys(const std::string& name) :
  RDS_SplitterSys(name),
  _sumKminU(),
  _sumKmin(),
  _invK(),
  _uInflow(),
  _uDiff(),
  _temp(),
  _tempMat(),
  _sumKplus(),
  _betaLDA()
{
}

//////////////////////////////////////////////////////////////////////////////

NSchemeSys::~NSchemeSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeSys::setup()
{
  RDS_SplitterSys::setup();

  _sumKminU.resize(_nbEquations);
  _sumKmin.resize(_nbEquations, _nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
  _uInflow.resize(_nbEquations);
  _uDiff.resize(_nbEquations);
  _temp.resize(_nbEquations);
  _tempMat.resize(_nbEquations, _nbEquations);
  _sumKplus.resize(_nbEquations, _nbEquations);
  _betaLDA.resize(_nbEquations, _nbEquations);

}

//////////////////////////////////////////////////////////////////////////////

void NSchemeSys::distribute(vector<RealVector>& residual)
{
  // for efficiency reason iState = 0 is treated out of the loop
  // this allows to avoid resetting _sumKminU and _sumKmin to 0.0
  // after every loop. That operation would be done
  // nbTimeSteps*nbStatesPerCell*nbCells times!!
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  _sumKminU = (*_kMin[0])*(*tStates[0]);
  _sumKmin  = *_kMin[0];

  const CFuint nbStatesInCell = _nbStatesInCell;
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKminU += (*_kMin[iState])*(*tStates[iState]);
    _sumKmin  += *_kMin[iState];
  }

  _inverter->invert(_sumKmin, _invK);

  _uInflow = _invK * _sumKminU;

//   CFout << "sk_u " << _sumKminU << "\n";
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {

//     CFout << "kmin   " << *_kMin[iState] << "\n";
//     CFout << "kplus  " << *_kPlus[iState] << "\n";
//     CFout << "states " << (*tStates[iState]) << "\n";

    residual[iState] = *(_kPlus[iState]) * (*tStates[iState] - _uInflow);
  }

  if (getMethodData().getDistributionData().computeBetas) {
    _sumKplus = *_kPlus[0];
    for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
      _sumKplus += *_kPlus[iState];
    }
    // invert the sum of K+ matrix
    _inverter->invert(_sumKplus, _invK);
    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
      (*getMethodData().getDistributionData().currBetaMat)[iState] =
	(*_kPlus[iState])*_invK;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeSys::distributePart(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  
  _sumKminU.slice(0, _nbEquations) = (*_kMin[0]) * tStates[0]->slice(_firstVarID, _nbEquations);
  _sumKmin  = *_kMin[0];

  const CFuint nbStatesInCell = _nbStatesInCell;
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    _sumKminU.slice(0, _nbEquations) += (*_kMin[iState])*tStates[iState]->slice(_firstVarID, _nbEquations);
    _sumKmin  += *_kMin[iState];
  }

  _inverter->invert(_sumKmin, _invK);

  _uInflow = _invK*_sumKminU;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    residual[iState].slice(_firstVarID, _nbEquations) =
      (*_kPlus[iState])*(tStates[iState]->slice(_firstVarID, _nbEquations) - _uInflow.slice(0, _nbEquations));
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeSys::computePicardJacob(vector<RealMatrix*>& jacob)
{
  // carefull with the signs !!!
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    const CFuint nStart = iState*_nbStatesInCell;
    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      RealMatrix *const block = jacob[nStart + jState];

      _tempMat = _invK*(*_kMin[jState]);
      if (iState == jState) {
	for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
	  _tempMat(iEq,iEq) -= 1.0;
	}
      }

      _tempMat *= -1.0;

      (*block) = (*_kPlus[iState])*_tempMat;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
