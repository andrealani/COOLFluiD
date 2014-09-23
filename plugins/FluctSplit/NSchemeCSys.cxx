#include "NSchemeCSys.hh"
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

MethodStrategyProvider<NSchemeCSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
ncSchemeSysProvider("SysNC");

//////////////////////////////////////////////////////////////////////////////

NSchemeCSys::NSchemeCSys(const std::string& name) :
  RDS_SplitterSys(name),
  _sumKplusU(),
  _sumKplus(),
  _invK(),
  _uInflow(),
  _uDiff(),
  _temp(),
  _tempBkp(),
  _tempMat(),
  _tmp(),
  _sumKU()
{
}

//////////////////////////////////////////////////////////////////////////////

NSchemeCSys::~NSchemeCSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCSys::setup()
{
  RDS_SplitterSys::setup();
  
  _sumKplusU.resize(_nbEquations);
  _sumKplus.resize(_nbEquations, _nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
  _uInflow.resize(_nbEquations);
  _uDiff.resize(_nbEquations);
  _temp.resize(_nbEquations);
  _tempBkp.resize(_nbEquations);
  _tempMat.resize(_nbEquations, _nbEquations);
  _tmp.resize(_nbEquations,_nbEquations);
  _sumKU.resize(_nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCSys::distribute(vector<RealVector>& residual)
{
  // AL: OLD implementation left here for the moment (performance comparison needed)
  //  _sumKplusU = (*_kPlus[0]) * (*tStates[0]);
  //   _sumKplus = *_kPlus[0];
  //   for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
  //     _sumKplusU += (*_kPlus[iState])*(*tStates[iState]);
  //     _sumKplus += *_kPlus[iState];
  //   }

  //   CFLogDebugMax( "sumKplusU = " << _sumKplusU << "\n");
  //   CFLogDebugMax( "sumKplus = " << _sumKplus << "\n");

  //   _inverter->invert(_sumKplus, _invK);

  //   CFLogDebugMax( "invK = " << "\n" <<_invK << "\n");

  //   _uInflow = _invK * (_sumKplusU - phiT);

  //   CFLogDebugMax( "uInflow = " << "\n" << _uInflow << "\n");

  //   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
  //     residual[iState] =
  //       (*_kPlus[iState])*(*tStates[iState] - _uInflow);
  //   }

  // beta LW
  //   RealMatrix sumAbsK(_nbEquations,_nbEquations);
  //   RealMatrix invSumAbsK(_nbEquations,_nbEquations);
  //   RealMatrix k(_nbEquations,_nbEquations);
  //   const CFuint dimPlus1 = PhysicalModelStack::getActive()->getDim() + 1;

  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;

  _tmp = (*_kPlus[0])  + (*_kMin[0]);

  //  for (CFuint i=0; i < tmp.size(); ++i) {
  //     sumAbsK[i] = std::abs(tmp[i]);
  //   }

  _sumKU = _tmp*(*tStates[0]);
  _sumKplusU = (*_kPlus[0]) * (*tStates[0]);
  _sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _sumKplusU += (*_kPlus[iState])*(*tStates[iState]);
    _sumKplus += *_kPlus[iState];

    _tmp = (*_kPlus[iState])  + (*_kMin[iState]);
    _sumKU += _tmp*(*tStates[iState]);

    //     for (CFuint i=0; i < tmp.size(); ++i) {
    //       sumAbsK[i] += std::abs(tmp[i]);
    //     }
  }

  //  for (CFuint iState = 0; iState < _kPlus.size(); ++iState) {
  // for (CFuint i = 0; i < _nbEquations; ++i) {
  //   (*_kPlus[iState])(i,i) += 1e-8;
  //  }
  // }

  _inverter->invert(_sumKplus, _invK);

  // _inverter->invert(sumAbsK, invSumAbsK);

  CFLogDebugMax( "invK = " << "\n" <<_invK << "\n");

  _uInflow = _invK * (_sumKplusU - _sumKU);

  CFLogDebugMax( "uInflow = " << "\n" << _uInflow << "\n");

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    residual[iState] = (*_kPlus[iState])*(*tStates[iState] - _uInflow);

    RealMatrix& betaLDA = (*getMethodData().getDistributionData().currBetaMat)[iState];
    betaLDA = (*_kPlus[iState])*_invK;

    // L_W
    // k = (*_kPlus[iState])  + (*_kMin[iState]);
    //     tmp = k*invSumAbsK;
    //     tmp *= CFL::getInstance().getCFL();
    //     for (CFuint i=0; i < _nbEquations; ++i) {
    //       tmp(i,i) += (1./dimPlus1);
    //     }

    //    RealVector verr(_nbEquations);
    //verr = sumKU - phiT;

    //  cout << "sumKU = " << sumKU << endl;
    //     cout << "phiT  = " << phiT << endl;
    //     RealVector t(_nbEquations);
    //     t = tmp*(sumKU - phiT);
    //     cout << "t  = " << t << endl << endl;

    //    if (_isPerturb) {
    //  _temp = _tempBkp;
    //}
    //else {
    //  _temp = _tmp*(_sumKU - phiT);
    //  _tempBkp = _temp;
    // }
    residual[iState] -= betaLDA*(_sumKU - phiT);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCSys::distributePart(vector<RealVector>& residual)
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

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    residual[iState].slice(_firstVarID, _nbEquations) = (*_kPlus[iState]) *
      (tStates[iState]->slice(_firstVarID, _nbEquations) - _uInflow.slice(0, _nbEquations));
  }
}
  
//////////////////////////////////////////////////////////////////////////////

void NSchemeCSys::computePicardJacob(vector<RealMatrix*>& jacob)
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
