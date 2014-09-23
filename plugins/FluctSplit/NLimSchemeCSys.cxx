#include "NLimSchemeCSys.hh"
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

MethodStrategyProvider<NLimSchemeCSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
nlimcSchemeSysProvider("SysNLimC");

////////////////////////////////////////////////////////////////////////////

void NLimSchemeCSys::defineConfigOptions(Config::OptionList& options)
{
  ///@todo NV: for 3D we should put as input two angles
   options.addConfigOption<CFreal>("toto","Choose a direction to project the residual");
}
//////////////////////////////////////////////////////////////////////////////

NLimSchemeCSys::NLimSchemeCSys(const std::string& name) :
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
  _sumKU(),
  _normal(),
  _eValues(0)
{
  addConfigOptionsTo(this);

  m_angle = 10.0;
  setParameter("toto",&m_angle);
}

//////////////////////////////////////////////////////////////////////////////

NLimSchemeCSys::~NLimSchemeCSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void NLimSchemeCSys::setup()
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
  _eValues.resize(_nbEquations);
  // we size the direction to 2 because 3D is not implemented for this
  _normal.resize(2);
  _phy.resize(_nbEquations);
  _phyChar.resize(_nbEquations);
  _rightEigenVector.resize(_nbEquations, _nbEquations);
  _leftEigenVector.resize(_nbEquations, _nbEquations);
  _sumBeta.resize(_nbEquations);
  
CFuint maxNbStatesInCell = _kPlus.size();

  _residualChar.resize(maxNbStatesInCell);
for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _residualChar[iState].resize(_nbEquations);
  }

  _beta.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _beta[iState].resize(_nbEquations);
  }

  _betaLim.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
    _betaLim[iState].resize(_nbEquations);
  }
  _normal[XX] = cos(m_angle*3.1415/180.);
  _normal[YY] = sin(m_angle*3.1415/180.);
}

//////////////////////////////////////////////////////////////////////////////

void NLimSchemeCSys::distribute(vector<RealVector>& residual)
{
   const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
   const RealVector& phiT = getMethodData().getDistributionData().phi;

   _sumKplusU = (*_kPlus[0]) * (*tStates[0]);
   _sumKplus = *_kPlus[0];
   for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
      _sumKplusU += (*_kPlus[iState])*(*tStates[iState]);
      _sumKplus += *_kPlus[iState];
    }

    CFLogDebugMax( "sumKplusU = " << _sumKplusU << "\n");
    CFLogDebugMax( "sumKplus = " << _sumKplus << "\n");

    _inverter->invert(_sumKplus, _invK);

    CFLogDebugMax( "invK = " << "\n" <<_invK << "\n");

    _uInflow = _invK * (_sumKplusU - phiT);

    CFLogDebugMax( "uInflow = " << "\n" << _uInflow << "\n");

    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      residual[iState] =
        (*_kPlus[iState])*(*tStates[iState] - _uInflow);
    }


  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    residual[iState] = (*_kPlus[iState])*(*tStates[iState] - _uInflow);

    RealMatrix& betaLDA = (*getMethodData().getDistributionData().currBetaMat)[iState];
    betaLDA = (*_kPlus[iState])*_invK;

 
  }

  /// Compute the cell residual...
  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    _phy[iEq] = 0.0;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      _phy[iEq] += residual[iState][iEq];
    }
  }

  /// Transform into characteristic variables the residual, phy
  // Compute the EigenVectors
  getMethodData().getDistribVar()->computeEigenValuesVectors(_rightEigenVector,
							 _leftEigenVector,
							 _eValues,
							 _normal);

  //Transform the residual into characteristic variables
  _phyChar = _leftEigenVector *_phy;

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

      _residualChar[iState] = _leftEigenVector * residual[iState];

  }

  // Do the limiting just as in scalar case
  

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    if (std::abs(_phyChar[iEq])>MathTools::MathConsts::CFrealEps()){
      _sumBeta[iEq] = 0.;
      for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

        _beta[iState][iEq] = _residualChar[iState][iEq]/_phyChar[iEq];

        _sumBeta[iEq] += max(0.,_beta[iState][iEq]);

      }
    }
  }

  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
    if (std::abs(_phyChar[iEq])>MathTools::MathConsts::CFrealEps()){
      for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
        _betaLim[iState][iEq] = max(0.,_beta[iState][iEq])/_sumBeta[iEq];
        _residualChar[iState][iEq] = _betaLim[iState][iEq] * _phyChar[iEq];
        }
    }
    else{
      for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
        _residualChar[iState][iEq] = 0.;
      }
    }
  }

  //Transform the residual back into conservative variables
 _phy = _rightEigenVector * _phyChar;
 for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
     residual[iState] = _rightEigenVector * _residualChar[iState];
 }
}


//////////////////////////////////////////////////////////////////////////////

void NLimSchemeCSys::distributePart(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;

  _sumKplusU.slice(0, _nbEquations ) = (*_kPlus[0]) * tStates[0]->slice(_firstVarID, _nbEquations);
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

void NLimSchemeCSys::computePicardJacob(vector<RealMatrix*>& jacob)
{
   throw Common::NotImplementedException (FromHere(),getName() + "::ComputePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
