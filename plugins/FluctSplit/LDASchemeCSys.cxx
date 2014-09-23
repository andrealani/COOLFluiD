#include "FluctSplit/LDASchemeCSys.hh"
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

MethodStrategyProvider<LDASchemeCSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
ldacSchemeSysProvider("SysLDAC");

//////////////////////////////////////////////////////////////////////////////

LDASchemeCSys::LDASchemeCSys(const std::string& name) :
  RDS_SplitterSys(name),
  _sumKplus(),
  _invK(),
  _uTemp()
{
}

//////////////////////////////////////////////////////////////////////////////

LDASchemeCSys::~LDASchemeCSys()
{
}


//////////////////////////////////////////////////////////////////////////////

void LDASchemeCSys::setup()
{
  RDS_SplitterSys::setup();

  _uTemp.resize(_nbEquations);
  _sumKplus.resize(_nbEquations, _nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
  _beta.resize(_nbEquations, _nbEquations);
  _k.resize(_nbEquations, _nbEquations);


}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeCSys::distribute(vector<RealVector>& residual)
{
  DistributionData& ddata = getMethodData().getDistributionData();

  const RealVector& phiT = ddata.phi;
  
  _sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState)
  {
    _sumKplus  += *_kPlus[iState];
  }

  _inverter->invert(_sumKplus, _invK);

  _uTemp = _invK * phiT;
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
  {
    residual[iState] = (*_kPlus[iState]) * _uTemp;

    if (ddata.computeBetas)
    {
      (*ddata.currBetaMat)[iState] = (*_kPlus[iState]) * _invK;
    }
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void LDASchemeCSys::distributePart(vector<RealVector>& residual)
{
  _sumKplus = *_kPlus[0];
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    _sumKplus  += *_kPlus[iState];
  }
  
  _inverter->invert(_sumKplus, _invK);
  
  RealVector& phi = getMethodData().getDistributionData().phi;
  _uTemp.slice(0, _nbEquations) = _invK * phi.slice(_firstVarID, _nbEquations);
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    residual[iState].slice(_firstVarID, _nbEquations) =
      (*_kPlus[iState]) * _uTemp.slice(0, _nbEquations);
  }
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeCSys::computePicardJacob(vector<RealMatrix*>& jacob)
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
