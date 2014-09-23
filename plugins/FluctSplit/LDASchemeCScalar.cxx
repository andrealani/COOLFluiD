#include "LDASchemeCScalar.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitScalar.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LDASchemeCScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
ldacSchemeScalarProvider("ScalarLDAC");

//////////////////////////////////////////////////////////////////////////////

LDASchemeCScalar::LDASchemeCScalar(const std::string& name) :
  RDS_SplitterScalar(name),
  m_sumKplus(),
  m_uTemp()
{
}

//////////////////////////////////////////////////////////////////////////////

LDASchemeCScalar::~LDASchemeCScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeCScalar::setup()
{
  RDS_SplitterScalar::setup();

  m_sumKplus.resize(_nbEquations);
  m_uTemp.resize(_nbEquations);
}


//////////////////////////////////////////////////////////////////////////////

void LDASchemeCScalar::distribute(vector<RealVector>& residual)
{
  m_sumKplus = _kPlus[0];

  const CFuint nbStatesInCell = _nbStatesInCell;
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    m_sumKplus += _kPlus[iState];
  }
  
  const RealVector& phiT = getMethodData().getDistributionData().phi;
  if (!getMethodData().getDistributionData().computeBetas) {
    m_uTemp = phiT / m_sumKplus;
    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
      residual[iState] = _kPlus[iState]*m_uTemp;
    }
  }
  else {
    vector<RealMatrix>& betas = *getMethodData().getDistributionData().currBetaMat;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
      // AL: compute betas on the fly is more expensive because
      // it involves a matrix*matrix and a matrix*vector
      // while normal scheme requires only two matrix*vector operations
      
      betas[iState] = _kPlus[iState]/m_sumKplus;
      residual[iState] = betas[iState]*phiT;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeCScalar::distributePart(vector<RealVector>& residual)
{
  m_sumKplus = _kPlus[0];
  
  const CFuint nbStatesInCell = _nbStatesInCell;
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    m_sumKplus += _kPlus[iState];
  }
  RealVector& phiT = getMethodData().getDistributionData().phi;
  m_uTemp = phiT.slice(_firstVarID, _nbEquations) / m_sumKplus;
  
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    residual[iState].slice(_firstVarID, _nbEquations) =
      _kPlus[iState].slice(0, _nbEquations)*m_uTemp.slice(0, _nbEquations);
  }
}

//////////////////////////////////////////////////////////////////////////////

void LDASchemeCScalar::computePicardJacob(vector<RealMatrix*>& jacob)
{
    // carefull with the signs !!!
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    // beta coefficient
    CFreal _beta = _kPlus[iState][0]/m_sumKplus[0];
    const CFuint nStart = iState*_nbStatesInCell;

    for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
      RealMatrix *const block = jacob[nStart + jState];

       _k[0] =  _kPlus[jState][0] + _kMin[jState][0];
      (*block) = _beta*_k[0];
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
