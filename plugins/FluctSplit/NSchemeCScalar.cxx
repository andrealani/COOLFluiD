#include "NSchemeCScalar.hh"
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

MethodStrategyProvider<NSchemeCScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
ncSchemeScalarProvider("ScalarNC");

//////////////////////////////////////////////////////////////////////////////

NSchemeCScalar::NSchemeCScalar(const std::string& name) :
  RDS_SplitterScalar(name),
  m_sumKplusU(),
  m_sumKplus(),
  m_uTemp(),
  m_uMin(),
  m_temp()
{
}

//////////////////////////////////////////////////////////////////////////////

NSchemeCScalar::~NSchemeCScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCScalar::setup()
{
  RDS_SplitterScalar::setup();

  m_sumKplusU.resize(_nbEquations);
  m_sumKplus.resize(_nbEquations);
  m_uTemp.resize(_nbEquations);
  m_uMin.resize(_nbEquations);
  m_temp.resize(_nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCScalar::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;
  
  m_sumKplusU = _kPlus[0] * (*tStates[0]);
  m_sumKplus  = _kPlus[0];

  const CFuint nbStatesInCell = _nbStatesInCell;
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    m_sumKplusU += _kPlus[iState] * (*tStates[iState]);
    m_sumKplus  += _kPlus[iState];
  }
  m_temp = ( m_sumKplusU - phiT )/ m_sumKplus;

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    residual[iState] = _kPlus[iState]*(*tStates[iState] - m_temp); 
    
    if (getMethodData().getDistributionData().computeBetas) {
      (*getMethodData().getDistributionData().currBetaMat)[iState] = 
	_kPlus[iState]/m_sumKplus;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCScalar::distributePart(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;
  
  m_sumKplusU.slice(0, _nbEquations) = _kPlus[0].slice(0, _nbEquations) * tStates[0]->slice(_firstVarID, _nbEquations);
  m_sumKplus = _kPlus[0];

  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    m_sumKplusU.slice(0, _nbEquations) += _kPlus[iState].slice(0, _nbEquations) * tStates[iState]->slice(_firstVarID, _nbEquations);
    m_sumKplus  += _kPlus[iState];
  }
  
  RealVector& phi = const_cast<RealVector&>(phiT);
  
  m_temp.slice(0, _nbEquations) = (m_sumKplusU.slice(0, _nbEquations) - phi.slice(_firstVarID, _nbEquations));
  m_uTemp = m_temp / m_sumKplus;
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    m_uMin.slice(0, _nbEquations) =
      tStates[iState]->slice(_firstVarID, _nbEquations) - m_uTemp.slice(0, _nbEquations);
    
    residual[iState].slice(_firstVarID, _nbEquations) =
      _kPlus[iState].slice(0, _nbEquations) * m_uMin.slice(0, _nbEquations);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCScalar::computePicardJacob(vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"NSchemeCScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
