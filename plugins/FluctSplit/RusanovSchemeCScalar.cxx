#include "FluctSplit/RusanovSchemeCScalar.hh"
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

MethodStrategyProvider<RusanovSchemeCScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
rusanovSchemeScalarCProvider("ScalarRusanovC");

//////////////////////////////////////////////////////////////////////////////

RusanovSchemeCScalar::RusanovSchemeCScalar(const std::string& name) :
  RDS_SplitterScalar(name),
  m_alpha(),
  m_sumUmin()
{
}

//////////////////////////////////////////////////////////////////////////////

RusanovSchemeCScalar::~RusanovSchemeCScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void RusanovSchemeCScalar::setup()
{
  RDS_SplitterScalar::setup();

  m_alpha.resize(_nbEquations);
  m_sumUmin.resize(_nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void RusanovSchemeCScalar::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;
  const CFreal invnstates = 1.0 / _nbStatesInCell;

   // calculate the alpha artificial diffusion coefficient
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
     for (CFuint iEq = 0; iEq < _nbEquations; iEq++) {
        m_alpha[iEq] = std::max(m_alpha[iEq],std::abs(_k[iState][iEq]));
     }
  }

  // calculate the residual distribution
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

     // calculate the diffusion value
     m_sumUmin = (*tStates[iState] - *tStates[0]);
     for (CFuint jState = 1; jState < _nbStatesInCell; ++jState) {
       m_sumUmin += (*tStates[iState] - *tStates[jState]);
     }

     residual[iState] = m_alpha * m_sumUmin;
     residual[iState] += phiT;
     residual[iState] *= invnstates;
  }
}

//////////////////////////////////////////////////////////////////////////////

void RusanovSchemeCScalar::distributePart(vector<RealVector>& residual)
{
  throw Common::NotImplementedException (FromHere(),"RusanovSchemeCScalar::distributePart()");
}

//////////////////////////////////////////////////////////////////////////////

void RusanovSchemeCScalar::computePicardJacob(vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"RusanovSchemeCScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
