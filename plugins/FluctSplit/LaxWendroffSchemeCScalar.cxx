#include "FluctSplit/LaxWendroffSchemeCScalar.hh"
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

MethodStrategyProvider<LaxWendroffSchemeCScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
LaxWendroffSchemeScalarCProvider("ScalarLWC");

//////////////////////////////////////////////////////////////////////////////

LaxWendroffSchemeCScalar::LaxWendroffSchemeCScalar(const std::string& name) :
  RDS_SplitterScalar(name)
{
}

//////////////////////////////////////////////////////////////////////////////

LaxWendroffSchemeCScalar::~LaxWendroffSchemeCScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void LaxWendroffSchemeCScalar::setup()
{
  RDS_SplitterScalar::setup();

  m_sumKabs.resize(_nbEquations);
}

//////////////////////////////////////////////////////////////////////////////

void LaxWendroffSchemeCScalar::distribute(vector<RealVector>& residual)
{
  const CFreal nuCell = 1.0; // should this be like this ?
  const CFreal invnstates = 1.0 / _nbStatesInCell;
  
  // compute sumKabs for each scalar equation
  for (CFuint iEq = 0; iEq < _nbEquations; iEq++) {
      m_sumKabs[iEq] = std::abs(_k[0][iEq]);
      for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
         m_sumKabs[iEq] += std::abs(_k[iState][iEq]);
      }
      m_sumKabs[iEq] = 1.0 / std::max(MathTools::MathConsts::CFrealEps(),m_sumKabs[iEq]);
   }
  
  DistributionData& ddata = getMethodData().getDistributionData();
  // calculate the residual distribution
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    residual[iState] = _k[iState] * m_sumKabs;
    residual[iState] *= 0.5 * nuCell;
    residual[iState] += invnstates;
    residual[iState] *= ddata.phi;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LaxWendroffSchemeCScalar::distributePart(vector<RealVector>& residual)
{
  throw Common::NotImplementedException (FromHere(),"LaxWendroffSchemeCScalar::distributePart()");
}

//////////////////////////////////////////////////////////////////////////////

void LaxWendroffSchemeCScalar::computePicardJacob(vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"LaxWendroffSchemeCScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
