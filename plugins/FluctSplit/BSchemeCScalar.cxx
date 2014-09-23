#include "MathTools/MatrixInverter.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitScalar.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/BSchemeCScalar.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<BSchemeCScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
bcSchemeScalarProvider("ScalarBC");

//////////////////////////////////////////////////////////////////////////////

void BSchemeCScalar::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

BSchemeCScalar::BSchemeCScalar(const std::string& name) :
  BSchemeBase<RDS_SplitterScalar>(name),
  m_sumKplusU(),
  m_sumKplus(),
  m_phiLDA(),
  m_temp(),
  m_uMin()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

BSchemeCScalar::~BSchemeCScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeCScalar::configure ( Config::ConfigArgs& args )
{
  BSchemeBase<RDS_SplitterScalar>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeCScalar::setup()
{
  BSchemeBase<RDS_SplitterScalar>::setup();

  m_sumKplusU.resize(_nbEquations);
  m_sumKplus.resize(_nbEquations);
  m_phiLDA.resize(_nbEquations);
  m_theta.resize(_nbEquations);
  m_uTemp.resize(_nbEquations);
  m_uInFlow.resize(_nbEquations);
  m_temp.resize(_nbEquations);
  m_uMin.resize(_nbEquations);
}


//////////////////////////////////////////////////////////////////////////////

void BSchemeCScalar::distribute(vector<RealVector>& residual)
{
  DistributionData& distdata = getMethodData().getDistributionData();
  const RealVector& phi = distdata.phi;
  const vector<State*>& tStates = *distdata.tStates;
  const CFuint nbEqs = _nbEquations;
  
  // computation of the sum of Kplus
  m_sumKplusU = _kPlus[0] * (*tStates[0]);
  m_sumKplus  = _kPlus[0];
  const CFuint nbStatesInCell = _nbStatesInCell;
  for (CFuint iState = 1; iState < nbStatesInCell; ++iState) {
    m_sumKplusU += _kPlus[iState] * (*tStates[iState]);
    m_sumKplus  += _kPlus[iState];
  }

  // computation of inflow parameters
  m_temp = ( m_sumKplusU - phi );
  m_uInFlow = m_temp / m_sumKplus;

  // computation of the N residual and its sum
  m_sumPhiN = 0.0;
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    m_uMin = (*tStates[iState] - m_uInFlow);
    m_phiN[iState] = _kPlus[iState]*m_uMin;
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      m_sumPhiN[iEq] += std::abs(m_phiN[iState][iEq]);
    }
  }

  computeBlendingCoeff();
  
  if ( m_store_thetas ) storeThetas();
  
  // computation of LDA residual and the blending
  m_uTemp = phi / m_sumKplus;
  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    m_phiLDA = _kPlus[iState]*m_uTemp;
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      residual[iState][iEq] =  m_theta[iEq]*m_phiN[iState][iEq] +
	(1. - m_theta[iEq]) * m_phiLDA[iEq];
      
      if (getMethodData().getDistributionData().computeBetas) {
	(*getMethodData().getDistributionData().currBetaMat)[iState] = _kPlus[iState]/m_sumKplus;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeCScalar::distributePart(vector<RealVector>& residual)
{
  throw Common::NotImplementedException (FromHere(),"BSchemeCScalar::distributePart()");
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeCScalar::computePicardJacob(vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"BSchemeCScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeCScalar::computeBlendingCoeff()
{ 
  const RealVector& phi = getMethodData().getDistributionData().phi;
  const CFuint nbEqs = _nbEquations;
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
  {
    m_theta[iEq] = max ( std::abs(phi[iEq])/max(MathTools::MathConsts::CFrealEps(), m_sumPhiN[iEq]) , m_min_theta );
  }
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeCScalar::addExtraDissipation()
{
 throw Common::NotImplementedException (FromHere(),"BSchemeCScalar::addExtraDissipation()");
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
