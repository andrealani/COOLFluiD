#include "MathTools/MatrixInverter.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FluctSplit/FluctSplitSystem.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/BSchemeCSys.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider< BSchemeCSys,
			FluctuationSplitData,
                        Splitter,
                        FluctSplitSystemModule>
aBSchemeCSysProvider("SysBC");

//////////////////////////////////////////////////////////////////////////////

void BSchemeCSys::defineConfigOptions(Config::OptionList& options)
{ 
  options.addConfigOption< CFuint, Config::DynamicOption<> >("FirstOrder","Run first order.");
}

//////////////////////////////////////////////////////////////////////////////

BSchemeCSys::BSchemeCSys(const std::string& name) :
  BSchemeBase<NSchemeCSys>(name),
  m_phiLDA()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
  
  m_firstOrder = 0;
  setParameter("FirstOrder",&m_firstOrder); 
}
      
//////////////////////////////////////////////////////////////////////////////

BSchemeCSys::~BSchemeCSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeCSys::configure ( Config::ConfigArgs& args )
{
  BSchemeBase<NSchemeCSys>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeCSys::setup()
{
  CFAUTOTRACE;
  
  BSchemeBase<NSchemeCSys>::setup();
  
  m_phiLDA.resize(_nbEquations);
}
      
//////////////////////////////////////////////////////////////////////////////

void BSchemeCSys::distribute(vector<RealVector>& residual)
{
  DistributionData& distdata = getMethodData().getDistributionData();
  const bool isPerturb = distdata.isPerturb;
  
  if ((m_firstOrder == 0 && !m_firstOrderJacob) || 
      (m_firstOrder == 0 && m_firstOrderJacob && !isPerturb)) {
    
    const RealVector& phiT = distdata.phi;
    const vector<State*>& tStates = *distdata.tStates;
    const CFuint nbEqs = _nbEquations;
    const CFuint nbStates = _nbStatesInCell;
    
    _sumKplusU = (*_kPlus[0])*(*tStates[0]);
    _sumKplus  = *_kPlus[0];
    for (CFuint iState = 1; iState < nbStates; ++iState) {
      _sumKplusU += (*_kPlus[iState])*(*tStates[iState]);
      _sumKplus  += *_kPlus[iState];
    }
    
    _inverter->invert(_sumKplus, _invK);
    _sumKplusU -= phiT;
    _uInflow = _invK * _sumKplusU;
    m_uTemp = _invK*phiT;
    
    // computation of the N residual and its sum
    m_sumPhiN = 0.0;
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      _uDiff = *tStates[iState] - _uInflow;
      m_phiN[iState] = (*_kPlus[iState])*_uDiff;
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  	m_sumPhiN[iEq] += std::abs(m_phiN[iState][iEq]);
      }
    }
    
    computeBlendingCoeff();
    
    if ( m_store_thetas ) storeThetas();
    
    if ( m_addExtraDiss ) addExtraDissipation(residual);
    
    
    // computation of LDA residual and the blending
    for (CFuint iState = 0; iState < nbStates; ++iState)
    {
	  m_phiLDA = (*_kPlus[iState])*m_uTemp;
  	  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  	residual[iState][iEq] =
	     m_theta[iEq]*m_phiN[iState][iEq] + (1. - m_theta[iEq])*m_phiLDA[iEq];
	  }

	  if (distdata.computeBetas)
      {
	    (*distdata.currBetaMat)[iState] = (*_kPlus[iState])*_invK;
	  }
    }
  }
  else {
    assert(m_firstOrder == 1 || (m_firstOrderJacob && isPerturb));
    BSchemeBase<NSchemeCSys>::distribute(residual);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void BSchemeCSys::distributePart(vector<RealVector>& residual)
{  
 DistributionData& distdata = getMethodData().getDistributionData();
  const bool isPerturb = distdata.isPerturb;
  
  if ((m_firstOrder == 0 && !m_firstOrderJacob) || 
      (m_firstOrder == 0 && m_firstOrderJacob && !isPerturb)) {
    
    const CFuint nbEqs = _nbEquations;
    RealVector& phi = distdata.phi;
    const vector<State*>& tStates = *distdata.tStates;
    
    _sumKplusU.slice(0, nbEqs) = (*_kPlus[0]) * tStates[0]->slice(_firstVarID, nbEqs);
    _sumKplus  =*_kPlus[0];
    for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
      _sumKplusU.slice(0, nbEqs) +=
	(*_kPlus[iState]) *  tStates[iState]->slice(_firstVarID, nbEqs);
      _sumKplus  += *_kPlus[iState];
    }
    
    _inverter->invert(_sumKplus, _invK);
    
    _sumKplusU.slice(0, nbEqs) -= phi.slice(_firstVarID, nbEqs);
    _uInflow = _invK*_sumKplusU;
    m_uTemp.slice(0, nbEqs) = _invK * phi.slice(_firstVarID, nbEqs);
    
    m_sumPhiN = 0.0;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
      {
	m_phiN[iState].slice(0, nbEqs) = (*_kPlus[iState])*
	  (tStates[iState]->slice(_firstVarID, nbEqs) -  _uInflow.slice(0, nbEqs));
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  m_sumPhiN[iEq] += std::abs(m_phiN[iState][iEq]);
	}
      }
    
    computeBlendingCoeff();
    
    if ( m_store_thetas ) storeThetas();
    
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
      {
	m_phiLDA = (*_kPlus[iState])*m_uTemp;
	for (CFuint iEq = _firstVarID, jEq = 0; iEq < _lastVarID; ++iEq, ++jEq)
	  {
	    residual[iState][iEq] =
	      m_theta[jEq]*m_phiN[iState][jEq] +  (1. - m_theta[jEq])*m_phiLDA[jEq];
	  }
      }
  }
  else {
    assert(m_firstOrder == 1 || (m_firstOrderJacob && isPerturb));
    BSchemeBase<NSchemeCSys>::distributePart(residual);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeCSys::computeBlendingCoeff()
{
  const RealVector& phi = getMethodData().getDistributionData().phi;
  const CFuint nbEqs = _nbEquations;
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
  {
    m_theta[iEq] = max ( std::abs(phi[iEq])/max(MathTools::MathConsts::CFrealEps(), m_sumPhiN[iEq]) , m_min_theta );
  }
}

//////////////////////////////////////////////////////////////////////////////

void BSchemeCSys::addExtraDissipation(vector<RealVector>& residual)
{
 // std::cout << "\tnothing is done\n";
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

   } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
