#include "Framework/MethodStrategyProvider.hh"
#include "Framework/BaseTerm.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

RiemannFlux::RiemannFlux(const std::string& name) :
  FluxReconstructionSolverStrategy(name),
  m_maxNbrFlxPnts(),
  m_rFlux(),
  m_multiRFlux(),
  m_nbrEqs()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

RiemannFlux::~RiemannFlux()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void RiemannFlux::setup()
{
  CFAUTOTRACE;

  // get number of equations
  m_nbrEqs = Framework::PhysicalModelStack::getActive()->getNbEq();

  // get maximum number of flux points
  m_maxNbrFlxPnts = getMethodData().getMaxNbrRFluxPnts();

  m_rFlux.resize(m_nbrEqs);

  m_multiRFlux.resize(m_maxNbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {
    m_multiRFlux[iFlx].resize(m_nbrEqs);
  }
  
  // resize m_pData which will contains the physical data for the left and the right internal flux points
  m_pData.resize(2);
  Common::SafePtr<Framework::BaseTerm> convTerm = Framework::PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  for (CFuint i = 0; i < m_pData.size(); ++i) {
    convTerm->resizePhysicalData(m_pData[i]);
  }
  
  m_solStates.resize(2);
  m_updateStates.resize(2);
}

//////////////////////////////////////////////////////////////////////////////

void RiemannFlux::unsetup()
{
  CFAUTOTRACE;
}
    
//////////////////////////////////////////////////////////////////////////////

std::vector< RealVector >& RiemannFlux::computeFlux(std::vector< Framework::State* >& lState,
                                                    std::vector< Framework::State* >& rState,
                                                    const std::vector< RealVector >& normal,
                                                    const CFuint nbrFlxPnts)
{
  cf_assert(nbrFlxPnts <= lState.size());
  cf_assert(nbrFlxPnts <= rState.size());
  cf_assert(nbrFlxPnts <= normal.size());
  
  // compute the flux for every flux point
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx) {
    m_multiRFlux[iFlx] = computeFlux((*lState[iFlx]),(*rState[iFlx]),normal[iFlx]);
  }
  
  return m_multiRFlux;
}

//////////////////////////////////////////////////////////////////////////////

RealVector& RiemannFlux::computeFlux(Framework::State& lState,
                                     Framework::State& rState,
                                     const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"RiemannFlux::computeFlux must be overriden");
}
    
//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

