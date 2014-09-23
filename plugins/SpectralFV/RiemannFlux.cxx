#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/RiemannFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

RiemannFlux::RiemannFlux(const std::string& name) :
  SpectralFVMethodStrategy(name),
  m_maxNbrFlxPnts(),
  m_rFlux(),
  m_numDamping(),
  m_multiRFlux(),
  m_multiNumDamping(),
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
  m_numDamping.resize(m_nbrEqs);

  m_multiRFlux.resize(m_maxNbrFlxPnts);
  m_multiNumDamping.resize(m_maxNbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {
    m_multiRFlux[iFlx].resize(m_nbrEqs);
    m_multiNumDamping[iFlx].resize(m_nbrEqs);
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

