#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/NSFaceDiffusiveFluxLocalApproach.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<
    NSFaceDiffusiveFluxLocalApproach,SpectralFDMethodData,FaceDiffusiveFlux,SpectralFDNavierStokesModule >
NSFaceDiffusiveFluxLocalApproachProvider("NSLocalApproach");

//////////////////////////////////////////////////////////////////////////////

NSFaceDiffusiveFluxLocalApproach::NSFaceDiffusiveFluxLocalApproach(const std::string& name) :
  FaceDiffusiveFluxLocalApproach(name),
  m_navierStokesVarSet(CFNULL),
  m_flxPntLGradVars(),
  m_flxPntRGradVars()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NSFaceDiffusiveFluxLocalApproach::~NSFaceDiffusiveFluxLocalApproach()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& NSFaceDiffusiveFluxLocalApproach::computeAvgGradVars(vector< RealVector* >& lStates,
                                                                           vector< RealVector* >& rStates,
                                                                           const CFuint nbrFlxPnts)
{
  cf_assert(nbrFlxPnts <= lStates.size());
  cf_assert(nbrFlxPnts <= rStates.size());

  // compute the gradient variables in the flux points
  m_navierStokesVarSet->setGradientVars(lStates,m_flxPntLGradVars,nbrFlxPnts);
  m_navierStokesVarSet->setGradientVars(rStates,m_flxPntRGradVars,nbrFlxPnts);

  // compute the averaged gradient variables
  cf_assert(m_multiDiffFlux.size() >= nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      m_multiDiffFlux[iFlx][iGrad]  = m_beta*m_flxPntLGradVars(iGrad,iFlx);
      m_multiDiffFlux[iFlx][iGrad] += m_oEminusBeta*m_flxPntRGradVars(iGrad,iFlx);
    }
  }

  return m_multiDiffFlux;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& NSFaceDiffusiveFluxLocalApproach::computeDiffFlux(vector< vector< RealVector* >* >& lGrads,
                                                                        vector< vector< RealVector* >* >& rGrads,
                                                                        vector< RealVector* >& lStates,
                                                                        vector< RealVector* >& rStates,
                                                                        const std::vector< CFreal >& faceInvCharLength,
                                                                        const vector< RealVector >& normal,
                                                                        const CFuint nbrFlxPnts)
{
  cf_assert(nbrFlxPnts == lStates.size());
  cf_assert(nbrFlxPnts == rStates.size());

  // compute the gradient variables in the flux points
  m_navierStokesVarSet->setGradientVars(lStates,m_flxPntLGradVars,nbrFlxPnts);
  m_navierStokesVarSet->setGradientVars(rStates,m_flxPntRGradVars,nbrFlxPnts);

  // compute the diffusive fluxes in the flux points
  cf_assert(m_multiDiffFlux.size() >= nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    // damping factor
    const CFreal dampFactor = m_alpha*faceInvCharLength[iFlx];

    // dereference the states
    RealVector& lFlxPntSol = *lStates[iFlx];
    RealVector& rFlxPntSol = *rStates[iFlx];

    // dereference the gradients
    vector< RealVector* >& lFlxPntGrads = *lGrads[iFlx];
    vector< RealVector* >& rFlxPntGrads = *rGrads[iFlx];

    // compute averaged solution
    m_avgSol = 0.5*(lFlxPntSol+rFlxPntSol);

    // compute averaged (damped) gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      // compute damping term
      const RealVector dGradVarXNormal = (m_flxPntLGradVars(iGrad,iFlx) - m_flxPntRGradVars(iGrad,iFlx))*normal[iFlx];
      *m_avgGrads[iGrad] = m_oEminusBeta*(*lFlxPntGrads[iGrad]) + m_beta*(*rFlxPntGrads[iGrad]) - dampFactor*dGradVarXNormal;
    }

    // compute averaged diffusive flux
    m_multiDiffFlux[iFlx] = m_diffusiveVarSet->getFlux(m_avgSol,m_avgGrads,normal[iFlx],0);
  }

  return m_multiDiffFlux;
}

//////////////////////////////////////////////////////////////////////////////

void NSFaceDiffusiveFluxLocalApproach::setup()
{
  CFAUTOTRACE;

  FaceDiffusiveFluxLocalApproach::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  // create gradient variables
  m_flxPntLGradVars.resize(m_nbrEqs,m_maxNbrFlxPnts);
  m_flxPntRGradVars.resize(m_nbrEqs,m_maxNbrFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void NSFaceDiffusiveFluxLocalApproach::unsetup()
{
  CFAUTOTRACE;

  FaceDiffusiveFluxLocalApproach::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
