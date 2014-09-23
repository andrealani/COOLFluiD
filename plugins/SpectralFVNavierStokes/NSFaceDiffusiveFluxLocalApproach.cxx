#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFVNavierStokes/SpectralFVNavierStokes.hh"
#include "SpectralFVNavierStokes/NSFaceDiffusiveFluxLocalApproach.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<
    NSFaceDiffusiveFluxLocalApproach,SpectralFVMethodData,FaceDiffusiveFlux,SpectralFVNavierStokesModule >
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

/// @warning Not entirely happy with this formulation, think that states and gradients should be averaged first,
/// instead of the diffusive fluxes afterward... Cannot implement this because of the way DiffusiveVarSet is
/// implemented, getFlux() asks the states instead of the gradient variables...
vector< RealVector >& NSFaceDiffusiveFluxLocalApproach::computeDiffFlux(vector< vector< RealVector* >* >& lGrads,
                                                                        vector< vector< RealVector* >* >& rGrads,
                                                                        vector< RealVector* >& lStates,
                                                                        vector< RealVector* >& rStates,
                                                                        const CFreal& lCellVol,
                                                                        const CFreal& rCellVol,
                                                                        const CFreal& faceSurf,
                                                                        const vector< RealVector >& normal,
                                                                        const CFuint nbrFlxPnts)
{
  cf_assert(nbrFlxPnts == lStates.size());
  cf_assert(nbrFlxPnts == rStates.size());

  // inverse face length scale factor X alpha
  // not sure wether this is a good formulation...
  const CFreal dampFactorL = m_alpha/pow(faceSurf,m_dimM1Inv);
  const CFreal dampFactorR = dampFactorL;
//   const CFreal dampFactorL = m_alpha*faceSurf/lCellVol;
//   const CFreal dampFactorR = m_alpha*faceSurf/rCellVol;

  // compute the gradient variables in the flux points
  m_navierStokesVarSet->setGradientVars(lStates,m_flxPntLGradVars,nbrFlxPnts);
  m_navierStokesVarSet->setGradientVars(rStates,m_flxPntRGradVars,nbrFlxPnts);

  // compute the diffusive fluxes in the flux points
  cf_assert(m_multiDiffFlux.size() >= nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    // dereference the states
    RealVector& lFlxPntSol = *lStates[iFlx];
    RealVector& rFlxPntSol = *rStates[iFlx];

    // dereference the gradients
    vector< RealVector* >& lFlxPntGrads = *lGrads[iFlx];
    vector< RealVector* >& rFlxPntGrads = *rGrads[iFlx];

    // add damping term
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      // compute damping term
      const RealVector dGradVarXNormal
	= (m_flxPntLGradVars(iGrad,iFlx) - m_flxPntRGradVars(iGrad,iFlx))*normal[iFlx];
      m_dampTermL[iGrad] = dampFactorL*dGradVarXNormal;
      m_dampTermR[iGrad] = dampFactorR*dGradVarXNormal;

      // add damping term
      *lFlxPntGrads[iGrad] -= m_dampTermL[iGrad];
      *rFlxPntGrads[iGrad] -= m_dampTermR[iGrad];
    }

    // compute averaged diffusive flux
    m_multiDiffFlux[iFlx]  = m_oEminusBeta*m_diffusiveVarSet->getFlux(lFlxPntSol,lFlxPntGrads,normal[iFlx],0.);
    m_multiDiffFlux[iFlx] += m_beta        *m_diffusiveVarSet->getFlux(rFlxPntSol,rFlxPntGrads,normal[iFlx],0.);

    // remove damping term
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      *lFlxPntGrads[iGrad] += m_dampTermL[iGrad];
      *rFlxPntGrads[iGrad] += m_dampTermR[iGrad];
    }
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

  }  // namespace SpectralFV

}  // namespace COOLFluiD
