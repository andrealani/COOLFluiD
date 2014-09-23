#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/FaceDiffusiveFluxBR2Approach.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<
    FaceDiffusiveFluxBR2Approach,SpectralFDMethodData,FaceDiffusiveFlux,SpectralFDModule >
FaceDiffusiveFluxBR2ApproachProvider("BR2Approach");

//////////////////////////////////////////////////////////////////////////////

FaceDiffusiveFluxBR2Approach::FaceDiffusiveFluxBR2Approach(const std::string& name) :
  FaceDiffusiveFluxCompactApproach(name),
  m_solPolyOrderFac(),
  m_avgSol(),
  m_avgGrads()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

FaceDiffusiveFluxBR2Approach::~FaceDiffusiveFluxBR2Approach()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& FaceDiffusiveFluxBR2Approach::computeAvgGradVars(vector< RealVector* >& lStates,
                                                                       vector< RealVector* >& rStates,
                                                                       const CFuint nbrFlxPnts)
{
  cf_assert(nbrFlxPnts <= lStates.size());
  cf_assert(nbrFlxPnts <= rStates.size());

  // compute the averaged gradient variables
  cf_assert(m_multiDiffFlux.size() >= nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      m_multiDiffFlux[iFlx][iGrad]  = 0.5*((*lStates[iFlx])[iGrad] +
                                           (*rStates[iFlx])[iGrad]);
    }
  }

  return m_multiDiffFlux;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& FaceDiffusiveFluxBR2Approach::computeAvgGradVar(const CFuint iVar,
                                                                      vector< RealVector* >& lStates,
                                                                      vector< RealVector* >& rStates,
                                                                      const CFuint nbrFlxPnts)
{
  cf_assert(nbrFlxPnts <= lStates.size());
  cf_assert(nbrFlxPnts <= rStates.size());

  // compute the averaged gradient variables
  cf_assert(m_multiDiffFlux.size() >= nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_multiDiffFlux[iFlx][iVar]  = 0.5*((*lStates[iFlx])[iVar] +
                                        (*rStates[iFlx])[iVar]);
  }

  return m_multiDiffFlux;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& FaceDiffusiveFluxBR2Approach::computeDiffFlux(vector< vector< RealVector* >* >& lGrads,
                                                                    vector< vector< RealVector* >* >& rGrads,
                                                                    vector< RealVector* >& lStates,
                                                                    vector< RealVector* >& rStates,
                                                                    const std::vector< CFreal >& faceInvCharLength,
                                                                    const vector< RealVector >& normal,
                                                                    const CFuint nbrFlxPnts)
{
  cf_assert(nbrFlxPnts == lStates.size());
  cf_assert(nbrFlxPnts == rStates.size());

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

    // compute averaged solution
    m_avgSol = 0.5*(lFlxPntSol+rFlxPntSol);

    // compute averaged gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      *m_avgGrads[iGrad] = 0.5*((*lFlxPntGrads[iGrad]) + (*rFlxPntGrads[iGrad]));
//       CF_DEBUG_OBJ(*m_avgGrads[iGrad]);
    }

    // compute averaged diffusive flux
    m_multiDiffFlux[iFlx] = m_diffusiveVarSet->getFlux(m_avgSol,m_avgGrads,normal[iFlx],0);
//     CF_DEBUG_OBJ(m_multiDiffFlux[iFlx]);
  }

  return m_multiDiffFlux;
}

//////////////////////////////////////////////////////////////////////////////

void FaceDiffusiveFluxBR2Approach::setup()
{
  CFAUTOTRACE;

  FaceDiffusiveFluxCompactApproach::setup();

  // number of variables and dimensionality
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  // resize m_avgSol
  m_avgSol.resize(m_nbrEqs);

  // create m_avgGrads
  m_avgGrads.resize(m_nbrEqs);
  for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
  {
    m_avgGrads[iGrad] = new RealVector(dim);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FaceDiffusiveFluxBR2Approach::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iGrad = 0; iGrad < m_avgGrads.size(); ++iGrad)
  {
    deletePtr(m_avgGrads[iGrad]);
  }
  m_avgGrads.resize(0);

  FaceDiffusiveFluxCompactApproach::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

