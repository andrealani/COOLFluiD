#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/FaceDiffusiveFluxIPApproach.hh"
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
    FaceDiffusiveFluxIPApproach,SpectralFDMethodData,FaceDiffusiveFlux,SpectralFDModule >
FaceDiffusiveFluxIPApproachProvider("IPApproach");

//////////////////////////////////////////////////////////////////////////////

void FaceDiffusiveFluxIPApproach::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Alpha","Parameter of Interior Penalty approach determining the amount of damping.");
}

//////////////////////////////////////////////////////////////////////////////

FaceDiffusiveFluxIPApproach::FaceDiffusiveFluxIPApproach(const std::string& name) :
  FaceDiffusiveFluxCompactApproach(name),
  m_alpha(),
  m_solPolyOrderFac(),
  m_avgSol(),
  m_avgGrads()/*,
  m_dampTermL(),
  m_dampTermR()*/
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_alpha = 10.0;
  setParameter("Alpha",&m_alpha);
}

//////////////////////////////////////////////////////////////////////////////

FaceDiffusiveFluxIPApproach::~FaceDiffusiveFluxIPApproach()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& FaceDiffusiveFluxIPApproach::computeAvgGradVars(vector< RealVector* >& lStates,
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

vector< RealVector >& FaceDiffusiveFluxIPApproach::computeAvgGradVar(const CFuint iVar,
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

vector< RealVector >& FaceDiffusiveFluxIPApproach::computeDiffFlux(vector< vector< RealVector* >* >& lGrads,
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
      const RealVector dGradVarXNormal
          = (lFlxPntSol[iGrad] - rFlxPntSol[iGrad])*normal[iFlx];
      *m_avgGrads[iGrad] = 0.5*((*lFlxPntGrads[iGrad]) + (*rFlxPntGrads[iGrad])) - dampFactor*dGradVarXNormal;
    }

    // compute averaged diffusive flux
    m_multiDiffFlux[iFlx] = m_diffusiveVarSet->getFlux(m_avgSol,m_avgGrads,normal[iFlx],0);
  }

  return m_multiDiffFlux;
}

//////////////////////////////////////////////////////////////////////////////

void FaceDiffusiveFluxIPApproach::setup()
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

//   // resize m_dampTermL and m_dampTermR
//   m_dampTermL.resize(m_nbrEqs);
//   m_dampTermR.resize(m_nbrEqs);
//   for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
//   {
//     m_dampTermL[iVar].resize(dim);
//     m_dampTermR[iVar].resize(dim);
//   }

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // number of solution points
  const CFuint polyOrder = sdLocalData[0]->getPolyOrder();

  // set m_alpha
  if (polyOrder == 0)
  {
    m_alpha = 1.0;
  }
  else
  {
//     m_alpha = m_alpha*(polyOrder + 1)*(polyOrder + 1);
    m_alpha = m_alpha*(polyOrder*polyOrder + polyOrder + 1);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FaceDiffusiveFluxIPApproach::unsetup()
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

