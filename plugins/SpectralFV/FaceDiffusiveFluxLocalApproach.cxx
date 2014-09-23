#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/FaceDiffusiveFluxLocalApproach.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<
    FaceDiffusiveFluxLocalApproach,SpectralFVMethodData,FaceDiffusiveFlux,SpectralFVModule >
  FaceDiffusiveFluxLocalApproachProvider("LocalApproach");

//////////////////////////////////////////////////////////////////////////////

void FaceDiffusiveFluxLocalApproach::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Alpha","Parameter of local approach determining the amount of damping.");
  options.addConfigOption< CFreal >("Beta","Parameter of local approach determining the wheight of left and right states.");
}

//////////////////////////////////////////////////////////////////////////////

FaceDiffusiveFluxLocalApproach::FaceDiffusiveFluxLocalApproach(const std::string& name) :
  FaceDiffusiveFlux(name),
  m_alpha(),
  m_beta(),
  m_oEminusBeta(),
  m_dimM1Inv(),/*
  m_avgSol(),
  m_avgGrads(),*/
  m_dampTermL(),
  m_dampTermR()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_alpha = 1.0;
  setParameter("Alpha",&m_alpha);

  m_beta = 0.5;
  setParameter("Beta" ,&m_beta );
}

//////////////////////////////////////////////////////////////////////////////

FaceDiffusiveFluxLocalApproach::~FaceDiffusiveFluxLocalApproach()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& FaceDiffusiveFluxLocalApproach::computeAvgGradVars(vector< RealVector* >& lStates,
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
      m_multiDiffFlux[iFlx][iGrad]  = m_beta        *(*lStates[iFlx])[iGrad];
      m_multiDiffFlux[iFlx][iGrad] += m_oEminusBeta*(*rStates[iFlx])[iGrad];
    }
  }

  return m_multiDiffFlux;
}

//////////////////////////////////////////////////////////////////////////////

/// @warning Not entirely happy with this formulation, think that states and gradients should be averaged first,
/// instead of the diffusive fluxes afterward... Cannot implement this because of the way DiffusiveVarSet is
/// implemented, getFlux() asks the states instead of the gradient variables...
vector< RealVector >& FaceDiffusiveFluxLocalApproach::computeDiffFlux(vector< vector< RealVector* >* >& lGrads,
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
          = ((*lStates[iFlx])[iGrad] - (*rStates[iFlx])[iGrad])*normal[iFlx];
      m_dampTermL[iGrad] = dampFactorL*dGradVarXNormal;
      m_dampTermR[iGrad] = dampFactorR*dGradVarXNormal;

      // add damping term
      *lFlxPntGrads[iGrad] -= m_dampTermL[iGrad];
      *rFlxPntGrads[iGrad] -= m_dampTermR[iGrad];
    }

    // compute averaged diffusive flux
    m_multiDiffFlux[iFlx]  = m_oEminusBeta*m_diffusiveVarSet->getFlux(lFlxPntSol,lFlxPntGrads,normal[iFlx],0);
    m_multiDiffFlux[iFlx] += m_beta        *m_diffusiveVarSet->getFlux(rFlxPntSol,rFlxPntGrads,normal[iFlx],0);

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

void FaceDiffusiveFluxLocalApproach::setup()
{
  CFAUTOTRACE;

  FaceDiffusiveFlux::setup();

  // compute wheight factor
  m_oEminusBeta = 1.0 - m_beta;

  // number of variables and dimensionality
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  m_dimM1Inv = 1.0/(dim - 1);

//   // resize m_avgSol
//   m_avgSol.resize(m_nbrEqs);
//
//   // create m_avgGrads
//   const CFuint m_nbrEqs = m_nbrEqs;
//   m_avgGrads.resize(m_nbrEqs);
//   for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
//   {
//     m_avgGrads[iGrad] = new RealVector(dim);
//   }

  // resize m_dampTermL and m_dampTermR
  m_dampTermL.resize(m_nbrEqs);
  m_dampTermR.resize(m_nbrEqs);
  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
  {
    m_dampTermL[iVar].resize(dim);
    m_dampTermR[iVar].resize(dim);
  }
}

//////////////////////////////////////////////////////////////////////////////

void FaceDiffusiveFluxLocalApproach::unsetup()
{
  CFAUTOTRACE;

//   for (CFuint iGrad = 0; iGrad < m_avgGrads.size(); ++iGrad)
//   {
//     deletePtr(m_avgGrads[iGrad]);
//   }
//   m_avgGrads.resize(0);

  FaceDiffusiveFlux::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

