#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/FaceDiffusiveFluxLocalApproach.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<
    FaceDiffusiveFluxLocalApproach,SpectralFDMethodData,FaceDiffusiveFlux,SpectralFDModule >
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
//   m_dimM1Inv(),
  m_avgSol(),
  m_avgGrads()/*,
  m_dampTermL(),
  m_dampTermR()*/
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

vector< RealVector >& FaceDiffusiveFluxLocalApproach::computeDiffFlux(vector< vector< RealVector* >* >& lGrads,
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
      const RealVector dGradVarXNormal
          = ((*lStates[iFlx])[iGrad] - (*rStates[iFlx])[iGrad])*normal[iFlx];
      *m_avgGrads[iGrad] = m_oEminusBeta*(*lFlxPntGrads[iGrad]) + m_beta*(*rFlxPntGrads[iGrad]) - dampFactor*dGradVarXNormal;
    }

    // compute averaged diffusive flux
    m_multiDiffFlux[iFlx] = m_diffusiveVarSet->getFlux(m_avgSol,m_avgGrads,normal[iFlx],0);
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

void FaceDiffusiveFluxLocalApproach::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iGrad = 0; iGrad < m_avgGrads.size(); ++iGrad)
  {
    deletePtr(m_avgGrads[iGrad]);
  }
  m_avgGrads.resize(0);

  FaceDiffusiveFlux::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

