#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFDLES/SpectralFDLES.hh"
#include "SpectralFDLES/FaceDiffusiveFluxLocalApproachLES.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<FaceDiffusiveFluxLocalApproachLES,SpectralFDMethodData,FaceDiffusiveFlux,SpectralFDLESModule >
FaceDiffusiveFluxLocalApproachLESProvider("LESLocalApproach");

//////////////////////////////////////////////////////////////////////////////

FaceDiffusiveFluxLocalApproachLES::FaceDiffusiveFluxLocalApproachLES(const std::string& name) :
  NSFaceDiffusiveFluxLocalApproach(name),
  m_lesVarSet(CFNULL),
  m_filterWidthVolumes(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

FaceDiffusiveFluxLocalApproachLES::~FaceDiffusiveFluxLocalApproachLES()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void FaceDiffusiveFluxLocalApproachLES::setup()
{
  CFAUTOTRACE;

  // Call setup of the parent class
  NSFaceDiffusiveFluxLocalApproach::setup();

  // Dynamic_cast the diffusive varset to a LESVarSet
  m_lesVarSet = m_diffusiveVarSet.d_castTo< LES::LESVarSet >();

}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& FaceDiffusiveFluxLocalApproachLES::computeDiffFlux(vector< vector< RealVector* >* >& lGrads,
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
      *m_avgGrads[iGrad] = m_o.minusBeta*(*lFlxPntGrads[iGrad]) + m_beta*(*rFlxPntGrads[iGrad]) - dampFactor*dGradVarXNormal;
    }

    /// @todo setVolume is probably should be called in an other way.
    /// Here setVolume sets the filter width which will be used in the
    /// LES calculation to compute the diffusive flux in one the flux point
    /// of one face.
    // Filter Width Volume to use in the LES calculation
    m_lesVarSet->setVolume((*m_filterWidthVolumes)[iFlx]);

    // compute averaged diffusive flux
    m_multiDiffFlux[iFlx] = m_diffusiveVarSet->getFlux(m_avgSol,m_avgGrads,normal[iFlx],0);
  }

  return m_multiDiffFlux;
}

//////////////////////////////////////////////////////////////////////////////

void FaceDiffusiveFluxLocalApproachLES::unsetup()
{
  CFAUTOTRACE;

  // Call unsetup of the parent class
  NSFaceDiffusiveFluxLocalApproach::unsetup();

}
//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
