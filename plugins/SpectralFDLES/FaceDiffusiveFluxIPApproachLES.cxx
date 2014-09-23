#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFDLES/SpectralFDLES.hh"
#include "SpectralFDLES/FaceDiffusiveFluxIPApproachLES.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<FaceDiffusiveFluxIPApproachLES,SpectralFDMethodData,FaceDiffusiveFlux,SpectralFDLESModule >
FaceDiffusiveFluxIPApproachLESProvider("LESIPApproach");

//////////////////////////////////////////////////////////////////////////////

FaceDiffusiveFluxIPApproachLES::FaceDiffusiveFluxIPApproachLES(const std::string& name) :
  NSFaceDiffusiveFluxIPApproach(name),
  m_lesVarSet(CFNULL),
  m_filterWidthVolumes(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

FaceDiffusiveFluxIPApproachLES::~FaceDiffusiveFluxIPApproachLES()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void FaceDiffusiveFluxIPApproachLES::setup()
{
  CFAUTOTRACE;

  // Call setup of the parent class
  NSFaceDiffusiveFluxIPApproach::setup();

  // dynamic_cast the diffusive varset to a LESVarSet
  m_lesVarSet = m_diffusiveVarSet.d_castTo< LES::LESVarSet >();

}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& FaceDiffusiveFluxIPApproachLES::computeDiffFlux(vector< vector< RealVector* >* >& lGrads,
                                                                     vector< vector< RealVector* >* >& rGrads,
                                                                     vector< RealVector* >& lStates,
                                                                     vector< RealVector* >& rStates,
                                                                     const std::vector< CFreal >& faceInvCharLength,
                                                                     const vector< RealVector >& normal,
                                                                     const CFuint nbrFlxPnts)
{
  cf_assert(nbrFlxPnts == lStates.size());
  cf_assert(nbrFlxPnts == rStates.size());

  // compute the gradient variable gradients
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

    // compute averaged state
    *m_avgStates[iFlx] = 0.5*(lFlxPntSol+rFlxPntSol);

    // compute averaged (damped) gradients
    vector<RealVector* >& stateGrad = m_stateGradients[iFlx];
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      // compute damping term
      const RealVector dGradVarXNormal
          = (lFlxPntSol[iGrad] - rFlxPntSol[iGrad])*normal[iFlx];
      *stateGrad[iGrad] = 0.5*((*lFlxPntGrads[iGrad]) + (*rFlxPntGrads[iGrad])) - dampFactor*dGradVarXNormal;
    }
  }

  // compute gradient variable gradients
  m_navierStokesVarSet->setGradientVarGradients(m_avgStates,
                                                m_stateGradients,
                                                m_gradVarGradients,
                                                nbrFlxPnts);

  // compute the diffusive fluxes in the flux points
  cf_assert(m_multiDiffFlux.size() >= nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    /// @todo setVolume is probably should be called in an other way.
    /// Here setVolume sets the filter width which will be used in the
    /// LES calculation to compute the diffusive flux in one the flux point
    /// of one face.
    // Filter Width Volume to use in the LES calculation
    m_lesVarSet->setVolume((*m_filterWidthVolumes)[iFlx]);

    m_multiDiffFlux[iFlx] = m_diffusiveVarSet->getFlux(*m_avgStates[iFlx],m_gradVarGradients[iFlx],
                                                       normal[iFlx],0);
  }

  return m_multiDiffFlux;
}


//////////////////////////////////////////////////////////////////////////////

void FaceDiffusiveFluxIPApproachLES::unsetup()
{
  CFAUTOTRACE;

  // Call unsetup of the parent class
  NSFaceDiffusiveFluxIPApproach::unsetup();

}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

