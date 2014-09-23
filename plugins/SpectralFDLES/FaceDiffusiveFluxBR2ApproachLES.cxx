#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFDLES/SpectralFDLES.hh"
#include "SpectralFDLES/FaceDiffusiveFluxBR2ApproachLES.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<FaceDiffusiveFluxBR2ApproachLES,SpectralFDMethodData,FaceDiffusiveFlux,SpectralFDLESModule >
FaceDiffusiveFluxBR2ApproachLESProvider("LESBR2Approach");

//////////////////////////////////////////////////////////////////////////////

FaceDiffusiveFluxBR2ApproachLES::FaceDiffusiveFluxBR2ApproachLES(const std::string& name) :
  NSFaceDiffusiveFluxBR2Approach(name),
  m_lesVarSet(CFNULL),
  m_filterWidthVolumes(CFNULL)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

FaceDiffusiveFluxBR2ApproachLES::~FaceDiffusiveFluxBR2ApproachLES()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void FaceDiffusiveFluxBR2ApproachLES::setup()
{
  CFAUTOTRACE;

  // Call setup of the parent class
  NSFaceDiffusiveFluxBR2Approach::setup();

  // dynamic_cast the diffusive varset to a LESVarSet
  m_lesVarSet = m_diffusiveVarSet.d_castTo< LES::LESVarSet >();

}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& FaceDiffusiveFluxBR2ApproachLES::computeDiffFlux(vector< vector< RealVector* >* >& lGrads,
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
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    // dereference the states
    RealVector& lFlxPntSol = *lStates[iFlx];
    RealVector& rFlxPntSol = *rStates[iFlx];

    // dereference the gradients
    vector< RealVector* >& lFlxPntGrads = *lGrads[iFlx];
    vector< RealVector* >& rFlxPntGrads = *rGrads[iFlx];

    // compute averaged state
    *m_avgStates[iFlx] = 0.5*(lFlxPntSol+rFlxPntSol);

    // compute averaged gradients
    vector<RealVector* >& stateGrad = m_stateGradients[iFlx];
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      *stateGrad[iGrad] = 0.5*((*lFlxPntGrads[iGrad]) + (*rFlxPntGrads[iGrad]));
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

void FaceDiffusiveFluxBR2ApproachLES::unsetup()
{
  CFAUTOTRACE;

  // Call unsetup of the parent class
  NSFaceDiffusiveFluxBR2Approach::unsetup();

}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

