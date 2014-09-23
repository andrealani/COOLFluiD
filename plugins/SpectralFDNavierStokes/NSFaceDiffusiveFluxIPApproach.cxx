#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"
#include "SpectralFDNavierStokes/NSFaceDiffusiveFluxIPApproach.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<
    NSFaceDiffusiveFluxIPApproach,SpectralFDMethodData,FaceDiffusiveFlux,SpectralFDNavierStokesModule >
NSFaceDiffusiveFluxIPApproachProvider("NSIPApproach");

//////////////////////////////////////////////////////////////////////////////

NSFaceDiffusiveFluxIPApproach::NSFaceDiffusiveFluxIPApproach(const std::string& name) :
  FaceDiffusiveFluxIPApproach(name),
  m_navierStokesVarSet(CFNULL),
//   m_bcStateComputer(CFNULL),
  m_flxPntLGradVars(),
  m_flxPntRGradVars(),
  m_avgStates(),
  m_stateGradients(),
  m_gradVarGradients()/*,
  m_ghostGradVarGradients()*/
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

NSFaceDiffusiveFluxIPApproach::~NSFaceDiffusiveFluxIPApproach()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector >& NSFaceDiffusiveFluxIPApproach::computeDiffFlux(vector< vector< RealVector* >* >& lGrads,
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
    m_multiDiffFlux[iFlx] = m_diffusiveVarSet->getFlux(*m_avgStates[iFlx],m_gradVarGradients[iFlx],
                                                       normal[iFlx],0);
  }

  return m_multiDiffFlux;
}

//////////////////////////////////////////////////////////////////////////////

// vector< RealVector >& NSFaceDiffusiveFluxIPApproach::computeBndDiffFlux(const vector< vector< RealVector* > >& intGrads,
//                                                                         const vector< RealVector* >& intStates,
//                                                                         const vector< RealVector* >& ghostStates,
//                                                                         const std::vector< CFreal >& faceInvCharLength,
//                                                                         const vector< RealVector >& normal,
//                                                                         const vector< RealVector >& flxPntCoords,
//                                                                         const CFuint nbrFlxPnts)
// {
//   cf_assert(nbrFlxPnts == intStates  .size());
//   cf_assert(nbrFlxPnts == ghostStates.size());
//
//   // compute internal gradient variable gradients
//   m_navierStokesVarSet->setGradientVarGradients(intStates,
//                                                 intGrads,
//                                                 m_gradVarGradients,
//                                                 nbrFlxPnts);
//
//   // compute ghost gradient variable gradients
//   m_bcStateComputer->computeGhostGradients(m_gradVarGradients,
//                                            m_ghostGradVarGradients,
//                                            normal,
//                                            flxPntCoords);
//
//
//   // compute the damping terms and average internal and ghost gradient variable gradients
//   for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
//   {
//     const CFreal dampFactor = m_alpha*faceInvCharLength[iFlx];
//
//     // dereference the states
//     RealVector& intFlxPntSol   = *intStates  [iFlx];
//     RealVector& ghostFlxPntSol = *ghostStates[iFlx];
//
//     // compute averaged state
//     *m_avgStates[iFlx] = 0.5*(intFlxPntSol+ghostFlxPntSol);
//
//     // compute damping term and average gradients
//     vector<RealVector* >& stateGrad = m_stateGradients[iFlx];
//     for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
//     {
//       // damping term
//       const RealVector dGradVarXNormal
//           = (intFlxPntSol[iGrad] - ghostFlxPntSol[iGrad])*normal[iFlx];
//       *stateGrad[iGrad] = - dampFactor*dGradVarXNormal;
//
//       // averaging
//       *m_ghostGradVarGradients[iFlx][iGrad] += *m_gradVarGradients[iFlx][iGrad];
//       *m_ghostGradVarGradients[iFlx][iGrad] *= 0.5;
//     }
//   }
//
//   // convert damping terms
//   m_navierStokesVarSet->setGradientVarGradients(m_avgStates,
//                                                 m_stateGradients,
//                                                 m_gradVarGradients,
//                                                 nbrFlxPnts);
//
//   // compute the diffusive fluxes in the flux points
//   cf_assert(m_multiDiffFlux.size() >= nbrFlxPnts);
//   for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
//   {
//     // add damping terms to averaged gradient variable gradients
//     vector<RealVector* >& gradVarGrad    = m_gradVarGradients     [iFlx];
//     vector<RealVector* >& avgGradVarGrad = m_ghostGradVarGradients[iFlx];
//     for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
//     {
//       *gradVarGrad[iGrad] += *avgGradVarGrad[iGrad];
//     }
//
//     // compute averaged diffusive flux
//     m_multiDiffFlux[iFlx] = m_diffusiveVarSet->getFlux(*m_avgStates[iFlx],m_gradVarGradients[iFlx],
//                                                        normal[iFlx],0);
//   }
//
//   return m_multiDiffFlux;
// }

//////////////////////////////////////////////////////////////////////////////

void NSFaceDiffusiveFluxIPApproach::setup()
{
  CFAUTOTRACE;

  FaceDiffusiveFluxIPApproach::setup();

  // dynamic_cast the diffusive varset to a navier-stokes varset
  m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  // create gradient variables
  m_flxPntLGradVars.resize(m_nbrEqs,m_maxNbrFlxPnts);
  m_flxPntRGradVars.resize(m_nbrEqs,m_maxNbrFlxPnts);

  // create averaged states
  m_avgStates.resize(0);
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {
    m_avgStates.push_back(new RealVector(m_nbrEqs));
  }

  // dimensionality
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  // create gradient variable gradients
  m_stateGradients       .resize(0);
  m_gradVarGradients     .resize(0);
//   m_ghostGradVarGradients.resize(0);
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {
    vector<RealVector*> flxPntStateGrads       ;
    vector<RealVector*> flxPntGradVarGrads     ;
//     vector<RealVector*> flxPntGhostGradVarGrads;
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      flxPntStateGrads       .push_back(new RealVector(dim));
      flxPntGradVarGrads     .push_back(new RealVector(dim));
//       flxPntGhostGradVarGrads.push_back(new RealVector(dim));
    }
    m_stateGradients       .push_back(flxPntStateGrads       );
    m_gradVarGradients     .push_back(flxPntGradVarGrads     );
//     m_ghostGradVarGradients.push_back(flxPntGhostGradVarGrads);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NSFaceDiffusiveFluxIPApproach::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iFlx = 0; iFlx < m_avgStates.size(); ++iFlx)
  {
    deletePtr(m_avgStates[iFlx]);
  }
  m_avgStates.resize(0);

  for (CFuint iFlx = 0; iFlx < m_gradVarGradients.size(); ++iFlx)
  {
    for (CFuint iEq = 0; iEq < m_gradVarGradients[iFlx].size(); ++iEq)
    {
      deletePtr(m_stateGradients       [iFlx][iEq]);
      deletePtr(m_gradVarGradients     [iFlx][iEq]);
//       deletePtr(m_ghostGradVarGradients[iFlx][iEq]);
    }
    m_stateGradients       [iFlx].resize(0);
    m_gradVarGradients     [iFlx].resize(0);
//     m_ghostGradVarGradients[iFlx].resize(0);
  }
  m_stateGradients       .resize(0);
  m_gradVarGradients     .resize(0);
//   m_ghostGradVarGradients.resize(0);

  FaceDiffusiveFluxIPApproach::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
