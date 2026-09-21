// Copyright (C) 2026 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/PhysicalModel.hh"
#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionTurb/GammaAlphaTransitionCriterion.hh"
#include "FluxReconstructionMethod/GradientVariables.hh"
#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/NSTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

GammaAlphaTransitionCriterion::GammaAlphaTransitionCriterion() :
  m_nbrEqs(),
  m_dim(),
  m_nbrSolPnts(),
  m_nbrSolSolDep(),
  m_solSolDep(CFNULL),
  m_solPolyDerivAtSolPnts(CFNULL),
  m_flxSolDep(CFNULL),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_eulerVarSet(CFNULL),
  m_diffusiveVarSet(CFNULL),
  m_navierStokesVarSet(CFNULL),
  m_pData(),
  m_gradVarStatePtrs(),
  m_gradVarsSolPnts(),
  m_flxPntGradVarsStore(),
  m_flxPntGradVars(),
  m_bndGradVarsStore(),
  m_bndGradVars(),
  m_cellFluxProjVects(),
  m_projectedCorr(),
  m_gradsSolPnts(),
  m_gradsFlxPntStore(),
  m_gradsFlxPnt()
{
}

//////////////////////////////////////////////////////////////////////////////

void GammaAlphaTransitionCriterion::setup(FluxReconstructionElementData& frData,
                                          const CFuint nbrEqs,
                                          const CFuint dim,
                                          SafePtr< ConvectiveVarSet > updateVarSet,
                                          SafePtr< DiffusiveVarSet > diffusiveVarSet)
{
  m_nbrEqs = nbrEqs;
  m_dim = dim;
  m_nbrSolPnts = frData.getNbrOfSolPnts();

  m_solSolDep = frData.getSolPntSolDependency();
  m_nbrSolSolDep = ((*m_solSolDep)[0]).size();
  m_solPolyDerivAtSolPnts = frData.getCoefSolPolyDerivInSolPnts();
  m_flxSolDep = frData.getFlxPntSolDependency();
  m_solPolyValsAtFlxPnts = frData.getCoefSolPolyInFlxPnts();

  // the largest number of flux points a face has
  const CFuint nbrFaceFlxPntsMax = frData.getFaceFlxPntsFaceLocalCoords()->size();

  m_eulerVarSet = updateVarSet.d_castTo< EulerVarSet >();
  m_diffusiveVarSet = diffusiveVarSet;
  m_navierStokesVarSet = diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  m_eulerVarSet->getModel()->resizePhysicalData(m_pData);

  m_gradVarStatePtrs.resize(m_nbrSolPnts);
  m_gradVarsSolPnts.resize(m_nbrEqs,m_nbrSolPnts);

  m_flxPntGradVarsStore.resize(nbrFaceFlxPntsMax);
  m_flxPntGradVars.resize(nbrFaceFlxPntsMax);
  m_bndGradVarsStore.resize(nbrFaceFlxPntsMax);
  m_bndGradVars.resize(nbrFaceFlxPntsMax);
  m_gradsFlxPntStore.resize(nbrFaceFlxPntsMax);
  m_gradsFlxPnt.resize(nbrFaceFlxPntsMax);
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPntsMax; ++iFlx)
  {
    m_flxPntGradVarsStore[iFlx].resize(m_nbrEqs);
    m_flxPntGradVars[iFlx] = &m_flxPntGradVarsStore[iFlx];
    m_bndGradVarsStore[iFlx].resize(m_nbrEqs);
    m_bndGradVars[iFlx] = &m_bndGradVarsStore[iFlx];
    m_gradsFlxPntStore[iFlx].resize(m_nbrEqs);
    m_gradsFlxPnt[iFlx].resize(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradsFlxPntStore[iFlx][iEq].resize(m_dim);
      m_gradsFlxPnt[iFlx][iEq] = &m_gradsFlxPntStore[iFlx][iEq];
    }
  }

  m_cellFluxProjVects.resize(m_dim);
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    m_cellFluxProjVects[iDim].resize(m_nbrSolPnts);
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_cellFluxProjVects[iDim][iSol].resize(m_dim);
    }
  }

  m_projectedCorr.resize(m_dim);

  m_gradsSolPnts.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_gradsSolPnts[iSol].resize(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradsSolPnts[iSol][iEq].resize(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void GammaAlphaTransitionCriterion::setTransitionFlags(BCStateComputer& bc,
                                                       const vector< State* >& cellStates,
                                                       const vector< State* >& cellStatesFlxPnt,
                                                       vector< State* >& ghostStates,
                                                       const vector< RealVector >& unitNormals,
                                                       const vector< RealVector >& flxPntCoords,
                                                       const vector< CFreal >& faceJacobVecSizeFlxPnts,
                                                       const vector< CFuint >& faceFlxPntConn,
                                                       const CFuint nbrFaceFlxPnts,
                                                       const vector< vector< CFreal > >& corrFctDiv,
                                                       const DataHandle< CFreal >& solPntNormals,
                                                       const DataHandle< CFreal >& volumes)
{
  // ghost states with every flag false, the boundary value rule reads them
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  {
    bc.setTransitionCriterion(iFlx,false);
  }
  bc.computeGhostStates(cellStatesFlxPnt,ghostStates,unitNormals,flxPntCoords);

  // gradient variables at the solution points and extrapolated to the flux points
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_gradVarStatePtrs[iSol] = cellStates[iSol]->getData();
  }
  m_diffusiveVarSet->setGradientVars(m_gradVarStatePtrs,m_gradVarsSolPnts,m_nbrSolPnts);
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdx = faceFlxPntConn[iFlx];
    extrapolateGradVarsToFlxPnt(m_gradVarsSolPnts,(*m_flxSolDep)[flxIdx],(*m_solPolyValsAtFlxPnts)[flxIdx],m_nbrEqs,*m_flxPntGradVars[iFlx]);
  }

  // the boundary value the boundary condition lifts them to
  bc.computeBndGradVars(m_flxPntGradVars,cellStatesFlxPnt,ghostStates,unitNormals,flxPntCoords,m_bndGradVars);

  // metric of the cell at the solution points
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const CFuint solID = cellStates[iSol]->getLocalID();
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      for (CFuint jDim = 0; jDim < m_dim; ++jDim)
      {
        m_cellFluxProjVects[iDim][iSol][jDim] = solPntNormals[solID*m_dim*m_dim+iDim*m_dim+jDim];
      }
    }
  }

  // gradient corrected with this face only: volume term plus the lifting of the jump g_b - g^D_f
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradsSolPnts[iSol][iEq] = 0.0;
    }
  }
  addGradVarsVolumeTerm(m_gradVarsSolPnts,m_cellFluxProjVects,*m_solSolDep,*m_solPolyDerivAtSolPnts,m_nbrSolSolDep,
                        m_nbrSolPnts,m_nbrEqs,m_dim,m_projectedCorr,m_gradsSolPnts);
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdx = faceFlxPntConn[iFlx];
    const CFuint nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFreal gradVarsJump = bndGradVarsJump((*m_flxPntGradVars[iFlx])[iEq],(*m_bndGradVars[iFlx])[iEq]);
      addGradVarsLifting(gradVarsJump,faceJacobVecSizeFlxPnts[iFlx],unitNormals[iFlx],
                         (*m_flxSolDep)[flxIdx],nbrSolDep,corrFctDiv,flxIdx,iEq,m_projectedCorr,m_gradsSolPnts);
    }
  }
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const CFreal invJacobDet = 1.0/volumes[cellStates[iSol]->getLocalID()];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradsSolPnts[iSol][iEq] *= invJacobDet;
    }
  }
  extrapolateGradsToFlxPnts(m_gradsSolPnts,faceFlxPntConn,nbrFaceFlxPnts,*m_flxSolDep,*m_solPolyValsAtFlxPnts,m_nbrEqs,m_gradsFlxPnt);

  // the flag at every flux point
  const CFreal pRef = m_eulerVarSet->getModel()->getPressRef();
  const CFreal TRef = m_eulerVarSet->getModel()->getTempRef();
  const CFreal muRef = (m_navierStokesVarSet->getModel().getReferencePhysicalData())[NSTerm::MU];
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  {
    const State& state = *cellStatesFlxPnt[iFlx];
    m_eulerVarSet->computePhysicalData(state,m_pData);
    const CFreal rho = m_navierStokesVarSet->getDensity(state);
    const CFreal mu = m_navierStokesVarSet->getModel().getDynViscosityDim(m_pData[EulerTerm::P]*pRef,m_pData[EulerTerm::T]*TRef)/muRef;
    const CFreal tau = wallShear(unitNormals[iFlx],m_gradsFlxPnt[iFlx],mu);
    const CFreal tauCrit = tau/sqrt(rho*mu);
    bc.setTransitionCriterion(iFlx,tauCrit <= state[5+m_dim]);
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal GammaAlphaTransitionCriterion::wallShear(const RealVector& normal,
                                                const vector< RealVector* >& grads,
                                                const CFreal mu) const
{
  // normal derivatives of the velocity components
  const CFreal dUdn = MathFunctions::innerProd(*grads[1],normal);
  const CFreal dVdn = MathFunctions::innerProd(*grads[2],normal);

  if (m_dim == 2)
  {
    // tangent (ny, -nx)
    return mu*(normal[YY]*dUdn - normal[XX]*dVdn);
  }

  const CFreal dWdn = MathFunctions::innerProd(*grads[3],normal);
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  CFreal tauT1 = 0.0;
  CFreal tauT2 = 0.0;
  if (fabs(nz) <= fabs(nx))
  {
    // tangents (ny, -nx, 0) and (nx nz, ny nz, -(nx^2 + ny^2)), normalised
    const CFreal nsize1 = sqrt(ny*ny + nx*nx);
    tauT1 = mu*(ny*dUdn - nx*dVdn)/nsize1;
    const CFreal nx2 = nx*nz;
    const CFreal ny2 = ny*nz;
    const CFreal nz2 = -(nx*nx + ny*ny);
    const CFreal nsize2 = sqrt(nx2*nx2 + ny2*ny2 + nz2*nz2);
    tauT2 = mu*(nx2*dUdn + ny2*dVdn + nz2*dWdn)/nsize2;
  }
  else
  {
    // tangents (0, -nz, ny) and (ny^2 + nz^2, -nx ny, -nx nz), normalised
    const CFreal nsize1 = sqrt(nz*nz + ny*ny);
    tauT1 = mu*(-nz*dVdn + ny*dWdn)/nsize1;
    const CFreal nx2 = ny*ny + nz*nz;
    const CFreal ny2 = -nx*ny;
    const CFreal nz2 = -nx*nz;
    const CFreal nsize2 = sqrt(nx2*nx2 + ny2*ny2 + nz2*nz2);
    tauT2 = mu*(nx2*dUdn + ny2*dVdn + nz2*dWdn)/nsize2;
  }
  return sqrt(tauT1*tauT1 + tauT2*tauT2);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
