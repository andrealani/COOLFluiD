// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <algorithm>
#include <cmath>

#include "Common/NotImplementedException.hh"
#include "Common/BadValueException.hh"
#include "Common/StringOps.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/GradientVariables.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/DiffRHSJacobFluxReconstructionNSSubcellBlending.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffRHSJacobFluxReconstructionNSSubcellBlending,
                       FluxReconstructionSolverData,
                       FluxReconstructionNavierStokesModule >
diffRHSJacobNSSubcellBlendingProvider("DiffRHSJacobNSSubcellBlending");

//////////////////////////////////////////////////////////////////////////////

namespace {

/// exp(x) with the argument capped below the overflow of double precision
inline CFreal cappedExp(const CFreal x)
{
  return std::exp(std::min(x, 700.0));
}

}

//////////////////////////////////////////////////////////////////////////////

DiffRHSJacobFluxReconstructionNSSubcellBlending::DiffRHSJacobFluxReconstructionNSSubcellBlending(const std::string& name) :
  DiffRHSJacobFluxReconstructionNS(name),
  socket_alpha("alpha"),
  m_scData(),
  m_nbrSolPnts1D(0),
  m_solPnts1D(),
  m_widths1D(),
  m_closestSolToFlx(CFNULL),
  m_nbrCellFlxPnts(0),
  m_lineEndFlx(),
  m_faceFactors(),
  m_intfSolL(),
  m_intfSolR(),
  m_intfWidthL(),
  m_intfWidthR(),
  m_intfOfSol(),
  m_isLogGradVar(),
  m_cellOfState(),
  m_faceSamples(),
  m_sampleStamp(),
  m_subGrads(),
  m_subGradStamp(),
  m_cellData(),
  m_stamp(0),
  m_blendGradVars(),
  m_sideGradVars(),
  m_blendProjVects(),
  m_twoStatePtrs(),
  m_twoGradVars(),
  m_twoGrad(),
  m_twoDir(),
  m_twoFlux(),
  m_twoFluxPert()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
DiffRHSJacobFluxReconstructionNSSubcellBlending::needsSockets()
{
  std::vector< SafePtr< BaseDataSocketSink > > result = DiffRHSJacobFluxReconstructionNS::needsSockets();
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::execute()
{
  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();

  bool anyAlpha = false;
  for (CFuint i = 0; i < alpha.size(); ++i)
  {
    if (alpha[i] > 0.0) { anyAlpha = true; break; }
  }

  // subcell data of the alpha > 0 cells from the unperturbed states, used by the
  // virtual functions the parent calls
  if (anyAlpha)
  {
    ++m_stamp;
    computeFaceSamples();
    prepareAlphaCells();
  }

  DiffRHSJacobFluxReconstructionNS::execute();
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::computeFaceSamples()
{
  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size()-1;

  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder->getDataGE();
  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;

  for (CFuint orient = 0; orient < nbrFaceOrients; ++orient)
  {
    // cell flux point indices of the face flux points, [side][face flux point]
    const vector< vector< CFuint > >& flxConn = (*m_faceFlxPntConnPerOrient)[orient];
    const CFuint nbrFaceFlx = flxConn[LEFT].size();

    for (CFuint faceID = innerFacesStartIdxs[orient]; faceID < innerFacesStartIdxs[orient+1]; ++faceID)
    {
      geoDataFace.idx = faceID;
      GeometricEntity* face = m_faceBuilder->buildGE();

      vector< State* >* states[2];
      CFuint cellIDs[2];
      for (CFuint iSide = 0; iSide < 2; ++iSide)
      {
        GeometricEntity* cell = face->getNeighborGeo(iSide);
        states[iSide] = cell->getStates();
        cellIDs[iSide] = cell->getID();
      }

      if (alpha[(*states[LEFT])[0]->getLocalID()] > 0.0 || alpha[(*states[RIGHT])[0]->getLocalID()] > 0.0)
      {
        computeCellGradVars(*states[LEFT ],m_sideGradVars[LEFT ]);
        computeCellGradVars(*states[RIGHT],m_sideGradVars[RIGHT]);

        for (CFuint iFlx = 0; iFlx < nbrFaceFlx; ++iFlx)
        {
          const CFuint flxIdxL = flxConn[LEFT ][iFlx];
          const CFuint flxIdxR = flxConn[RIGHT][iFlx];
          const CFuint solL = (*m_closestSolToFlx)[flxIdxL];
          const CFuint solR = (*m_closestSolToFlx)[flxIdxR];

          const CFuint idxL = cellIDs[LEFT ]*m_nbrCellFlxPnts + flxIdxL;
          const CFuint idxR = cellIDs[RIGHT]*m_nbrCellFlxPnts + flxIdxR;

          for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
          {
            const CFreal sample = 0.5*(m_sideGradVars[LEFT](iEq,solL) + m_sideGradVars[RIGHT](iEq,solR));
            m_faceSamples[idxL*m_nbrEqs + iEq] = sample;
            m_faceSamples[idxR*m_nbrEqs + iEq] = sample;
          }

          m_sampleStamp[idxL] = m_stamp;
          m_sampleStamp[idxR] = m_stamp;
        }
      }

      m_faceBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal DiffRHSJacobFluxReconstructionNSSubcellBlending::lineDerivative(const CFuint dir, const CFuint j, const CFuint i,
                                                                       const CFuint iEq, const CFuint cellID) const
{
  const CFuint n = m_nbrSolPnts1D;

  // solution point at position k of line j is base + k*stride (iSol = iKsi*n + iEta)
  const CFuint base   = (dir == KSI) ? j : j*n;
  const CFuint stride = (dir == KSI) ? n : 1;

  CFreal uLow, xLow, uHigh, xHigh;

  if (i > 0)
  {
    uLow = m_blendGradVars(iEq,base + (i-1)*stride);
    xLow = m_solPnts1D[i-1];
  }
  else
  {
    const CFuint idx = cellID*m_nbrCellFlxPnts + m_lineEndFlx[dir][0][j];
    if (m_sampleStamp[idx] == m_stamp)
    {
      uLow = m_faceSamples[idx*m_nbrEqs + iEq];
      xLow = -1.0;
    }
    else
    {
      uLow = m_blendGradVars(iEq,base + i*stride);
      xLow = m_solPnts1D[i];
    }
  }

  if (i + 1 < n)
  {
    uHigh = m_blendGradVars(iEq,base + (i+1)*stride);
    xHigh = m_solPnts1D[i+1];
  }
  else
  {
    const CFuint idx = cellID*m_nbrCellFlxPnts + m_lineEndFlx[dir][1][j];
    if (m_sampleStamp[idx] == m_stamp)
    {
      uHigh = m_faceSamples[idx*m_nbrEqs + iEq];
      xHigh = 1.0;
    }
    else
    {
      uHigh = m_blendGradVars(iEq,base + i*stride);
      xHigh = m_solPnts1D[i];
    }
  }

  return (uHigh - uLow)/(xHigh - xLow);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::prepareAlphaCells()
{
  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  CellToFaceGEBuilder::GeoData& geoDataCell = m_cellBuilders[LEFT]->getDataGE();
  geoDataCell.trs = cells;

  const CFuint n = m_nbrSolPnts1D;
  const CFuint nbrCells = cells->getLocalNbGeoEnts();

  for (CFuint cellID = 0; cellID < nbrCells; ++cellID)
  {
    geoDataCell.idx = cellID;
    GeometricEntity* cell = m_cellBuilders[LEFT]->buildGE();
    vector< State* >* states = cell->getStates();

    if (alpha[(*states)[0]->getLocalID()] > 0.0)
    {
      // subcell gradients of the points of the cell (central differences, chain-rule form of FR)
      computeCellGradVars(*states,m_blendGradVars);

      const std::valarray< CFreal > jacobDet = cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
        m_blendProjVects[iDim] = cell->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
      }

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        const CFuint iKsi = iSol/n;
        const CFuint iEta = iSol%n;
        const CFuint stateID = (*states)[iSol]->getLocalID();
        const CFreal invJacobDet = 1.0/jacobDet[iSol];

        vector< RealVector >& subGrad = m_subGrads[stateID];
        if (subGrad.size() != m_nbrEqs)
        {
          subGrad.assign(m_nbrEqs, RealVector(m_dim));
        }

        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          const CFreal dKsi = lineDerivative(KSI,iEta,iKsi,iEq,cellID);
          const CFreal dEta = lineDerivative(ETA,iKsi,iEta,iEq,cellID);
          for (CFuint iDim = 0; iDim < m_dim; ++iDim)
          {
            subGrad[iEq][iDim] = (m_blendProjVects[KSI][iSol][iDim]*dKsi + m_blendProjVects[ETA][iSol][iDim]*dEta)*invJacobDet;
          }
        }

        m_subGradStamp[stateID] = m_stamp;
      }

      // telescoped subcell normals (geometry, computed once per cell)
      CellSubcellData& cellData = m_cellData[cellID];
      if (cellData.normals.size() == 0)
      {
        m_scData.computeCellNormals(cell);
        cellData.normals = m_scData.getIntfNormals();
      }

      // unperturbed internal fluxes and R_sub of the cell
      computeCellSubcellRes(cellID,*states);
    }

    m_cellBuilders[LEFT]->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::computeCellSubcellRes(const CFuint cellID,
                                                                             const std::vector< State* >& states)
{
  CellSubcellData& cellData = m_cellData[cellID];
  const CFuint nbrIntf = m_intfSolL.size();

  if (cellData.intfFlux.size() != nbrIntf)
  {
    cellData.intfFlux.assign(nbrIntf, RealVector(m_nbrEqs));
    cellData.subRes.assign(m_nbrSolPnts, RealVector(m_nbrEqs));
  }

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    cellData.subRes[iSol] = 0.0;
  }

  // R_sub,i = (F_(i+1/2) - F_(i-1/2))/w_i: the flux enters the low side with + and the high side with -
  for (CFuint iIntf = 0; iIntf < nbrIntf; ++iIntf)
  {
    computeIntfFlux(cellID,states,iIntf,cellData.intfFlux[iIntf]);

    const RealVector& flux = cellData.intfFlux[iIntf];
    const CFuint solL = m_intfSolL[iIntf];
    const CFuint solR = m_intfSolR[iIntf];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      cellData.subRes[solL][iEq] += flux[iEq]/m_intfWidthL[iIntf];
      cellData.subRes[solR][iEq] -= flux[iEq]/m_intfWidthR[iIntf];
    }
  }

  cellData.stamp = m_stamp;
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::computeIntfFlux(const CFuint cellID,
                                                                       const std::vector< State* >& states,
                                                                       const CFuint iIntf, RealVector& flux)
{
  State* stateL = states[m_intfSolL[iIntf]];
  State* stateR = states[m_intfSolR[iIntf]];

  // radius of the interface (axisymmetric): mean of the two points
  const CFreal radius = m_axisymmetric ? 0.5*(stateL->getCoordinates()[YY] + stateR->getCoordinates()[YY]) : 0.;

  computeTwoPointFlux(*stateL,*stateR,m_subGrads[stateL->getLocalID()],m_subGrads[stateR->getLocalID()],
                      m_cellData[cellID].normals[iIntf],radius,flux);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::computeTwoPointFlux(State& sL, State& sR,
                                                                           const std::vector< RealVector >& gradL,
                                                                           const std::vector< RealVector >& gradR,
                                                                           const RealVector& normal, const CFreal radius,
                                                                           RealVector& flux)
{
  // face state: average of the update variables
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_avgSol[iEq] = 0.5*(sL[iEq] + sR[iEq]);
  }

  // gradient variables of the two points and of the face state
  m_twoStatePtrs[0] = sL.getData();
  m_twoStatePtrs[1] = sR.getData();
  m_twoStatePtrs[2] = &m_avgSol;
  m_diffusiveVarSet->setGradientVars(m_twoStatePtrs,m_twoGradVars,3);

  // unit vector e from L to R and distance d
  const RealVector& coordL = sL.getCoordinates();
  const RealVector& coordR = sR.getCoordinates();
  CFreal dist2 = 0.0;
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    m_twoDir[iDim] = coordR[iDim] - coordL[iDim];
    dist2 += m_twoDir[iDim]*m_twoDir[iDim];
  }
  const CFreal dist = std::sqrt(dist2);
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    m_twoDir[iDim] /= dist;
  }

  // g_f = gbar + (g2p - gbar.e) e
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    const CFreal uL = m_twoGradVars(iEq,0);
    const CFreal uR = m_twoGradVars(iEq,1);
    const CFreal uF = m_twoGradVars(iEq,2);

    const CFreal g2p = m_isLogGradVar[iEq] ? (cappedExp(uR - uF) - cappedExp(uL - uF))/dist : (uR - uL)/dist;

    CFreal gbarDotDir = 0.0;
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      gbarDotDir += 0.5*(gradL[iEq][iDim] + gradR[iEq][iDim])*m_twoDir[iDim];
    }

    RealVector& grad = *(m_twoGrad[iEq]);
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      grad[iDim] = 0.5*(gradL[iEq][iDim] + gradR[iEq][iDim]) + (g2p - gbarDotDir)*m_twoDir[iDim];
    }
  }

  prepareFluxComputation();
  computeFlux(m_avgSol,m_twoGrad,normal,radius,flux);
}

//////////////////////////////////////////////////////////////////////////////

const std::vector< RealVector >& DiffRHSJacobFluxReconstructionNSSubcellBlending::getPointGrad(State& state,
                                                                                               const std::vector< RealVector >& br2Grad) const
{
  const CFuint stateID = state.getLocalID();
  return (m_subGradStamp[stateID] == m_stamp) ? m_subGrads[stateID] : br2Grad;
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::computeInterfaceFlxCorrection()
{
  // gradient variables at the solution points of both cells and extrapolated to the flux points of the face
  computeCellGradVars(*m_states[LEFT],m_gradVarsSolPnts[LEFT]);
  computeCellGradVars(*m_states[RIGHT],m_gradVarsSolPnts[RIGHT]);

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlx];
    const CFuint flxIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx];

    extrapolateGradVarsToFlxPnt(m_gradVarsSolPnts[LEFT],(*m_flxSolDep)[flxIdxL],(*m_solPolyValsAtFlxPnts)[flxIdxL],m_nbrEqs,m_flxPntGradVars);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradVarsFlxPnt[LEFT](iEq,iFlx) = m_flxPntGradVars[iEq];
    }

    extrapolateGradVarsToFlxPnt(m_gradVarsSolPnts[RIGHT],(*m_flxSolDep)[flxIdxR],(*m_solPolyValsAtFlxPnts)[flxIdxR],m_nbrEqs,m_flxPntGradVars);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradVarsFlxPnt[RIGHT](iEq,iFlx) = m_flxPntGradVars[iEq];
    }
  }

  // compact BR2 gradients of the face, which depend on the two cells of the face only
  computeCompactBR2FaceGradients(m_gradVarsSolPnts[LEFT],m_gradVarsSolPnts[RIGHT],m_gradVarsFlxPnt[LEFT],m_gradVarsFlxPnt[RIGHT]);

  // weight of the two-point flux, single valued at the face
  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  const CFreal alphaF = std::max(alpha[(*m_states[LEFT])[0]->getLocalID()], alpha[(*m_states[RIGHT])[0]->getLocalID()]);

  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // average state and average compact face gradient
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;

      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0;
    }

    prepareFlxPntFluxComputation(iFlxPnt);

    // compute the common diffusive flux
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],m_flxPntRadius[iFlxPnt],m_flxPntRiemannFlux[iFlxPnt]);

    // blend with the two-point flux between the solution points closest to the flux point
    if (alphaF > 0.0)
    {
      const CFuint solL = (*m_closestSolToFlx)[(*m_faceFlxPntConnPerOrient)[m_orient][LEFT ][iFlxPnt]];
      const CFuint solR = (*m_closestSolToFlx)[(*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt]];
      State& stateL = *(*m_states[LEFT ])[solL];
      State& stateR = *(*m_states[RIGHT])[solR];

      computeTwoPointFlux(stateL,stateR,getPointGrad(stateL,*m_cellGrads[LEFT][solL]),getPointGrad(stateR,*m_cellGrads[RIGHT][solR]),
                          m_unitNormalFlxPnts[iFlxPnt],m_flxPntRadius[iFlxPnt],m_twoFlux);

      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_flxPntRiemannFlux[iFlxPnt][iEq] = (1.0 - alphaF)*m_flxPntRiemannFlux[iFlxPnt][iEq] + alphaF*m_twoFlux[iEq];
      }
    }

    // compute FI in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::computeDivDiscontFlx(std::vector< RealVector >& residuals)
{
  DiffRHSJacobFluxReconstructionNS::computeDivDiscontFlx(residuals);

  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  const CFuint stateID0 = (*m_cellStates)[0]->getLocalID();
  const CFreal a = alpha[stateID0];
  if (a <= 0.0) return;

  const CellSubcellData& cellData = m_cellData[m_cellOfState[stateID0]];
  cf_assert(cellData.stamp == m_stamp);

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      residuals[iSol][iEq] = (1.0 - a)*residuals[iSol][iEq] + a*cellData.subRes[iSol][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::computeDivDiscontFlxNeighb(RealVector& residuals, const CFuint side)
{
  DiffRHSJacobFluxReconstructionNS::computeDivDiscontFlxNeighb(residuals,side);

  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  const std::vector< State* >& states = *m_states[side];
  const CFuint stateID0 = states[0]->getLocalID();
  const CFreal a = alpha[stateID0];
  if (a <= 0.0) return;

  const CFuint cellID = m_cellOfState[stateID0];
  const CellSubcellData& cellData = m_cellData[cellID];
  cf_assert(cellData.stamp == m_stamp);

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFuint idx = m_nbrEqs*iSol + iEq;
      residuals[idx] = (1.0 - a)*residuals[idx] + a*cellData.subRes[iSol][iEq];
    }
  }

  // the perturbed point of this cell: change of the internal fluxes it touches
  if (side == m_pertSide)
  {
    const std::vector< CFuint >& pertIntfs = m_intfOfSol[m_pertSol];
    for (CFuint k = 0; k < pertIntfs.size(); ++k)
    {
      const CFuint iIntf = pertIntfs[k];
      computeIntfFlux(cellID,states,iIntf,m_twoFluxPert);

      const CFuint solL = m_intfSolL[iIntf];
      const CFuint solR = m_intfSolR[iIntf];
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        const CFreal dFlux = a*(m_twoFluxPert[iEq] - cellData.intfFlux[iIntf][iEq]);
        residuals[m_nbrEqs*solL + iEq] += dFlux/m_intfWidthL[iIntf];
        residuals[m_nbrEqs*solR + iEq] -= dFlux/m_intfWidthR[iIntf];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::computeCorrection(CFuint side, vector< RealVector >& corrections)
{
  // FR correction, sets m_cellStates to the states of side
  DiffRHSJacobFluxReconstructionNS::computeCorrection(side,corrections);

  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  const CFreal a = alpha[(*m_cellStates)[0]->getLocalID()];
  if (a <= 0.0) return;

  // alpha F_f (c_f/w_f at the closest point - divh(., f)) for each face flux point f
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlx];
    const RealVector& faceFlux = m_cellFlx[side][iFlx];

    const CFuint closest = (*m_closestSolToFlx)[flxIdx];
    const CFreal factor = a*m_faceFactors[flxIdx];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      corrections[closest][iEq] += factor*faceFlux[iEq];
    }

    const std::vector< CFuint >& solDep = (*m_flxSolDep)[flxIdx];
    for (CFuint k = 0; k < solDep.size(); ++k)
    {
      const CFuint solIdx = solDep[k];
      const CFreal divh = a*m_corrFctDiv[solIdx][flxIdx];
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        corrections[solIdx][iEq] -= divh*faceFlux[iEq];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::computeFaceSubcellFactors(FluxReconstructionElementData* frData,
                                                                                 const std::vector< std::vector< CFreal > >& corrFctDiv,
                                                                                 std::vector< CFreal >& factors)
{
  const CFuint n = frData->getSolPntsLocalCoord1D()->size();
  const std::vector< CFreal > widths = SubcellBlendingQuadData::computeSubcellWidths1D(frData);
  SafePtr< vector< CFuint > > closestSol = frData->getClosestSolToFlxIdx();
  SafePtr< vector< CFuint > > flxPntFlxDim = frData->getFluxPntFluxDim();
  const CFuint nbrFlxPnts = closestSol->size();
  const CFuint nbrSolPnts = n*n;

  factors.resize(nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    const CFuint sol  = (*closestSol)[iFlx];
    const CFuint dir  = (*flxPntFlxDim)[iFlx];
    const CFuint pos  = (dir == KSI) ? sol/n : sol%n;
    const CFuint line = (dir == KSI) ? sol%n : sol/n;

    // integral of the correction along the line of the flux point, and its size elsewhere
    CFreal lineIntegral = 0.0;
    CFreal offLineMax = 0.0;
    CFreal maxAbs = 0.0;
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      const CFreal divh = corrFctDiv[iSol][iFlx];
      const CFuint solPos  = (dir == KSI) ? iSol/n : iSol%n;
      const CFuint solLine = (dir == KSI) ? iSol%n : iSol/n;
      maxAbs = std::max(maxAbs, std::abs(divh));
      if (solLine == line)
      {
        lineIntegral += widths[solPos]*divh;
      }
      else
      {
        offLineMax = std::max(offLineMax, std::abs(divh));
      }
    }

    if (std::abs(std::abs(lineIntegral) - 1.0) > 1e-8 || offLineMax > 1e-10*maxAbs)
    {
      throw BadValueException(FromHere(), "DiffRHSJacobNSSubcellBlending: the correction function of a flux point is not "
                              "a 1D correction of its line with unit integral (tensor-product quads expected).");
    }

    factors[iFlx] = lineIntegral/widths[pos];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::setup()
{
  CFAUTOTRACE;

  DiffRHSJacobFluxReconstructionNS::setup();

  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() == 1);
  FluxReconstructionElementData* frData = frLocalData[0];

  if (frData->getShape() != CFGeoShape::QUAD)
  {
    throw NotImplementedException(FromHere(),
      "DiffRHSJacobNSSubcellBlending: subcell diffusion is only implemented for quadrilateral elements.");
  }

  // subcell grid with the telescoped normals of the convective subcell blending
  m_scData.setup(frData,m_dim,m_nbrEqs);
  m_scData.setCorrectionFunction(m_corrFctDiv);

  m_solPnts1D = *(frData->getSolPntsLocalCoord1D());
  m_nbrSolPnts1D = m_solPnts1D.size();
  m_widths1D = m_scData.getWidths1D();
  const CFuint n = m_nbrSolPnts1D;

  if (n < 2)
  {
    throw BadValueException(FromHere(), "DiffRHSJacobNSSubcellBlending: needs at least two solution points per direction (P1).");
  }

  // the line indexing assumes iSol = iKsi*n + iEta
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    if (std::abs((*m_solPntsLocalCoords)[iSol][KSI] - m_solPnts1D[iSol/n]) > 1e-12 ||
        std::abs((*m_solPntsLocalCoords)[iSol][ETA] - m_solPnts1D[iSol%n]) > 1e-12)
    {
      throw BadValueException(FromHere(), "DiffRHSJacobNSSubcellBlending: solution points are not numbered iKsi*n + iEta.");
    }
  }

  // cell flux point at each end of each reference line
  m_closestSolToFlx = frData->getClosestSolToFlxIdx();
  SafePtr< vector< CFuint > > flxPntFlxDim = frData->getFluxPntFluxDim();
  SafePtr< vector< RealVector > > flxPntCoords = frData->getFlxPntsLocalCoords();
  m_nbrCellFlxPnts = flxPntCoords->size();

  m_lineEndFlx.assign(2, vector< vector< CFuint > >(2, vector< CFuint >(n, m_nbrCellFlxPnts)));
  for (CFuint iFlx = 0; iFlx < m_nbrCellFlxPnts; ++iFlx)
  {
    const CFuint sol = (*m_closestSolToFlx)[iFlx];
    const CFuint dir = (*flxPntFlxDim)[iFlx];
    const CFuint end = ((*flxPntCoords)[iFlx][dir] > 0.0) ? 1 : 0;
    const CFuint pos  = (dir == KSI) ? sol/n : sol%n;
    const CFuint line = (dir == KSI) ? sol%n : sol/n;

    if (dir > ETA || pos != ((end == 0) ? 0 : n-1))
    {
      throw BadValueException(FromHere(), "DiffRHSJacobNSSubcellBlending: closest solution point of a flux point is not at the end of its line.");
    }
    m_lineEndFlx[dir][end][line] = iFlx;
  }

  // subcell boundary factor c_f/w_f of each flux point
  computeFaceSubcellFactors(frData,m_corrFctDiv,m_faceFactors);

  // internal interfaces, in the order of SubcellBlendingQuadData
  const CFuint nbrIntfPerDir = (n-1)*n;
  m_intfSolL.resize(2*nbrIntfPerDir);
  m_intfSolR.resize(2*nbrIntfPerDir);
  m_intfWidthL.resize(2*nbrIntfPerDir);
  m_intfWidthR.resize(2*nbrIntfPerDir);
  for (CFuint iKsi = 0; iKsi + 1 < n; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < n; ++iEta)
    {
      const CFuint iIntf = iKsi*n + iEta;
      m_intfSolL[iIntf] = iKsi*n + iEta;
      m_intfSolR[iIntf] = m_intfSolL[iIntf] + n;
      m_intfWidthL[iIntf] = m_widths1D[iKsi];
      m_intfWidthR[iIntf] = m_widths1D[iKsi+1];
    }
  }
  for (CFuint iKsi = 0; iKsi < n; ++iKsi)
  {
    for (CFuint iEta = 0; iEta + 1 < n; ++iEta)
    {
      const CFuint iIntf = nbrIntfPerDir + iKsi*(n-1) + iEta;
      m_intfSolL[iIntf] = iKsi*n + iEta;
      m_intfSolR[iIntf] = m_intfSolL[iIntf] + 1;
      m_intfWidthL[iIntf] = m_widths1D[iEta];
      m_intfWidthR[iIntf] = m_widths1D[iEta+1];
    }
  }
  m_intfOfSol.assign(m_nbrSolPnts, vector< CFuint >());
  for (CFuint iIntf = 0; iIntf < m_intfSolL.size(); ++iIntf)
  {
    m_intfOfSol[m_intfSolL[iIntf]].push_back(iIntf);
    m_intfOfSol[m_intfSolR[iIntf]].push_back(iIntf);
  }

  // gradient variables that are logarithms: ln y_i, ln T, ln Tv of LogRhoivLogTTv
  m_isLogGradVar.assign(m_nbrEqs, false);
  if (getMethodData().getUpdateVarStr() == "LogRhoivLogTTv")
  {
    SafePtr< MultiScalarTerm< EulerTerm > > term =
      PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm().d_castTo< MultiScalarTerm< EulerTerm > >();
    const CFuint nbSpecies = term->getNbScalarVars(0);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_isLogGradVar[iEq] = (iEq < nbSpecies) || (iEq >= nbSpecies + m_dim);
    }
  }
  std::string logVars;
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    if (m_isLogGradVar[iEq]) logVars += " " + StringOps::to_str(iEq);
  }
  CFLog(INFO, "DiffRHSJacobNSSubcellBlending: two-point gradients of logarithms for the gradient variables" << logVars << "\n");

  // cell of each local state
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbrCells = cells->getLocalNbGeoEnts();
  CellToFaceGEBuilder::GeoData& geoDataCell = m_cellBuilders[LEFT]->getDataGE();
  geoDataCell.trs = cells;
  for (CFuint cellID = 0; cellID < nbrCells; ++cellID)
  {
    geoDataCell.idx = cellID;
    GeometricEntity* cell = m_cellBuilders[LEFT]->buildGE();
    const vector< State* >& states = *cell->getStates();
    for (CFuint iSol = 0; iSol < states.size(); ++iSol)
    {
      const CFuint stateID = states[iSol]->getLocalID();
      if (stateID >= m_cellOfState.size()) m_cellOfState.resize(stateID+1, 0);
      m_cellOfState[stateID] = cellID;
    }
    m_cellBuilders[LEFT]->releaseGE();
  }
  const CFuint nbrStates = m_cellOfState.size();

  m_faceSamples.assign(nbrCells*m_nbrCellFlxPnts*m_nbrEqs, 0.0);
  m_sampleStamp.assign(nbrCells*m_nbrCellFlxPnts, 0);
  m_subGrads.assign(nbrStates, vector< RealVector >());
  m_subGradStamp.assign(nbrStates, 0);
  m_cellData.assign(nbrCells, CellSubcellData());
  for (CFuint cellID = 0; cellID < nbrCells; ++cellID)
  {
    m_cellData[cellID].stamp = 0;
  }
  m_stamp = 0;

  m_blendGradVars.resize(m_nbrEqs,m_nbrSolPnts);
  m_sideGradVars.resize(2);
  m_sideGradVars[LEFT ].resize(m_nbrEqs,m_nbrSolPnts);
  m_sideGradVars[RIGHT].resize(m_nbrEqs,m_nbrSolPnts);
  m_blendProjVects.resize(m_dim);

  m_twoStatePtrs.resize(3);
  m_twoGradVars.resize(m_nbrEqs,3);
  m_twoGrad.resize(m_nbrEqs);
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_twoGrad[iEq] = new RealVector(m_dim);
  }
  m_twoDir.resize(m_dim);
  m_twoFlux.resize(m_nbrEqs);
  m_twoFluxPert.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstructionNSSubcellBlending::unsetup()
{
  for (CFuint iEq = 0; iEq < m_twoGrad.size(); ++iEq)
  {
    deletePtr(m_twoGrad[iEq]);
  }
  m_twoGrad.clear();

  DiffRHSJacobFluxReconstructionNS::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
