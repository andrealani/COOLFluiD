// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>

#include "Common/BadValueException.hh"
#include "Common/CFLog.hh"
#include "Common/NotImplementedException.hh"

#include "FluxReconstructionMethod/SubcellBlendingQuadData.hh"
#include "FluxReconstructionMethod/TensorProductGaussIntegrator.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

SubcellBlendingQuadData::SubcellBlendingQuadData() :
  m_dim(0),
  m_nbrEqs(0),
  m_nbrSolPnts1D(0),
  m_widths1D(),
  m_intfSolL(),
  m_intfSolR(),
  m_intfWidthL(),
  m_intfWidthR(),
  m_intfCoords(),
  m_intfPlaneIdx(),
  m_intfNormals(),
  m_intfFlux(),
  m_intfOfSol(),
  m_subcellRes(),
  m_intfFluxPert(),
  m_intfFluxDiff(),
  m_closestSolToFlx(CFNULL),
  m_flxPntFlxDim(CFNULL),
  m_unitNormal(),
  m_closeSubcells(false),
  m_warnedOpenElement(false),
  m_intfTransWidth(),
  m_extCoords(),
  m_extPlaneIdx(),
  m_cosine(),
  m_eigenvalue(),
  m_closureDefect(),
  m_closurePotential(),
  m_closureCoef()
{
}

//////////////////////////////////////////////////////////////////////////////

SubcellBlendingQuadData::~SubcellBlendingQuadData()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< CFreal > SubcellBlendingQuadData::computeSubcellWidths1D(FluxReconstructionElementData* frData)
{
  SafePtr< vector< CFreal > > solPnts1D = frData->getSolPntsLocalCoord1D();
  const CFuint nbrSolPnts1D = solPnts1D->size();

  // Gauss quadrature on [-1,1], exact for the degree p Lagrange polynomials
  TensorProductGaussIntegrator tpIntegrator(DIM_1D, frData->getPolyOrder());
  vector< RealVector > nodeCoord(2);
  nodeCoord[0].resize(1);
  nodeCoord[0][KSI] = -1.0;
  nodeCoord[1].resize(1);
  nodeCoord[1][KSI] = +1.0;
  const vector< RealVector > quadPntCoords  = tpIntegrator.getQuadPntsCoords  (nodeCoord);
  const vector< CFreal >     quadPntWeights = tpIntegrator.getQuadPntsWheights(nodeCoord);
  const CFuint nbrQPnts = quadPntCoords.size();
  cf_assert(quadPntWeights.size() == nbrQPnts);

  // width of subcell i = integral over [-1,1] of the Lagrange polynomial of solution point i
  vector< CFreal > widths(nbrSolPnts1D, 0.0);
  CFreal widthSum = 0.0;
  for (CFuint iSol = 0; iSol < nbrSolPnts1D; ++iSol)
  {
    const CFreal ksiSol = (*solPnts1D)[iSol];
    for (CFuint iQPnt = 0; iQPnt < nbrQPnts; ++iQPnt)
    {
      const CFreal ksiQPnt = quadPntCoords[iQPnt][KSI];
      CFreal lagrangeVal = 1.0;
      for (CFuint iFac = 0; iFac < nbrSolPnts1D; ++iFac)
      {
        if (iFac != iSol)
        {
          const CFreal ksiFac = (*solPnts1D)[iFac];
          lagrangeVal *= (ksiQPnt - ksiFac)/(ksiSol - ksiFac);
        }
      }
      widths[iSol] += quadPntWeights[iQPnt]*lagrangeVal;
    }

    if (widths[iSol] <= 0.0)
    {
      throw BadValueException(FromHere(),
        "Subcell blending: the solution point distribution has a non-positive quadrature weight, "
        "so it cannot define subcells. Use GaussLegendre or Lobatto solution points.");
    }
    widthSum += widths[iSol];
  }

  if (std::abs(widthSum - 2.0) > 1.0e-10)
  {
    throw BadValueException(FromHere(),
      "Subcell blending: the 1D subcell widths do not sum to the reference length 2.");
  }

  return widths;
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::setup(FluxReconstructionElementData* frData,
                                    const CFuint dim, const CFuint nbrEqs)
{
  if (frData->getShape() != CFGeoShape::QUAD)
  {
    throw NotImplementedException(FromHere(),
      "SubcellBlendingQuadData: subcell blending is only implemented for quadrilateral elements.");
  }

  m_dim    = dim;
  m_nbrEqs = nbrEqs;

  m_closestSolToFlx = frData->getClosestSolToFlxIdx();
  m_flxPntFlxDim    = frData->getFluxPntFluxDim();

  SafePtr< vector< CFreal > > solPnts1D = frData->getSolPntsLocalCoord1D();
  m_nbrSolPnts1D = solPnts1D->size();
  m_widths1D     = computeSubcellWidths1D(frData);

  const CFuint nbr1D     = m_nbrSolPnts1D;
  const CFuint nbrSolPnts = nbr1D*nbr1D;

  // subcell boundaries in 1D: -1, -1 + w_0, -1 + w_0 + w_1, ..., +1
  vector< CFreal > bnds1D(nbr1D + 1);
  bnds1D[0] = -1.0;
  for (CFuint i = 0; i < nbr1D; ++i)
  {
    bnds1D[i+1] = bnds1D[i] + m_widths1D[i];
  }

  const CFuint nbrIntfPerDir = (nbr1D - 1)*nbr1D;
  const CFuint nbrIntf       = 2*nbrIntfPerDir;

  m_intfSolL.resize(nbrIntf);
  m_intfSolR.resize(nbrIntf);
  m_intfWidthL.resize(nbrIntf);
  m_intfWidthR.resize(nbrIntf);
  m_intfPlaneIdx.resize(nbrIntf);
  m_intfCoords.assign(nbrIntf, RealVector(2));
  m_intfNormals.assign(nbrIntf, RealVector(m_dim));
  m_intfFlux.assign(nbrIntf, RealVector(m_nbrEqs));
  m_subcellRes.assign(nbrSolPnts, RealVector(m_nbrEqs));

  // interfaces normal to ksi, between solution points (iKsi, iEta) and (iKsi+1, iEta)
  for (CFuint iKsi = 0; iKsi + 1 < nbr1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbr1D; ++iEta)
    {
      const CFuint iIntf = iKsi*nbr1D + iEta;
      m_intfSolL[iIntf]     = iKsi*nbr1D + iEta;
      m_intfSolR[iIntf]     = m_intfSolL[iIntf] + nbr1D;
      m_intfWidthL[iIntf]   = m_widths1D[iKsi];
      m_intfWidthR[iIntf]   = m_widths1D[iKsi+1];
      m_intfPlaneIdx[iIntf] = KSI;
      m_intfCoords[iIntf][KSI] = bnds1D[iKsi+1];
      m_intfCoords[iIntf][ETA] = (*solPnts1D)[iEta];
    }
  }

  // interfaces normal to eta, between solution points (iKsi, iEta) and (iKsi, iEta+1)
  for (CFuint iKsi = 0; iKsi < nbr1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta + 1 < nbr1D; ++iEta)
    {
      const CFuint iIntf = nbrIntfPerDir + iKsi*(nbr1D-1) + iEta;
      m_intfSolL[iIntf]     = iKsi*nbr1D + iEta;
      m_intfSolR[iIntf]     = m_intfSolL[iIntf] + 1;
      m_intfWidthL[iIntf]   = m_widths1D[iEta];
      m_intfWidthR[iIntf]   = m_widths1D[iEta+1];
      m_intfPlaneIdx[iIntf] = ETA;
      m_intfCoords[iIntf][KSI] = (*solPnts1D)[iKsi];
      m_intfCoords[iIntf][ETA] = bnds1D[iEta+1];
    }
  }

  // interfaces adjacent to each solution point, at most two per reference direction
  m_intfOfSol.assign(nbrSolPnts, vector< CFuint >());
  CFuint maxIntfPerSol = 0;
  for (CFuint iIntf = 0; iIntf < nbrIntf; ++iIntf)
  {
    m_intfOfSol[m_intfSolL[iIntf]].push_back(iIntf);
    m_intfOfSol[m_intfSolR[iIntf]].push_back(iIntf);
  }
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    maxIntfPerSol = std::max(maxIntfPerSol, static_cast<CFuint>(m_intfOfSol[iSol].size()));
  }

  m_intfFluxPert.assign(std::max(maxIntfPerSol, static_cast<CFuint>(1)), RealVector(m_nbrEqs));
  m_intfFluxDiff.resize(m_nbrEqs);
  m_unitNormal.resize(m_dim);

  // --- data for closing the subcells on curved cells (see closeSubcells) ---

  // transverse width of each internal interface
  m_intfTransWidth.resize(nbrIntf);
  for (CFuint iIntf = 0; iIntf < nbrIntf; ++iIntf)
  {
    const CFuint solL = m_intfSolL[iIntf];
    const CFuint transIdx = (m_intfPlaneIdx[iIntf] == KSI) ? solL%nbr1D : solL/nbr1D;
    m_intfTransWidth[iIntf] = m_widths1D[transIdx];
  }

  // the exterior subcell faces are the element flux points: this needs them at the
  // transverse solution point coordinates (same 1D distribution for both)
  SafePtr< vector< CFreal > > flxPnts1D = frData->getFlxPntsLocalCoord1D();
  m_closeSubcells = (flxPnts1D->size() == nbr1D);
  for (CFuint i = 0; m_closeSubcells && i < nbr1D; ++i)
  {
    m_closeSubcells = std::abs((*flxPnts1D)[i] - (*solPnts1D)[i]) < 1.0e-12;
  }
  if (!m_closeSubcells)
  {
    CFLog(WARN, "SubcellBlendingQuadData: flux and solution points differ in 1D, "
          "the subcells are not closed on curved cells\n");
  }

  // exterior faces ksi = -1, ksi = +1, eta = -1, eta = +1, at the transverse solution points
  m_extCoords.assign(4*nbr1D, RealVector(2));
  m_extPlaneIdx.resize(4*nbr1D);
  for (CFuint i = 0; i < nbr1D; ++i)
  {
    m_extCoords[i][KSI]           = -1.0; m_extCoords[i][ETA]           = (*solPnts1D)[i];
    m_extCoords[nbr1D+i][KSI]     = +1.0; m_extCoords[nbr1D+i][ETA]     = (*solPnts1D)[i];
    m_extCoords[2*nbr1D+i][KSI]   = (*solPnts1D)[i]; m_extCoords[2*nbr1D+i][ETA] = -1.0;
    m_extCoords[3*nbr1D+i][KSI]   = (*solPnts1D)[i]; m_extCoords[3*nbr1D+i][ETA] = +1.0;
    m_extPlaneIdx[i] = m_extPlaneIdx[nbr1D+i] = KSI;
    m_extPlaneIdx[2*nbr1D+i] = m_extPlaneIdx[3*nbr1D+i] = ETA;
  }

  // eigenvectors of the path graph Laplacian with nbr1D nodes: cos(pi k (i+1/2)/n),
  // eigenvalues 2 - 2 cos(pi k/n). The subcell grid Laplacian is the tensor product.
  const CFreal pi = std::acos(-1.0);
  m_cosine.assign(nbr1D, vector< CFreal >(nbr1D));
  m_eigenvalue.resize(nbr1D);
  for (CFuint k = 0; k < nbr1D; ++k)
  {
    const CFreal norm = std::sqrt((k == 0 ? 1.0 : 2.0)/nbr1D);
    m_eigenvalue[k] = 2.0 - 2.0*std::cos(pi*k/nbr1D);
    for (CFuint i = 0; i < nbr1D; ++i)
    {
      m_cosine[k][i] = norm*std::cos(pi*k*(i + 0.5)/nbr1D);
    }
  }

  m_closureDefect.assign(nbrSolPnts, RealVector(m_dim));
  m_closurePotential.assign(nbrSolPnts, RealVector(m_dim));
  m_closureCoef.resize(m_dim);
}

//////////////////////////////////////////////////////////////////////////////

CFreal SubcellBlendingQuadData::getFaceSubcellWidth(const CFuint flxIdx) const
{
  const CFuint solIdx = (*m_closestSolToFlx)[flxIdx];
  const CFuint flxDim = (*m_flxPntFlxDim)[flxIdx];

  // solution points are numbered iSol = iKsi*nbrSolPnts1D + iEta
  const CFuint idx1D = (flxDim == KSI) ? solIdx/m_nbrSolPnts1D : solIdx%m_nbrSolPnts1D;

  return m_widths1D[idx1D];
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::computeCellNormals(GeometricEntity* cell)
{
  // the plane index is given per point, so both reference directions are done in one call
  m_intfNormals = cell->computeMappedCoordPlaneNormalAtMappedCoords(m_intfPlaneIdx, m_intfCoords);

  // on a curved cell these sampled normals do not close the subcells
  if (m_closeSubcells && m_nbrSolPnts1D > 1) closeSubcells(cell);
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::closeSubcells(GeometricEntity* cell)
{
  // Area vector of a subcell face = metric normal x transverse width. Subcell s is closed
  // if its outward area vectors sum to zero. With B the interface/subcell incidence
  // (+1 on the low side, -1 on the high side), a the internal area vectors and b the
  // exterior ones, the defect is r = B a + b. The smallest change of a that gives
  // B a + b = 0 is a -= B^T phi with B B^T phi = r, B B^T being the Laplacian of the
  // subcell grid, solved exactly with its cosine eigenvectors.
  const CFuint n = m_nbrSolPnts1D;
  const CFuint nbrSol = n*n;
  const CFuint nbrIntf = m_intfSolL.size();

  // exterior faces: the same metric as the element flux points
  const vector< RealVector > extNormals =
    cell->computeMappedCoordPlaneNormalAtMappedCoords(m_extPlaneIdx, m_extCoords);

  for (CFuint iSol = 0; iSol < nbrSol; ++iSol) m_closureDefect[iSol] = 0.0;

  // exterior contributions, outward: minus the plane normal on the low faces
  CFreal extScale = 0.0;
  for (CFuint i = 0; i < n; ++i)
  {
    const CFreal w = m_widths1D[i];
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_closureDefect[i][iDim]             -= w*extNormals[i][iDim];        // ksi = -1
      m_closureDefect[(n-1)*n+i][iDim]     += w*extNormals[n+i][iDim];      // ksi = +1
      m_closureDefect[i*n][iDim]           -= w*extNormals[2*n+i][iDim];    // eta = -1
      m_closureDefect[i*n+n-1][iDim]       += w*extNormals[3*n+i][iDim];    // eta = +1
      extScale += w*(std::abs(extNormals[i][iDim]) + std::abs(extNormals[n+i][iDim]) +
                     std::abs(extNormals[2*n+i][iDim]) + std::abs(extNormals[3*n+i][iDim]));
    }
  }

  // internal contributions: outward for the low side, inward for the high side
  for (CFuint iIntf = 0; iIntf < nbrIntf; ++iIntf)
  {
    const CFreal w = m_intfTransWidth[iIntf];
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_closureDefect[m_intfSolL[iIntf]][iDim] += w*m_intfNormals[iIntf][iDim];
      m_closureDefect[m_intfSolR[iIntf]][iDim] -= w*m_intfNormals[iIntf][iDim];
    }
  }

  // the internal terms cancel in the sum over subcells, so the total is the closure of the
  // element boundary itself; internal changes cannot fix that part (constant mode, skipped below)
  if (!m_warnedOpenElement)
  {
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      CFreal total = 0.0;
      for (CFuint iSol = 0; iSol < nbrSol; ++iSol) total += m_closureDefect[iSol][iDim];
      if (std::abs(total) > 1.0e-10*extScale)
      {
        CFLog(WARN, "SubcellBlendingQuadData: the exterior faces of an element do not close "
              "(relative " << std::abs(total)/extScale << "), its subcells are closed up to that\n");
        m_warnedOpenElement = true;
        break;
      }
    }
  }

  // phi = sum over modes (k,l) != (0,0) of <r, e_kl>/(lambda_k + lambda_l) e_kl
  for (CFuint iSol = 0; iSol < nbrSol; ++iSol) m_closurePotential[iSol] = 0.0;
  for (CFuint k = 0; k < n; ++k)
  {
    for (CFuint l = 0; l < n; ++l)
    {
      if (k == 0 && l == 0) continue;

      m_closureCoef = 0.0;
      for (CFuint i = 0; i < n; ++i)
      {
        for (CFuint j = 0; j < n; ++j)
        {
          const CFreal e = m_cosine[k][i]*m_cosine[l][j];
          for (CFuint iDim = 0; iDim < m_dim; ++iDim)
          {
            m_closureCoef[iDim] += e*m_closureDefect[i*n+j][iDim];
          }
        }
      }
      m_closureCoef /= (m_eigenvalue[k] + m_eigenvalue[l]);

      for (CFuint i = 0; i < n; ++i)
      {
        for (CFuint j = 0; j < n; ++j)
        {
          const CFreal e = m_cosine[k][i]*m_cosine[l][j];
          for (CFuint iDim = 0; iDim < m_dim; ++iDim)
          {
            m_closurePotential[i*n+j][iDim] += e*m_closureCoef[iDim];
          }
        }
      }
    }
  }

  // a -= B^T phi, back from area vector to normal by the transverse width
  for (CFuint iIntf = 0; iIntf < nbrIntf; ++iIntf)
  {
    const CFuint solL = m_intfSolL[iIntf];
    const CFuint solR = m_intfSolR[iIntf];
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_intfNormals[iIntf][iDim] -=
        (m_closurePotential[solL][iDim] - m_closurePotential[solR][iDim])/m_intfTransWidth[iIntf];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::computeIntfFlux(const CFuint iIntf,
                                              const std::vector< State* >& states,
                                              RiemannFlux& riemannFlux,
                                              RealVector& result)
{
  const RealVector& metricNormal = m_intfNormals[iIntf];

  CFreal normalSize2 = 0.0;
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    normalSize2 += metricNormal[iDim]*metricNormal[iDim];
  }
  const CFreal normalSize = std::sqrt(normalSize2);

  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    m_unitNormal[iDim] = metricNormal[iDim]/normalSize;
  }

  const RealVector& flux = riemannFlux.computeFlux(*(states[m_intfSolL[iIntf]]),
                                                   *(states[m_intfSolR[iIntf]]),
                                                   m_unitNormal);

  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    result[iEq] = flux[iEq]*normalSize;
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::computeSubcellRes(const CFreal alpha,
                                                const std::vector< State* >& states,
                                                RiemannFlux& riemannFlux)
{
  const CFuint nbrSolPnts = m_subcellRes.size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_subcellRes[iSol] = 0.0;
  }

  const CFuint nbrIntf = m_intfSolL.size();
  for (CFuint iIntf = 0; iIntf < nbrIntf; ++iIntf)
  {
    computeIntfFlux(iIntf, states, riemannFlux, m_intfFlux[iIntf]);

    const RealVector& flux = m_intfFlux[iIntf];
    const CFuint solL = m_intfSolL[iIntf];
    const CFuint solR = m_intfSolR[iIntf];
    const CFreal factorL = alpha/m_intfWidthL[iIntf];
    const CFreal factorR = alpha/m_intfWidthR[iIntf];

    // the flux leaves the subcell on the low side and enters the one on the high side
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_subcellRes[solL][iEq] -= factorL*flux[iEq];
      m_subcellRes[solR][iEq] += factorR*flux[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::addSubcellResDelta(const CFreal alpha, const CFuint pertSol,
                                                 const std::vector< State* >& states,
                                                 RiemannFlux& riemannFlux,
                                                 RealVector& res)
{
  const std::vector< CFuint >& pertIntfs = m_intfOfSol[pertSol];
  const CFuint nbrPertIntfs = pertIntfs.size();

  for (CFuint k = 0; k < nbrPertIntfs; ++k)
  {
    const CFuint iIntf = pertIntfs[k];

    computeIntfFlux(iIntf, states, riemannFlux, m_intfFluxPert[k]);

    const RealVector& fluxNew = m_intfFluxPert[k];
    const RealVector& fluxOld = m_intfFlux[iIntf];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_intfFluxDiff[iEq] = fluxNew[iEq] - fluxOld[iEq];
    }

    const CFuint solL = m_intfSolL[iIntf];
    const CFuint solR = m_intfSolR[iIntf];
    const CFreal factorL = alpha/m_intfWidthL[iIntf];
    const CFreal factorR = alpha/m_intfWidthR[iIntf];

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      res[m_nbrEqs*solL+iEq] -= factorL*m_intfFluxDiff[iEq];
      res[m_nbrEqs*solR+iEq] += factorR*m_intfFluxDiff[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
