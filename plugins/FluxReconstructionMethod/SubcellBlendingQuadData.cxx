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
  m_telescope(false),
  m_derivMat1D(),
  m_lagrangeAtEnds(),
  m_corrDerivL(),
  m_corrDerivR(),
  m_metricCoords(),
  m_metricPlaneIdx(),
  m_faceSum(),
  m_metricDeriv(),
  m_endJumpL(),
  m_endJumpR(),
  m_linear(false),
  m_cellLinear(true),
  m_frozenPhi(CFNULL),
  m_frozenValid(CFNULL),
  m_recordPhi(false),
  m_limiter(0),
  m_limiterEps(0.),
  m_updateVarSet(CFNULL),
  m_solPnts1D(),
  m_bnds1D(),
  m_lineFlx(),
  m_slopes(),
  m_reconstructed(),
  m_pertSlopes(),
  m_cellSamples(CFNULL),
  m_recStateL(CFNULL),
  m_recStateR(CFNULL),
  m_recNodeL(CFNULL),
  m_recNodeR(CFNULL),
  m_testState(),
  m_toRec(CFNULL),
  m_fromRec(CFNULL),
  m_transInState(CFNULL),
  m_recPrev(),
  m_recCur(),
  m_recNext(),
  m_recFace(),
  m_pointSlope()
{
}

//////////////////////////////////////////////////////////////////////////////

SubcellBlendingQuadData::~SubcellBlendingQuadData()
{
  // the states own their coordinates (Node built with isOnMesh = false)
  delete m_recStateL;
  delete m_recStateR;
  delete m_transInState;
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
  m_faceSum.resize(m_dim);
  m_metricDeriv.resize(m_dim);
  m_endJumpL.resize(m_dim);
  m_endJumpR.resize(m_dim);

  // --- data for the telescoped subcell normals (see telescopeNormals) ---
  // the correction function arrives with setCorrectionFunction; until then the normals are sampled
  m_telescope = false;

  // derivative of the Lagrange polynomial of point m at point k, and the polynomials at -1, +1
  m_derivMat1D.assign(nbr1D, vector< CFreal >(nbr1D, 0.0));
  m_lagrangeAtEnds.assign(2, vector< CFreal >(nbr1D, 1.0));
  for (CFuint m = 0; m < nbr1D; ++m)
  {
    const CFreal xm = (*solPnts1D)[m];
    for (CFuint q = 0; q < nbr1D; ++q)
    {
      if (q == m) continue;
      m_lagrangeAtEnds[0][m] *= (-1.0 - (*solPnts1D)[q])/(xm - (*solPnts1D)[q]);
      m_lagrangeAtEnds[1][m] *= (+1.0 - (*solPnts1D)[q])/(xm - (*solPnts1D)[q]);
    }
    for (CFuint k = 0; k < nbr1D; ++k)
    {
      const CFreal xk = (*solPnts1D)[k];
      CFreal deriv = 0.0;
      for (CFuint l = 0; l < nbr1D; ++l)
      {
        if (l == m) continue;
        CFreal term = 1.0/(xm - (*solPnts1D)[l]);
        for (CFuint q = 0; q < nbr1D; ++q)
        {
          if (q != m && q != l) term *= (xk - (*solPnts1D)[q])/(xm - (*solPnts1D)[q]);
        }
        deriv += term;
      }
      m_derivMat1D[k][m] = deriv;
    }
  }

  // plane normals needed per cell, in one list (see m_metricCoords)
  m_metricCoords.assign(2*nbrSolPnts + 4*nbr1D, RealVector(2));
  m_metricPlaneIdx.resize(2*nbrSolPnts + 4*nbr1D);
  for (CFuint iKsi = 0; iKsi < nbr1D; ++iKsi)
  {
    for (CFuint iEta = 0; iEta < nbr1D; ++iEta)
    {
      const CFuint iSol = iKsi*nbr1D + iEta;
      m_metricCoords[iSol][KSI] = m_metricCoords[nbrSolPnts+iSol][KSI] = (*solPnts1D)[iKsi];
      m_metricCoords[iSol][ETA] = m_metricCoords[nbrSolPnts+iSol][ETA] = (*solPnts1D)[iEta];
      m_metricPlaneIdx[iSol] = KSI;
      m_metricPlaneIdx[nbrSolPnts+iSol] = ETA;
    }
  }
  for (CFuint i = 0; i < nbr1D; ++i)
  {
    const CFreal x = (*solPnts1D)[i];
    const CFuint start = 2*nbrSolPnts;
    m_metricCoords[start+i][KSI] = -1.0;         m_metricCoords[start+i][ETA] = x;
    m_metricCoords[start+nbr1D+i][KSI] = x;      m_metricCoords[start+nbr1D+i][ETA] = -1.0;
    m_metricCoords[start+2*nbr1D+i][KSI] = +1.0; m_metricCoords[start+2*nbr1D+i][ETA] = x;
    m_metricCoords[start+3*nbr1D+i][KSI] = x;    m_metricCoords[start+3*nbr1D+i][ETA] = +1.0;
    m_metricPlaneIdx[start+i] = m_metricPlaneIdx[start+2*nbr1D+i] = KSI;
    m_metricPlaneIdx[start+nbr1D+i] = m_metricPlaneIdx[start+3*nbr1D+i] = ETA;
  }

  // --- data for the linear reconstruction (see setReconstruction) ---

  m_solPnts1D.assign(solPnts1D->begin(), solPnts1D->end());
  m_bnds1D = bnds1D;

  // flux point ending each line. A flux point on a face normal to dir sits on the line
  // of its closest solution point, at the low end if that point is the first of the line.
  m_lineFlx.assign(2, vector< vector< CFuint > >(nbr1D, vector< CFuint >(2, 0)));
  const CFuint nbrFlxPnts = m_closestSolToFlx->size();
  vector< vector< vector< bool > > > found(2, vector< vector< bool > >(nbr1D, vector< bool >(2, false)));
  for (CFuint flxIdx = 0; flxIdx < nbrFlxPnts; ++flxIdx)
  {
    const CFuint sol = (*m_closestSolToFlx)[flxIdx];
    const CFuint dir = (*m_flxPntFlxDim)[flxIdx];
    const CFuint i = (dir == KSI) ? sol/nbr1D : sol%nbr1D;
    const CFuint j = (dir == KSI) ? sol%nbr1D : sol/nbr1D;
    // with one point per line the side follows from the face, not from i
    const CFuint side = (nbr1D == 1) ? (found[dir][j][0] ? 1 : 0) : (i == 0 ? 0 : 1);
    m_lineFlx[dir][j][side] = flxIdx;
    found[dir][j][side] = true;
  }
  for (CFuint dir = 0; dir < 2; ++dir)
  {
    for (CFuint j = 0; j < nbr1D; ++j)
    {
      cf_assert(found[dir][j][0] && found[dir][j][1]);
    }
  }

  m_slopes.assign(2, vector< RealVector >(nbrSolPnts, RealVector(m_nbrEqs)));
  m_reconstructed.assign(2, vector< bool >(nbrSolPnts, true));
  m_pertSlopes.assign(nbrSolPnts, RealVector(m_nbrEqs));
  m_testState.resize(m_nbrEqs);
  m_pointSlope.resize(m_nbrEqs);
  m_recPrev.resize(m_nbrEqs);
  m_recCur.resize(m_nbrEqs);
  m_recNext.resize(m_nbrEqs);
  m_recFace.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::setReconstruction(const std::string& limiter, const CFreal limiterEps,
                                               SafePtr< ConvectiveVarSet > updateVarSet)
{
  if      (limiter == "Minmod")    m_limiter = 0;
  else if (limiter == "VanAlbada") m_limiter = 1;
  else if (limiter == "None")      m_limiter = 2;
  else
  {
    throw BadValueException(FromHere(),
      "SubcellBlendingQuadData: SubcellLimiter must be Minmod, VanAlbada or None, not " + limiter);
  }

  if (m_nbrSolPnts1D < 2)
  {
    throw BadValueException(FromHere(),
      "SubcellBlendingQuadData: the linear reconstruction needs at least two solution points per direction.");
  }

  m_linear = true;
  m_limiterEps = limiterEps;
  m_updateVarSet = updateVarSet;

  if (m_recStateL == CFNULL)
  {
    RealVector dummyCoord(m_dim);
    dummyCoord = 0.0;
    m_recNodeL  = new Node(dummyCoord, false);
    m_recNodeR  = new Node(dummyCoord, false);
    m_recStateL = new State();
    m_recStateR = new State();
    m_recStateL->setSpaceCoordinates(m_recNodeL);
    m_recStateR->setSpaceCoordinates(m_recNodeR);
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::setReconstructionVars(SafePtr< VarSetTransformer > toRec,
                                                    SafePtr< VarSetTransformer > fromRec)
{
  m_toRec   = toRec;
  m_fromRec = fromRec;
  if (m_transInState == CFNULL)
  {
    RealVector dummyCoord(m_dim);
    dummyCoord = 0.0;
    m_transInState = new State();
    m_transInState->setSpaceCoordinates(new Node(dummyCoord, false));
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::toRecVars(const RealVector& update, RealVector& rec)
{
  if (m_toRec.isNull())
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) rec[iEq] = update[iEq];
    return;
  }
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) (*m_transInState)[iEq] = update[iEq];
  const State& out = *(m_toRec->transform(m_transInState));
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) rec[iEq] = out[iEq];
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::fromRecVars(const RealVector& rec, RealVector& update)
{
  if (m_fromRec.isNull())
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) update[iEq] = rec[iEq];
    return;
  }
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) (*m_transInState)[iEq] = rec[iEq];
  const State& out = *(m_fromRec->transform(m_transInState));
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) update[iEq] = out[iEq];
}

//////////////////////////////////////////////////////////////////////////////

CFreal SubcellBlendingQuadData::limitSlope(const CFreal a, const CFreal b, const CFreal eps2) const
{
  if (m_limiter == 2) return 0.5*(a + b);

  if (m_limiter == 0)
  {
    // minmod: zero at an extremum (secants of opposite sign), else the smaller one
    if (a*b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
  }

  // smooth van Albada, (max(ab, 0)(a + b) + eps^2 (a + b))/(a^2 + b^2 + 2 eps^2): for ab > 0
  // the usual ((b^2 + eps^2) a + (a^2 + eps^2) b)/(a^2 + b^2 + 2 eps^2). Secants of opposite
  // sign much larger than eps (an overshoot) give about zero; secants much smaller than eps
  // (a smooth extremum) give their average.
  const CFreal den = a*a + b*b + 2.0*eps2;
  return (den > 0.0) ? (std::max(a*b, 0.0) + eps2)*(a + b)/den : 0.0;
}

//////////////////////////////////////////////////////////////////////////////

bool SubcellBlendingQuadData::computePointSlope(const CFuint dir, const CFuint j, const CFuint i,
                                                const std::vector< State* >& states,
                                                const CFreal* lowSample, const CFreal* highSample,
                                                RealVector& slope)
{
  const CFuint n = m_nbrSolPnts1D;
  const CFreal x = m_solPnts1D[i];

  // previous and next neighbour along the line: a solution point, or the face sample.
  // All three in reconstruction variables.
  const bool lowEnd  = (i == 0);
  const bool highEnd = (i + 1 == n);
  const CFreal xPrev = lowEnd  ? -1.0 : m_solPnts1D[i-1];
  const CFreal xNext = highEnd ? +1.0 : m_solPnts1D[i+1];

  toRecVars(*(states[lineSol(dir, j, i)]), m_recCur);
  if (lowEnd)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) m_testState[iEq] = lowSample[iEq];
    toRecVars(m_testState, m_recPrev);
  }
  else
  {
    toRecVars(*(states[lineSol(dir, j, i-1)]), m_recPrev);
  }
  if (highEnd)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) m_testState[iEq] = highSample[iEq];
    toRecVars(m_testState, m_recNext);
  }
  else
  {
    toRecVars(*(states[lineSol(dir, j, i+1)]), m_recNext);
  }

  // frozen limiter factors of this point, if any (see setLimiterFreeze)
  const CFuint stateID = states[lineSol(dir, j, i)]->getLocalID();
  const CFuint frozenIdx = 2*stateID + dir;
  const bool frozen = (m_frozenValid != CFNULL) && (m_frozenValid[frozenIdx] == 1);
  const bool record = (m_frozenValid != CFNULL) && !frozen && m_recordPhi;

  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    const CFreal u = m_recCur[iEq];
    const CFreal uPrev = m_recPrev[iEq];
    const CFreal uNext = m_recNext[iEq];
    const CFreal a = (u - uPrev)/(x - xPrev);
    const CFreal b = (uNext - u)/(xNext - x);
    const CFreal avg = 0.5*(a + b);
    if (frozen)
    {
      slope[iEq] = m_frozenPhi[m_nbrEqs*frozenIdx + iEq]*avg;
      continue;
    }
    // smoothing size: a fraction of the size of the values, per unit reference length
    const CFreal scale = m_limiterEps*std::max(std::abs(u), std::max(std::abs(uPrev), std::abs(uNext)));
    slope[iEq] = limitSlope(a, b, scale*scale);
    if (record)
    {
      m_frozenPhi[m_nbrEqs*frozenIdx + iEq] = (avg != 0.0) ? slope[iEq]/avg : 0.0;
    }
  }
  if (record) m_frozenValid[frozenIdx] = 1;

  // admissibility of the two reconstructed subcell face states, in update variables
  const CFreal xFaces[2] = {m_bnds1D[i], m_bnds1D[i+1]};
  for (CFuint iFace = 0; iFace < 2; ++iFace)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_recFace[iEq] = m_recCur[iEq] + slope[iEq]*(xFaces[iFace] - x);
    }
    fromRecVars(m_recFace, m_testState);
    if (!m_updateVarSet->isValid(m_testState))
    {
      slope = 0.0;
      return false;
    }
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::computeLineSlopes(const CFuint dir, const CFuint j,
                                                const std::vector< State* >& states,
                                                std::vector< RealVector >& lineSlopes)
{
  const CFreal* lowSample  = m_cellSamples + m_nbrEqs*m_lineFlx[dir][j][0];
  const CFreal* highSample = m_cellSamples + m_nbrEqs*m_lineFlx[dir][j][1];
  for (CFuint i = 0; i < m_nbrSolPnts1D; ++i)
  {
    const CFuint sol = lineSol(dir, j, i);
    m_reconstructed[dir][sol] = computePointSlope(dir, j, i, states, lowSample, highSample, lineSlopes[sol]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::computeCellSlopes(const std::vector< State* >& states,
                                                const CFreal* cellSamples)
{
  m_cellSamples = cellSamples;
  for (CFuint dir = 0; dir < 2; ++dir)
  {
    for (CFuint j = 0; j < m_nbrSolPnts1D; ++j)
    {
      computeLineSlopes(dir, j, states, m_slopes[dir]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::reconstructAtElementFace(const std::vector< State* >& states,
                                                       const CFuint flxIdx, RealVector& sample,
                                                       State& result)
{
  const CFuint sol = (*m_closestSolToFlx)[flxIdx];
  const CFuint dir = (*m_flxPntFlxDim)[flxIdx];
  const CFuint n   = m_nbrSolPnts1D;
  const CFuint i   = (dir == KSI) ? sol/n : sol%n;
  const CFuint j   = (dir == KSI) ? sol%n : sol/n;
  const bool atLow = (m_lineFlx[dir][j][0] == flxIdx);

  // only the sample on this face is read: for n > 1 the point is an end point on this side
  const CFreal* samplePtr = sample.ptr();
  computePointSlope(dir, j, i, states, atLow ? samplePtr : CFNULL, atLow ? CFNULL : samplePtr,
                    m_pointSlope);

  const CFreal dx = (atLow ? -1.0 : +1.0) - m_solPnts1D[i];
  toRecVars(*(states[sol]), m_recCur);
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_recFace[iEq] = m_recCur[iEq] + m_pointSlope[iEq]*dx;
  }
  fromRecVars(m_recFace, result);
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
  if (m_telescope)
  {
    telescopeNormals(cell);
    return;
  }

  // plane normals sampled at the faces; the plane index is given per point, so both
  // reference directions are done in one call
  m_intfNormals = cell->computeMappedCoordPlaneNormalAtMappedCoords(m_intfPlaneIdx, m_intfCoords);
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::setCorrectionFunction(const std::vector< std::vector< CFreal > >& corrFctDiv)
{
  const CFuint n = m_nbrSolPnts1D;
  m_telescope = false;
  if (n < 2) return;

  // On a quad, the correction function of the flux point that ends the line (dir, j) at the
  // low (high) side is gL (gR) along the line times the Lagrange polynomial of eta_j (xi_j)
  // across it, so its divergence is +-gL'(x_i) on the points of that line and zero elsewhere.
  // The sign depends on the face orientation convention and is fixed with the quadrature of
  // the derivative: sum_i w_i gL'(x_i) = gL(+1) - gL(-1) = -1, and +1 for gR.
  const CFreal tol = 1.0e-8;
  m_corrDerivL.assign(n, 0.0);
  m_corrDerivR.assign(n, 0.0);
  bool ok = true;
  for (CFuint dir = 0; ok && dir < 2; ++dir)
  {
    for (CFuint j = 0; ok && j < n; ++j)
    {
      for (CFuint side = 0; ok && side < 2; ++side)
      {
        const CFuint flxIdx = m_lineFlx[dir][j][side];
        CFreal quad = 0.0;
        for (CFuint i = 0; i < n; ++i) quad += m_widths1D[i]*corrFctDiv[lineSol(dir, j, i)][flxIdx];
        const CFreal target = (side == 0) ? -1.0 : 1.0;
        if (std::abs(std::abs(quad) - 1.0) > tol) { ok = false; break; }
        const CFreal sign = target/quad;
        vector< CFreal >& deriv = (side == 0) ? m_corrDerivL : m_corrDerivR;
        for (CFuint iSol = 0; iSol < n*n; ++iSol)
        {
          // position of iSol on the line, or not on it
          const CFuint iLine = (dir == KSI) ? iSol/n : iSol%n;
          const bool onLine = (lineSol(dir, j, iLine) == iSol);
          const CFreal value = sign*corrFctDiv[iSol][flxIdx];
          if (!onLine)
          {
            if (std::abs(value) > tol) ok = false;
          }
          else if (dir == KSI && j == 0)
          {
            deriv[iLine] = value;
          }
          else if (std::abs(value - deriv[iLine]) > tol*(1.0 + std::abs(deriv[iLine])))
          {
            ok = false;
          }
        }
      }
    }
  }

  // this form also requires the flux points at the transverse solution point coordinates,
  // so that the face normals at the line ends are those of the FR face fluxes
  if (!ok)
  {
    CFLog(WARN, "SubcellBlendingQuadData: the correction function divergence does not have the "
          "tensor-product form expected on quads, the subcell normals are sampled at the faces\n");
    return;
  }
  m_telescope = true;
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::telescopeNormals(GeometricEntity* cell)
{
  // n_(i+1/2) = S(-1) + sum_(k<=i) w_k dS_k along each line, dS_k the FR derivative of the plane
  // normal at point k including the correction terms (see the class comment)
  const CFuint n = m_nbrSolPnts1D;
  const CFuint nbrSol = n*n;
  const vector< RealVector > metric =
    cell->computeMappedCoordPlaneNormalAtMappedCoords(m_metricPlaneIdx, m_metricCoords);

  for (CFuint dir = 0; dir < 2; ++dir)
  {
    // offsets of the plane normals of this direction in metric
    const CFuint solOffset = (dir == KSI) ? 0 : nbrSol;
    const CFuint lowOffset = 2*nbrSol + ((dir == KSI) ? 0 : n);
    const CFuint highOffset = 2*nbrSol + ((dir == KSI) ? 2*n : 3*n);

    for (CFuint j = 0; j < n; ++j)
    {
      // face minus extrapolated plane normal at both line ends
      m_endJumpL = metric[lowOffset + j];
      m_endJumpR = metric[highOffset + j];
      for (CFuint m = 0; m < n; ++m)
      {
        const RealVector& Sm = metric[solOffset + lineSol(dir, j, m)];
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {
          m_endJumpL[iDim] -= m_lagrangeAtEnds[0][m]*Sm[iDim];
          m_endJumpR[iDim] -= m_lagrangeAtEnds[1][m]*Sm[iDim];
        }
      }

      m_faceSum = metric[lowOffset + j];
      for (CFuint i = 0; i + 1 < n; ++i)
      {
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {
          m_metricDeriv[iDim] = m_corrDerivL[i]*m_endJumpL[iDim] + m_corrDerivR[i]*m_endJumpR[iDim];
        }
        for (CFuint m = 0; m < n; ++m)
        {
          const RealVector& Sm = metric[solOffset + lineSol(dir, j, m)];
          for (CFuint iDim = 0; iDim < m_dim; ++iDim)
          {
            m_metricDeriv[iDim] += m_derivMat1D[i][m]*Sm[iDim];
          }
        }
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {
          m_faceSum[iDim] += m_widths1D[i]*m_metricDeriv[iDim];
        }
        // interface between points i and i+1 of the line (layout of setup)
        const CFuint iIntf = (dir == KSI) ? i*n + j : (n-1)*n + j*(n-1) + i;
        m_intfNormals[iIntf] = m_faceSum;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::computeIntfFlux(const CFuint iIntf,
                                              const std::vector< State* >& states,
                                              const std::vector< RealVector >& slopes,
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

  const CFuint solL = m_intfSolL[iIntf];
  const CFuint solR = m_intfSolR[iIntf];

  if (!m_linear || !m_cellLinear)
  {
    const RealVector& flux = riemannFlux.computeFlux(*(states[solL]), *(states[solR]), m_unitNormal);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      result[iEq] = flux[iEq]*normalSize;
    }
    return;
  }

  // states reconstructed to the interface along its reference direction
  const CFuint dir = m_intfPlaneIdx[iIntf];
  const CFuint n   = m_nbrSolPnts1D;
  const CFuint iL  = (dir == KSI) ? solL/n : solL%n;
  const CFreal xF  = m_bnds1D[iL+1];
  const CFreal dxL = xF - m_solPnts1D[iL];
  const CFreal dxR = xF - m_solPnts1D[iL+1];
  toRecVars(*(states[solL]), m_recCur);
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_recFace[iEq] = m_recCur[iEq] + slopes[solL][iEq]*dxL;
  }
  fromRecVars(m_recFace, *m_recStateL);
  toRecVars(*(states[solR]), m_recCur);
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_recFace[iEq] = m_recCur[iEq] + slopes[solR][iEq]*dxR;
  }
  fromRecVars(m_recFace, *m_recStateR);

  const RealVector& flux = riemannFlux.computeFlux(*m_recStateL, *m_recStateR, m_unitNormal);

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
    const vector< RealVector >& slopes = m_slopes[m_intfPlaneIdx[iIntf]];
    computeIntfFlux(iIntf, states, slopes, riemannFlux, m_intfFlux[iIntf]);

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
  if (!m_linear || !m_cellLinear)
  {
    const std::vector< CFuint >& pertIntfs = m_intfOfSol[pertSol];
    for (CFuint k = 0; k < pertIntfs.size(); ++k)
    {
      addIntfFluxDelta(alpha, pertIntfs[k], states, m_slopes[0], riemannFlux, res);
    }
    return;
  }

  // the perturbation changes the slopes of the points of its two lines, so every
  // interface of those lines. The face samples stay frozen.
  const CFuint n = m_nbrSolPnts1D;
  for (CFuint dir = 0; dir < 2; ++dir)
  {
    const CFuint j = (dir == KSI) ? pertSol%n : pertSol/n;
    computeLineSlopes(dir, j, states, m_pertSlopes);
    for (CFuint i = 0; i + 1 < n; ++i)
    {
      const CFuint iIntf = (dir == KSI) ? i*n + j : (n-1)*n + j*(n-1) + i;
      addIntfFluxDelta(alpha, iIntf, states, m_pertSlopes, riemannFlux, res);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubcellBlendingQuadData::addIntfFluxDelta(const CFreal alpha, const CFuint iIntf,
                                               const std::vector< State* >& states,
                                               const std::vector< RealVector >& slopes,
                                               RiemannFlux& riemannFlux,
                                               RealVector& res)
{
  computeIntfFlux(iIntf, states, slopes, riemannFlux, m_intfFluxPert[0]);

  const RealVector& fluxNew = m_intfFluxPert[0];
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

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
