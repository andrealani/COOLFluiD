// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/ConvRHSJacobFluxReconstructionSubcellBlending.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvRHSJacobFluxReconstructionSubcellBlending, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvRHSJacobFluxReconstructionSubcellBlendingProvider("ConvRHSJacobSubcellBlending");

//////////////////////////////////////////////////////////////////////////////

ConvRHSJacobFluxReconstructionSubcellBlending::ConvRHSJacobFluxReconstructionSubcellBlending(const std::string& name) :
  ConvRHSJacobFluxReconstruction(name),
  socket_alpha("alpha"),
  m_scData(),
  m_currFaceAlphaF(0.0),
  m_currCellAlpha(0.0),
  m_flxPntFaceConn(CFNULL),
  m_blendedFaceFlux()
{
  addConfigOptionsTo(this);

  m_faceFluxBlending = true;
  setParameter("FaceFluxBlending", &m_faceFluxBlending);
}

//////////////////////////////////////////////////////////////////////////////

ConvRHSJacobFluxReconstructionSubcellBlending::~ConvRHSJacobFluxReconstructionSubcellBlending()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionSubcellBlending::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("FaceFluxBlending",
    "Blend the face Riemann flux with the first-order flux between the adjacent solution points, "
    "weighted by max(alpha_L, alpha_R). With alpha = 1 the cell then runs a pure subcell P0 scheme.");
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionSubcellBlending::configure(Config::ConfigArgs& args)
{
  ConvRHSJacobFluxReconstruction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
ConvRHSJacobFluxReconstructionSubcellBlending::needsSockets()
{
  std::vector< SafePtr< BaseDataSocketSink > > result = ConvRHSJacobFluxReconstruction::needsSockets();
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvRHSJacobFluxReconstructionSubcellBlending::getCellAlpha(const std::vector< State* >& states)
{
  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  return alpha[states[0]->getLocalID()];
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvRHSJacobFluxReconstructionSubcellBlending::computeFaceAlphaF()
{
  if (!m_faceFluxBlending) return 0.0;

  // both neighbours evaluate the same expression, so the face flux stays single-valued
  return std::max(getCellAlpha(*m_states[LEFT]), getCellAlpha(*m_states[RIGHT]));
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionSubcellBlending::computeInterfaceFlxCorrection()
{
  // face Riemann flux from the extrapolated states, unscaled in m_flxPntRiemannFlux
  // and scaled by the face Jacobian in m_cellFlx
  ConvRHSJacobFluxReconstruction::computeInterfaceFlxCorrection();

  m_currFaceAlphaF = computeFaceAlphaF();
  if (m_currFaceAlphaF <= 0.0) return;

  const CFreal oneMinusAlphaF = 1.0 - m_currFaceAlphaF;

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT ][iFlxPnt];
    const CFuint flxIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];
    State& solStateL = *((*m_states[LEFT ])[m_scData.getClosestSol(flxIdxL)]);
    State& solStateR = *((*m_states[RIGHT])[m_scData.getClosestSol(flxIdxR)]);

    // first-order flux between the adjacent solution point states
    const RealVector& loFlux = m_riemannFluxComputer->computeFlux(solStateL, solStateR,
                                                                 m_unitNormalFlxPnts[iFlxPnt]);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_flxPntRiemannFlux[iFlxPnt][iEq] = oneMinusAlphaF*m_flxPntRiemannFlux[iFlxPnt][iEq]
                                        + m_currFaceAlphaF*loFlux[iEq];
    }

    m_cellFlx[LEFT ][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT ];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionSubcellBlending::computeCorrection(CFuint side, vector< RealVector >& corrections)
{
  // FR correction -(F* divh) of all solution points
  ConvRHSJacobFluxReconstruction::computeCorrection(side, corrections);

  const CFreal alpha = getCellAlpha(*m_states[side]);
  if (alpha <= 0.0) return;

  const CFreal oneMinusAlpha = 1.0 - alpha;
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    corrections[iSolPnt] *= oneMinusAlpha;
  }

  // each face flux point is the boundary of the subcell of its closest solution point
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx  = (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt];
    const CFuint solIdx  = m_scData.getClosestSol(flxIdx);
    const CFint  faceDir = (*m_faceLocalDir)[(*m_flxPntFaceConn)[flxIdx]];

    const CFreal factor = alpha*static_cast<CFreal>(faceDir)/m_scData.getFaceSubcellWidth(flxIdx);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      corrections[solIdx][iEq] -= factor*m_cellFlx[side][iFlxPnt][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionSubcellBlending::computeDivDiscontFlx(vector< RealVector >& residuals)
{
  // FR volume term -div(F_D) + div(h F_D) of all solution points
  ConvRHSJacobFluxReconstruction::computeDivDiscontFlx(residuals);

  const CFreal alpha = m_currCellAlpha;
  if (alpha <= 0.0) return;

  const CFreal oneMinusAlpha = 1.0 - alpha;
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    residuals[iSolPnt] *= oneMinusAlpha;
  }

  // fills the subcell residual buffer and stores the unperturbed interface fluxes,
  // both reused by computePertDivDiscontFlx for this cell
  m_scData.computeSubcellRes(alpha, *m_cellStates, *m_riemannFluxComputer);
  const vector< RealVector >& subcellRes = m_scData.getSubcellRes();

  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      residuals[iSolPnt][iEq] += subcellRes[iSolPnt][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionSubcellBlending::setCellData()
{
  // flux projection vectors at the solution points
  ConvRHSJacobFluxReconstruction::setCellData();

  m_currCellAlpha = getCellAlpha(*m_cellStates);
  if (m_currCellAlpha <= 0.0) return;

  m_scData.computeCellNormals(m_cell);
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionSubcellBlending::computePertInterfaceFlxCorrection()
{
  // the parent stores only the flux scaled by the face Jacobian, so the blend is
  // applied here on the unscaled flux before scaling
  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    const CFuint iFlx = m_influencedFlxPnts[iFlxPnt];

    m_blendedFaceFlux = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT ][iFlx]),
                                                          *(m_cellStatesFlxPnt[RIGHT][iFlx]),
                                                          m_unitNormalFlxPnts[iFlx]);

    if (m_currFaceAlphaF > 0.0)
    {
      const CFuint flxIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT ][iFlx];
      const CFuint flxIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx];

      // the closest solution point states carry the perturbation directly
      const RealVector& loFlux = m_riemannFluxComputer->computeFlux(
          *((*m_states[LEFT ])[m_scData.getClosestSol(flxIdxL)]),
          *((*m_states[RIGHT])[m_scData.getClosestSol(flxIdxR)]),
          m_unitNormalFlxPnts[iFlx]);

      const CFreal oneMinusAlphaF = 1.0 - m_currFaceAlphaF;
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_blendedFaceFlux[iEq] = oneMinusAlphaF*m_blendedFaceFlux[iEq] + m_currFaceAlphaF*loFlux[iEq];
      }
    }

    m_cellFlx[LEFT ][iFlx] = m_blendedFaceFlux*m_faceJacobVecSizeFlxPnts[iFlx][LEFT ];
    m_cellFlx[RIGHT][iFlx] = m_blendedFaceFlux*m_faceJacobVecSizeFlxPnts[iFlx][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionSubcellBlending::computePertCorrection(CFuint side, RealVector& corrections)
{
  // the parent resets the whole vector and fills only the influenced flux points,
  // so scaling the whole vector leaves the untouched entries at zero
  ConvRHSJacobFluxReconstruction::computePertCorrection(side, corrections);

  const CFreal alpha = getCellAlpha(*m_states[side]);
  if (alpha <= 0.0) return;

  corrections *= (1.0 - alpha);

  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    const CFuint iFlx    = m_influencedFlxPnts[iFlxPnt];
    const CFuint flxIdx  = (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlx];
    const CFuint solIdx  = m_scData.getClosestSol(flxIdx);
    const CFint  faceDir = (*m_faceLocalDir)[(*m_flxPntFaceConn)[flxIdx]];

    const CFreal factor = alpha*static_cast<CFreal>(faceDir)/m_scData.getFaceSubcellWidth(flxIdx);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      corrections[m_nbrEqs*solIdx+iEq] -= factor*m_cellFlx[side][iFlx][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionSubcellBlending::computePertDivDiscontFlx(RealVector& residuals)
{
  // the parent fills only the solution points in the dependency list of the perturbed one
  ConvRHSJacobFluxReconstruction::computePertDivDiscontFlx(residuals);

  const CFreal alpha = m_currCellAlpha;
  if (alpha <= 0.0) return;

  const CFreal oneMinusAlpha = 1.0 - alpha;
  const vector< RealVector >& subcellRes = m_scData.getSubcellRes();

  // scale the FR part and add back the unperturbed subcell contribution, on the same
  // solution points the parent filled
  const CFuint nbrSolSolDep = ((*m_solSolDep)[m_pertSol]).size();
  for (CFuint iSolPnt = 0; iSolPnt < nbrSolSolDep; ++iSolPnt)
  {
    const CFuint iSolIdx = (*m_solSolDep)[m_pertSol][iSolPnt];

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      residuals[m_nbrEqs*iSolIdx+iEq] = oneMinusAlpha*residuals[m_nbrEqs*iSolIdx+iEq]
                                      + subcellRes[iSolIdx][iEq];
    }
  }

  // correct the interfaces that touch the perturbed solution point. Both of their
  // endpoints share a reference line with it, so this stays inside the region above.
  m_scData.addSubcellResDelta(alpha, m_pertSol, *m_cellStates, *m_riemannFluxComputer, residuals);
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionSubcellBlending::setup()
{
  CFAUTOTRACE;

  ConvRHSJacobFluxReconstruction::setup();

  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() == 1);

  m_scData.setup(frLocalData[0], m_dim, m_nbrEqs);
  cf_assert(m_scData.getNbrSolPnts1D()*m_scData.getNbrSolPnts1D() == m_nbrSolPnts);

  m_flxPntFaceConn = frLocalData[0]->getFlxPntFaceConn();

  m_blendedFaceFlux.resize(m_nbrEqs);

  CFLog(INFO, "ConvRHSJacobFluxReconstructionSubcellBlending: FaceFluxBlending = " << m_faceFluxBlending
        << ", 1D subcell widths =");
  const vector< CFreal >& widths = m_scData.getWidths1D();
  for (CFuint i = 0; i < widths.size(); ++i)
  {
    CFLog(INFO, " " << widths[i]);
  }
  CFLog(INFO, "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionSubcellBlending::unsetup()
{
  CFAUTOTRACE;

  ConvRHSJacobFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
