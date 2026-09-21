// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/ConvRHSFluxReconstructionSubcellBlending.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvRHSFluxReconstructionSubcellBlending, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvRHSFluxReconstructionSubcellBlendingProvider("ConvRHSSubcellBlending");

//////////////////////////////////////////////////////////////////////////////

ConvRHSFluxReconstructionSubcellBlending::ConvRHSFluxReconstructionSubcellBlending(const std::string& name) :
  ConvRHSFluxReconstruction(name),
  socket_alpha("alpha"),
  m_scData(),
  m_currFaceAlphaF(0.0),
  m_currCellAlpha(0.0),
  m_flxPntFaceConn(CFNULL)
{
  addConfigOptionsTo(this);

  m_faceFluxBlending = true;
  setParameter("FaceFluxBlending", &m_faceFluxBlending);
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("FaceFluxBlending",
    "Blend the face Riemann flux with the first-order flux between the adjacent solution points, "
    "weighted by max(alpha_L, alpha_R). With alpha = 1 the cell then runs a pure subcell P0 scheme.");
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::configure(Config::ConfigArgs& args)
{
  ConvRHSFluxReconstruction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
ConvRHSFluxReconstructionSubcellBlending::needsSockets()
{
  std::vector< SafePtr< BaseDataSocketSink > > result = ConvRHSFluxReconstruction::needsSockets();
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvRHSFluxReconstructionSubcellBlending::getCellAlpha(const std::vector< State* >& states)
{
  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  return alpha[states[0]->getLocalID()];
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvRHSFluxReconstructionSubcellBlending::computeFaceAlphaF()
{
  if (!m_faceFluxBlending) return 0.0;

  // both neighbours evaluate the same expression, so the face flux stays single-valued
  return std::max(getCellAlpha(*m_states[LEFT]), getCellAlpha(*m_states[RIGHT]));
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::computeInterfaceFlxCorrection()
{
  // face Riemann flux from the extrapolated states, unscaled in m_flxPntRiemannFlux
  // and scaled by the face Jacobian in m_cellFlx
  ConvRHSFluxReconstruction::computeInterfaceFlxCorrection();

  m_currFaceAlphaF = computeFaceAlphaF();
  if (m_currFaceAlphaF <= 0.0) return;

  const CFreal oneMinusAlphaF = 1.0 - m_currFaceAlphaF;

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // solution points adjacent to this flux point in the left and right cell
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

    // interface flux in the mapped coordinate frame of each neighbour
    m_cellFlx[LEFT ][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT ];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::addSubcellFaceFlux(const CFuint side, const CFreal alpha,
                                                                  vector< RealVector >& corrections)
{
  // each face flux point is the boundary of the subcell of its closest solution point.
  // m_cellFlx holds the face flux times the face Jacobian, oriented along the positive
  // mapped direction of the face normal.
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

void ConvRHSFluxReconstructionSubcellBlending::computeCorrection(CFuint side, vector< RealVector >& corrections)
{
  // FR correction -(F* divh) of all solution points
  ConvRHSFluxReconstruction::computeCorrection(side, corrections);

  const CFreal alpha = getCellAlpha(*m_states[side]);
  if (alpha <= 0.0) return;

  const CFreal oneMinusAlpha = 1.0 - alpha;
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    corrections[iSolPnt] *= oneMinusAlpha;
  }

  addSubcellFaceFlux(side, alpha, corrections);
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::computeDivDiscontFlx(vector< RealVector >& residuals)
{
  // FR volume term -div(F_D) + div(h F_D) of all solution points
  ConvRHSFluxReconstruction::computeDivDiscontFlx(residuals);

  const CFreal alpha = m_currCellAlpha;
  if (alpha <= 0.0) return;

  const CFreal oneMinusAlpha = 1.0 - alpha;
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    residuals[iSolPnt] *= oneMinusAlpha;
  }

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

void ConvRHSFluxReconstructionSubcellBlending::setCellData()
{
  // flux projection vectors at the solution points
  ConvRHSFluxReconstruction::setCellData();

  m_currCellAlpha = getCellAlpha(*m_cellStates);
  if (m_currCellAlpha <= 0.0) return;

  m_scData.computeCellNormals(m_cell);
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::setup()
{
  CFAUTOTRACE;

  ConvRHSFluxReconstruction::setup();

  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() == 1);

  m_scData.setup(frLocalData[0], m_dim, m_nbrEqs);
  cf_assert(m_scData.getNbrSolPnts1D()*m_scData.getNbrSolPnts1D() == m_nbrSolPnts);

  m_flxPntFaceConn = frLocalData[0]->getFlxPntFaceConn();

  CFLog(INFO, "ConvRHSFluxReconstructionSubcellBlending: FaceFluxBlending = " << m_faceFluxBlending
        << ", 1D subcell widths =");
  const vector< CFreal >& widths = m_scData.getWidths1D();
  for (CFuint i = 0; i < widths.size(); ++i)
  {
    CFLog(INFO, " " << widths[i]);
  }
  CFLog(INFO, "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::unsetup()
{
  CFAUTOTRACE;

  ConvRHSFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
