// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlendingProvider("ConvBndCorrectionsRHSJacobSubcellBlending");

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending::ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending(const std::string& name) :
  ConvBndCorrectionsRHSJacobFluxReconstruction(name),
  socket_alpha("alpha"),
  m_scData(),
  m_scSolScaled()
{
}

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending::~ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending::needsSockets()
{
  std::vector< SafePtr< BaseDataSocketSink > > result = ConvBndCorrectionsRHSJacobFluxReconstruction::needsSockets();
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending::getCellAlpha()
{
  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  return alpha[(*m_cellStates)[0]->getLocalID()];
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending::addSubcellFaceFlux(const CFuint iFlxPnt,
                                                                                     const CFreal alpha,
                                                                                     const CFint faceDir,
                                                                                     vector< RealVector >& corrections)
{
  // the face flux point is the boundary of the subcell of its closest solution point.
  // m_flxPntRiemannFlux holds the boundary flux times the face Jacobian.
  const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
  const CFuint solIdx = m_scData.getClosestSol(flxIdx);

  const CFreal factor = alpha*static_cast<CFreal>(faceDir)/m_scData.getFaceSubcellWidth(flxIdx);

  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    corrections[solIdx][iEq] -= factor*m_flxPntRiemannFlux[iFlxPnt][iEq];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending::computeCorrection(vector< RealVector >& corrections)
{
  // FR correction -(F* divh) of all solution points
  ConvBndCorrectionsRHSJacobFluxReconstruction::computeCorrection(corrections);

  const CFreal alpha = getCellAlpha();
  if (alpha <= 0.0) return;

  const CFreal oneMinusAlpha = 1.0 - alpha;
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    corrections[iSolPnt] *= oneMinusAlpha;
  }

  // for boundary faces m_orient is the local face index of the cell
  const CFint faceDir = (*m_faceMappedCoordDir)[m_orient];

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    addSubcellFaceFlux(iFlxPnt, alpha, faceDir, corrections);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending::computePertCorrection(vector< RealVector >& corrections)
{
  // the parent resets and refills only the solution points that depend on an
  // influenced flux point, leaving the others at their unperturbed value
  ConvBndCorrectionsRHSJacobFluxReconstruction::computePertCorrection(corrections);

  const CFreal alpha = getCellAlpha();
  if (alpha <= 0.0) return;

  const CFreal oneMinusAlpha = 1.0 - alpha;

  // scale exactly the solution points the parent refilled, each one once
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    m_scSolScaled[iSolPnt] = false;
  }

  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][ m_influencedFlxPnts[iFlxPnt] ];

    m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSol];

      if (!m_scSolScaled[solIdx])
      {
        corrections[solIdx] *= oneMinusAlpha;
        m_scSolScaled[solIdx] = true;
      }
    }
  }

  // for boundary faces m_orient is the local face index of the cell
  const CFint faceDir = (*m_faceMappedCoordDir)[m_orient];

  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    addSubcellFaceFlux(m_influencedFlxPnts[iFlxPnt], alpha, faceDir, corrections);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending::setup()
{
  CFAUTOTRACE;

  ConvBndCorrectionsRHSJacobFluxReconstruction::setup();

  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() == 1);

  m_scData.setup(frLocalData[0], m_dim, m_nbrEqs);
  cf_assert(m_scData.getNbrSolPnts1D()*m_scData.getNbrSolPnts1D() == m_nbrSolPnts);

  m_scSolScaled.resize(m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionSubcellBlending::unsetup()
{
  CFAUTOTRACE;

  ConvBndCorrectionsRHSJacobFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
