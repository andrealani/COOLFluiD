// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstructionSubcellBlending.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvBndCorrectionsRHSFluxReconstructionSubcellBlending, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvBndCorrectionsRHSFluxReconstructionSubcellBlendingProvider("ConvBndCorrectionsRHSSubcellBlending");

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSFluxReconstructionSubcellBlending::ConvBndCorrectionsRHSFluxReconstructionSubcellBlending(const std::string& name) :
  ConvBndCorrectionsRHSFluxReconstruction(name),
  socket_alpha("alpha"),
  m_scData()
{
}

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSFluxReconstructionSubcellBlending::~ConvBndCorrectionsRHSFluxReconstructionSubcellBlending()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
ConvBndCorrectionsRHSFluxReconstructionSubcellBlending::needsSockets()
{
  std::vector< SafePtr< BaseDataSocketSink > > result = ConvBndCorrectionsRHSFluxReconstruction::needsSockets();
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvBndCorrectionsRHSFluxReconstructionSubcellBlending::getCellAlpha()
{
  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  return alpha[(*m_cellStates)[0]->getLocalID()];
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionSubcellBlending::computeCorrection(vector< RealVector >& corrections)
{
  // FR correction -(F* divh) of all solution points
  ConvBndCorrectionsRHSFluxReconstruction::computeCorrection(corrections);

  const CFreal alpha = getCellAlpha();
  if (alpha <= 0.0) return;

  const CFreal oneMinusAlpha = 1.0 - alpha;
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    corrections[iSolPnt] *= oneMinusAlpha;
  }

  // each face flux point is the boundary of the subcell of its closest solution point.
  // m_flxPntRiemannFlux holds the boundary flux times the face Jacobian. For boundary
  // faces m_orient is the local face index of the cell.
  const CFint faceDir = (*m_faceMappedCoordDir)[m_orient];

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    const CFuint solIdx = m_scData.getClosestSol(flxIdx);

    const CFreal factor = alpha*static_cast<CFreal>(faceDir)/m_scData.getFaceSubcellWidth(flxIdx);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      corrections[solIdx][iEq] -= factor*m_flxPntRiemannFlux[iFlxPnt][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionSubcellBlending::setup()
{
  CFAUTOTRACE;

  ConvBndCorrectionsRHSFluxReconstruction::setup();

  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() == 1);

  m_scData.setup(frLocalData[0], m_dim, m_nbrEqs);
  cf_assert(m_scData.getNbrSolPnts1D()*m_scData.getNbrSolPnts1D() == m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionSubcellBlending::unsetup()
{
  CFAUTOTRACE;

  ConvBndCorrectionsRHSFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
