// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "Framework/PhysicalModel.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstructionBlending.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvBndCorrectionsRHSFluxReconstructionBlending, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvBndCorrectionsRHSFluxReconstructionBlendingProvider("ConvBndCorrectionsRHSBlending");

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSFluxReconstructionBlending::ConvBndCorrectionsRHSFluxReconstructionBlending
(const std::string& name) :
  ConvBndCorrectionsRHSFluxReconstruction(name),
  socket_alpha("alpha"),
  socket_volumes("volumes"),
  m_p0Data(),
  m_cellAvgState(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
ConvBndCorrectionsRHSFluxReconstructionBlending::needsSockets()
{
  std::vector< SafePtr< BaseDataSocketSink > > result = ConvBndCorrectionsRHSFluxReconstruction::needsSockets();
  result.push_back(&socket_alpha);
  result.push_back(&socket_volumes);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvBndCorrectionsRHSFluxReconstructionBlending::getCellAlpha()
{
  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  return alpha[(*m_cellStates)[0]->getLocalID()];
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::computeCorrection(vector< RealVector >& corrections)
{
  // FR correction -(Fhat divh)
  ConvBndCorrectionsRHSFluxReconstruction::computeCorrection(corrections);

  const CFreal alpha = getCellAlpha();
  if (alpha <= 0.0) return;

  const CFreal oneMinusAlpha = 1.0 - alpha;
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    corrections[iSolPnt] *= oneMinusAlpha;
  }

  // -alpha sum_fp W_fp Fhat_fp, where m_flxPntRiemannFlux holds Fhat_fp times the signed face Jacobian
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFreal factor = alpha*m_p0Data.getW((*m_faceFlxPntConn)[m_orient][iFlxPnt]);

    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        corrections[iSolPnt][iEq] -= factor*m_flxPntRiemannFlux[iFlxPnt][iEq];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvBndCorrectionsRHSFluxReconstructionBlending::computeWaveSpeedP0()
{
  DataHandle< CFreal > volumes = socket_volumes.getDataHandle();
  m_p0Data.computeCellAverage(*m_cellStates, volumes, *m_cellAvgState);

  // face integral of the largest eigenvalue at the cell average; m_faceIntegrationCoefs
  // was set for this face type by computeWaveSpeedUpdates
  m_updateVarSet->computePhysicalData(*m_cellAvgState, m_pDataP0);

  CFreal waveSpeed = 0.0;
  CFreal sumW = 0.0;
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFreal jacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceIntegrationCoefs)[iFlx];
    waveSpeed += jacobXIntCoef*m_updateVarSet->getMaxAbsEigenValue(m_pDataP0, m_unitNormalFlxPnts[iFlx]);
    sumW += m_p0Data.getW((*m_faceFlxPntConn)[m_orient][iFlx]);
  }

  // face weight in the cell average relative to the reference value 0.5, floored at 1:
  // faces that drive the cell average more get a larger P0 wave speed
  const CFreal factor = std::max(1.0, std::abs(sumW)/0.5);

  return waveSpeed*factor*m_nbrSolPnts;
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::updateWaveSpeed()
{
  DataHandle< CFreal > updateCoeff = socket_updateCoeff.getDataHandle();

  const CFreal alpha = getCellAlpha();

  if (alpha > 0.0)
  {
    const CFreal waveSpeedFR = m_waveSpeedUpd*(2.0*m_order+1);
    const CFreal waveSpeed = (1.0 - alpha)*waveSpeedFR + alpha*computeWaveSpeedP0();

    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      updateCoeff[(*m_cellStates)[iSol]->getLocalID()] += waveSpeed;
    }
  }
  else
  {
    // same statement as the base class, so that an unblended cell gets the plain FR value to the last bit
    // (the compiler may fuse the multiply and the add of one statement)
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      updateCoeff[(*m_cellStates)[iSol]->getLocalID()] += m_waveSpeedUpd*(2.0*m_order+1);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::setup()
{
  CFAUTOTRACE;

  ConvBndCorrectionsRHSFluxReconstruction::setup();

  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() == 1);

  SafePtr< PhysicalModel > physModel = PhysicalModelStack::getActive();
  m_p0Data.setup(frLocalData[0], m_corrFctDiv, m_flxSolDep, m_nbrEqs,
                 getMethodData().getUpdateToSolutionVecTrans(), getMethodData().getSolutionVar(),
                 getMethodData().getUpdateVar());

  RealVector dummyCoord(m_dim);
  dummyCoord = 0.0;
  // physical data of the cell average
  physModel->getImplementor()->getConvectiveTerm()->resizePhysicalData(m_pDataP0);

  m_cellAvgState = new State();
  m_cellAvgState->setSpaceCoordinates(new Node(dummyCoord, false));
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::unsetup()
{
  CFAUTOTRACE;

  // the state owns its coordinates (Node built with isOnMesh = false)
  deletePtr(m_cellAvgState);
  m_p0Data.unsetup();

  ConvBndCorrectionsRHSFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
