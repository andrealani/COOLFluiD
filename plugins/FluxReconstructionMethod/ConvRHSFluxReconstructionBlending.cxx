// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/ConvRHSFluxReconstructionBlending.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvRHSFluxReconstructionBlending, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvRHSFluxReconstructionBlendingProvider("ConvRHSBlending");

//////////////////////////////////////////////////////////////////////////////

ConvRHSFluxReconstructionBlending::ConvRHSFluxReconstructionBlending(const std::string& name) :
  ConvRHSFluxReconstruction(name),
  socket_alpha("alpha"),
  m_p0Data(),
  m_cellAvgStates(),
  m_faceAvgsDone(false),
  m_currFaceAlphaF(0.0),
  m_currCellAlpha(0.0),
  m_flxPntFluxP0()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
ConvRHSFluxReconstructionBlending::needsSockets()
{
  std::vector< SafePtr< BaseDataSocketSink > > result = ConvRHSFluxReconstruction::needsSockets();
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvRHSFluxReconstructionBlending::getCellAlpha(const std::vector< State* >& states)
{
  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  return alpha[states[0]->getLocalID()];
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionBlending::computeFaceCellAverages()
{
  if (m_faceAvgsDone) return;

  DataHandle< CFreal > volumes = socket_volumes.getDataHandle();
  m_p0Data.computeCellAverage(*m_states[LEFT ], volumes, *m_cellAvgStates[LEFT ]);
  m_p0Data.computeCellAverage(*m_states[RIGHT], volumes, *m_cellAvgStates[RIGHT]);
  m_faceAvgsDone = true;
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionBlending::computeInterfaceFlxCorrection()
{
  // FR interface flux F_I from the extrapolated states, unscaled in m_flxPntRiemannFlux
  // and scaled by the face Jacobian in m_cellFlx, which keeps it for the FR part
  ConvRHSFluxReconstruction::computeInterfaceFlxCorrection();

  m_faceAvgsDone = false;

  // both neighbours evaluate the same expression, so the face flux stays single-valued
  m_currFaceAlphaF = std::max(getCellAlpha(*m_states[LEFT]), getCellAlpha(*m_states[RIGHT]));
  if (m_currFaceAlphaF <= 0.0) return;

  computeFaceCellAverages();

  // P0 flux F* of the cell averages with the normal of each flux point, used by addMeanFaceTerm
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    m_flxPntFluxP0[iFlxPnt] = m_riemannFluxComputer->computeFlux(*m_cellAvgStates[LEFT],
                                                                 *m_cellAvgStates[RIGHT],
                                                                 m_unitNormalFlxPnts[iFlxPnt]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionBlending::addMeanFaceTerm(const CFuint side, const CFreal alpha,
                                                        vector< RealVector >& corrections)
{
  // -sum_fp W_fp [alpha_f F* + (alpha_K - alpha_f) F_I], where m_cellFlx holds F_I times the
  // signed face Jacobian of this side; with the FR part (1 - alpha_K) of F_I the cell average
  // moves by -sum_fp W_fp [(1 - alpha_f) F_I + alpha_f F*]
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt];
    const CFreal w = m_p0Data.getW(flxIdx);
    const CFreal factorP0 = w*m_currFaceAlphaF*m_faceJacobVecSizeFlxPnts[iFlxPnt][side];
    const CFreal factorFR = w*(alpha - m_currFaceAlphaF);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFreal term = factorP0*m_flxPntFluxP0[iFlxPnt][iEq] + factorFR*m_cellFlx[side][iFlxPnt][iEq];
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
      {
        corrections[iSolPnt][iEq] -= term;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionBlending::computeCorrection(CFuint side, vector< RealVector >& corrections)
{
  // FR correction -(F_I divh)
  ConvRHSFluxReconstruction::computeCorrection(side, corrections);

  const CFreal alpha = getCellAlpha(*m_states[side]);

  if (alpha > 0.0)
  {
    const CFreal oneMinusAlpha = 1.0 - alpha;
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      corrections[iSolPnt] *= oneMinusAlpha;
    }
  }

  // a cell with alpha_K = 0 next to a blended cell also needs the mean correction
  if (m_currFaceAlphaF > 0.0) addMeanFaceTerm(side, alpha, corrections);
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionBlending::computeDivDiscontFlx(vector< RealVector >& residuals)
{
  // FR volume term -div(F_D) + div(h F_D); the P0 part has no volume term
  ConvRHSFluxReconstruction::computeDivDiscontFlx(residuals);

  const CFreal alpha = m_currCellAlpha;
  if (alpha <= 0.0) return;

  const CFreal oneMinusAlpha = 1.0 - alpha;
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    residuals[iSolPnt] *= oneMinusAlpha;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionBlending::setCellData()
{
  ConvRHSFluxReconstruction::setCellData();

  m_currCellAlpha = getCellAlpha(*m_cellStates);
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvRHSFluxReconstructionBlending::computeWaveSpeedP0(const CFuint side)
{
  computeFaceCellAverages();

  // face integral of the largest eigenvalue at the cell average; m_faceIntegrationCoefs
  // was set for this face type by computeWaveSpeedUpdates
  m_updateVarSet->computePhysicalData(*m_cellAvgStates[side], m_pData);

  CFreal waveSpeed = 0.0;
  CFreal sumW = 0.0;
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFreal jacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceIntegrationCoefs)[iFlx];
    waveSpeed += jacobXIntCoef*m_updateVarSet->getMaxAbsEigenValue(m_pData, m_unitNormalFlxPnts[iFlx]);
    sumW += m_p0Data.getW((*m_faceFlxPntConnPerOrient)[m_orient][side][iFlx]);
  }

  // face weight in the cell average relative to the reference value 0.5, floored at 1:
  // faces that drive the cell average more get a larger P0 wave speed
  const CFreal factor = std::max(1.0, std::abs(sumW)/0.5);

  return waveSpeed*factor*m_nbrSolPnts;
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionBlending::updateWaveSpeed()
{
  DataHandle< CFreal > updateCoeff = socket_updateCoeff.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFreal alpha = getCellAlpha(*m_states[iSide]);

    if (alpha > 0.0)
    {
      const CFreal waveSpeedFR = m_waveSpeedUpd[iSide]*(2.0*m_order+1);
      const CFreal waveSpeed = (1.0 - alpha)*waveSpeedFR + alpha*computeWaveSpeedP0(iSide);

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        updateCoeff[(*m_states[iSide])[iSol]->getLocalID()] += waveSpeed;
      }
    }
    else
    {
      // same statement as the base class, so that an unblended cell gets the plain FR value to the last bit
      // (the compiler may fuse the multiply and the add of one statement)
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        updateCoeff[(*m_states[iSide])[iSol]->getLocalID()] += m_waveSpeedUpd[iSide]*(2.0*m_order+1);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionBlending::setup()
{
  CFAUTOTRACE;

  ConvRHSFluxReconstruction::setup();

  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() == 1);

  m_p0Data.setup(frLocalData[0], m_corrFctDiv, m_flxSolDep, m_nbrEqs,
                 getMethodData().getUpdateToSolutionVecTrans(), getMethodData().getSolutionVar(),
                 getMethodData().getUpdateVar());

  RealVector dummyCoord(m_dim);
  dummyCoord = 0.0;
  for (CFuint side = 0; side < 2; ++side)
  {
    m_cellAvgStates.push_back(new State());
    m_cellAvgStates[side]->setSpaceCoordinates(new Node(dummyCoord, false));
  }

  m_flxPntFluxP0.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    m_flxPntFluxP0[iFlxPnt].resize(m_nbrEqs);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionBlending::unsetup()
{
  CFAUTOTRACE;

  for (CFuint side = 0; side < m_cellAvgStates.size(); ++side)
  {
    // the state owns its coordinates (Node built with isOnMesh = false)
    deletePtr(m_cellAvgStates[side]);
  }
  m_cellAvgStates.clear();
  m_p0Data.unsetup();

  ConvRHSFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
