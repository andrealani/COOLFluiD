// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/ConvRHSJacobFluxReconstructionBlending.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvRHSJacobFluxReconstructionBlending, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvRHSJacobFluxReconstructionBlendingProvider("ConvRHSJacobBlending");

//////////////////////////////////////////////////////////////////////////////

ConvRHSJacobFluxReconstructionBlending::ConvRHSJacobFluxReconstructionBlending(const std::string& name) :
  ConvRHSJacobFluxReconstruction(name),
  socket_alpha("alpha"),
  m_p0Data(),
  m_cellAvgStates(),
  m_cellAvgBackup(),
  m_faceAvgsDone(false),
  m_currFaceAlphaF(0.0),
  m_currCellAlpha(0.0),
  m_flxPntFluxP0(),
  m_flxPntFluxP0Backup(),
  m_meanTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
ConvRHSJacobFluxReconstructionBlending::needsSockets()
{
  std::vector< SafePtr< BaseDataSocketSink > > result = ConvRHSFluxReconstruction::needsSockets();
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "ConvRHSJacobFluxReconstructionBlending::execute()\n");

  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm() || getMethodData().hasArtificialViscosity();

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
  geoDataCell.trs = cells;

  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get number of face orientations
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size()-1;

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder->getDataGE();
  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;

  //// Loop over faces to calculate fluxes and interface fluxes in the flux points

  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
    m_nbrFaceFlxPnts = (*m_faceFlxPntConnPerOrient)[m_orient][0].size();

    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoDataFace.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      // get the neighbouring cells and their states
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();

      const bool leftUpdatable  = (*m_states[LEFT ])[0]->isParUpdatable();
      const bool rightUpdatable = (*m_states[RIGHT])[0]->isParUpdatable();

      // face data and states in the flux points
      if (leftUpdatable || rightUpdatable || hasDiffTerm)
      {
        setFaceData(m_face->getID());
        computeFlxPntStates();
      }

      // interface flux, wave speed and face corrections
      if (leftUpdatable || rightUpdatable)
      {
        computeInterfaceFlxCorrection();
        computeWaveSpeedUpdates(m_waveSpeedUpd);
        updateWaveSpeed();
        computeCorrection(LEFT, m_divContFlxL);
        computeCorrection(RIGHT, m_divContFlxR);
        updateRHSBothSides();
      }

      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
        computeGradientFaceCorrections();
      }

      // face Jacobian (skipped in Jacobian-free mode): dense on blended faces, FR otherwise
      if (getMethodData().doComputeJacobian())
      {
        if (m_currFaceAlphaF > 0.0)
        {
          if (leftUpdatable && rightUpdatable)
          {
            computeFaceJacobBlend(true, LEFT);
          }
          else if (leftUpdatable)
          {
            computeFaceJacobBlend(false, LEFT);
          }
          else if (rightUpdatable)
          {
            computeFaceJacobBlend(false, RIGHT);
          }
        }
        else
        {
          if (leftUpdatable && rightUpdatable)
          {
            computeBothJacobs();
          }
          else if (leftUpdatable)
          {
            computeOneJacob(LEFT);
          }
          else if (rightUpdatable)
          {
            computeOneJacob(RIGHT);
          }
        }
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }

  //// Loop over the elements to calculate the divergence of the continuous flux

  // loop over element types, for the moment there should only be one
  const CFuint nbrElemTypes = elemType->size();
  cf_assert(nbrElemTypes == 1);
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

    // create blockaccumulator
    m_acc.reset(m_lss->createBlockAccumulator(m_nbrSolPnts,m_nbrSolPnts,m_nbrEqs));

    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoDataCell.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
      {
        setCellData();
      }

      // volume term (-divFD+divhFD) scaled by (1 - alpha_K)
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        computeDivDiscontFlx(m_divContFlx);
        updateRHS();
      }

      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
        computeGradients();
      }

      // volume Jacobian: (1 - alpha_K) times the FR one, same sparsity
      if (getMethodData().doComputeJacobian() && (*m_cellStates)[0]->isParUpdatable())
      {
        computeJacobConvCorrection();
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvRHSJacobFluxReconstructionBlending::getCellAlpha(const std::vector< State* >& states)
{
  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  return alpha[states[0]->getLocalID()];
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::computeFaceCellAverages()
{
  if (m_faceAvgsDone) return;

  DataHandle< CFreal > volumes = socket_volumes.getDataHandle();
  m_p0Data.computeCellAverage(*m_states[LEFT ], volumes, *m_cellAvgStates[LEFT ]);
  m_p0Data.computeCellAverage(*m_states[RIGHT], volumes, *m_cellAvgStates[RIGHT]);
  m_faceAvgsDone = true;
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::computeFaceFluxP0()
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // P0 flux of the cell averages with the normal of this flux point
    m_flxPntFluxP0[iFlxPnt] = m_riemannFluxComputer->computeFlux(*m_cellAvgStates[LEFT],
                                                                 *m_cellAvgStates[RIGHT],
                                                                 m_unitNormalFlxPnts[iFlxPnt]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::computeInterfaceFlxCorrection()
{
  // FR interface flux from the extrapolated states, scaled by the face Jacobian in m_cellFlx
  ConvRHSFluxReconstruction::computeInterfaceFlxCorrection();

  m_faceAvgsDone = false;

  // both neighbours evaluate the same expression, so the face flux stays single-valued
  m_currFaceAlphaF = std::max(getCellAlpha(*m_states[LEFT]), getCellAlpha(*m_states[RIGHT]));
  if (m_currFaceAlphaF <= 0.0) return;

  computeFaceCellAverages();
  computeFaceFluxP0();
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::addMeanFaceTerm(const CFuint side, const CFreal alpha,
                                                             vector< RealVector >& corrections)
{
  // -sum_fp W_fp [alpha_f F* + (alpha_K - alpha_f) F_I], where m_cellFlx holds F_I times the
  // signed face Jacobian of this side; with the FR part (1 - alpha_K) of F_I the cell average
  // moves by -sum_fp W_fp Fhat_fp
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

void ConvRHSJacobFluxReconstructionBlending::computeCorrection(CFuint side, vector< RealVector >& corrections)
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

void ConvRHSJacobFluxReconstructionBlending::computeFullFaceCorrection(const CFuint side, RealVector& corrections)
{
  cf_assert(corrections.size() == m_nbrSolPnts*m_nbrEqs);

  corrections = 0.0;

  const CFreal alpha = getCellAlpha(*m_states[side]);
  const CFreal oneMinusAlpha = 1.0 - alpha;

  // FR part -(1 - alpha_K) sum_fp divh_i,fp F_I,fp, from every flux point of the face
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt];
    const RealVector& flx = m_cellFlx[side][iFlxPnt];

    const CFuint nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
    for (CFuint iSolPnt = 0; iSolPnt < nbrSolDep; ++iSolPnt)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];
      const CFreal factor = oneMinusAlpha*m_corrFctDiv[solIdx][flxIdx];

      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        corrections[m_nbrEqs*solIdx+iEq] -= factor*flx[iEq];
      }
    }
  }

  // mean part -sum_fp W_fp [alpha_f F* + (alpha_K - alpha_f) F_I], the same at every solution point
  m_meanTerm = 0.0;
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt];
    const CFreal w = m_p0Data.getW(flxIdx);
    const CFreal factorP0 = w*m_currFaceAlphaF*m_faceJacobVecSizeFlxPnts[iFlxPnt][side];
    const CFreal factorFR = w*(alpha - m_currFaceAlphaF);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_meanTerm[iEq] += factorP0*m_flxPntFluxP0[iFlxPnt][iEq] + factorFR*m_cellFlx[side][iFlxPnt][iEq];
    }
  }

  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      corrections[m_nbrEqs*iSolPnt+iEq] -= m_meanTerm[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::computeFaceJacobBlend(const bool bothSides, const CFuint side)
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_accFace;

  DataHandle< CFreal > volumes = socket_volumes.getDataHandle();

  // block rows and columns: the solution points of both cells
  CFuint solIdx = 0;
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[iSide])[iSol]->getLocalID());
    }
  }

  // sides whose rows are filled
  const CFuint firstSide = bothSides ? LEFT  : side;
  const CFuint lastSide  = bothSides ? RIGHT : side;

  // unperturbed face parts, with the operations of the perturbed ones
  for (CFuint iSide = firstSide; iSide <= lastSide; ++iSide)
  {
    computeFullFaceCorrection(iSide, m_resUpdates[iSide]);
  }

  // store backups of the states in the flux points, F_I and F*
  storeBackups();
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    m_flxPntFluxP0Backup[iFlxPnt] = m_flxPntFluxP0[iFlxPnt];
  }

  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

    // cell average of the perturbed side
    State& pertAvg = *m_cellAvgStates[m_pertSide];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_cellAvgBackup[iEq] = pertAvg[iEq];
    }

    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];

      // flux points of the face whose extrapolated state depends on the perturbed solution point
      m_influencedFlxPnts.resize(0);
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];

        for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
        {
          if (currFlxIdx == (*m_solFlxDep)[m_pertSol][jFlxPnt])
          {
            m_influencedFlxPnts.push_back(iFlxPnt);
          }
        }
      }
      m_NbInfluencedFlxPnts = m_influencedFlxPnts.size();

      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);

        // perturbed F_I at the flux points the solution point reaches
        extrapolatePerturbedState();
        computePertInterfaceFlxCorrection();

        // perturbed cell average and F* at every flux point of the face
        m_p0Data.computeCellAverage(*m_states[m_pertSide], volumes, pertAvg);
        computeFaceFluxP0();

        // every solution point of the filled sides depends on the perturbed one
        for (CFuint iSide2 = firstSide; iSide2 <= lastSide; ++iSide2)
        {
          const CFuint sideTerm = iSide2*m_nbrSolPnts;

          computeFullFaceCorrection(iSide2, m_pertResUpdates[iSide2]);

          // finite difference derivative of the face term, times the residual factor
          m_numJacob->computeDerivative(m_pertResUpdates[iSide2],m_resUpdates[iSide2],m_derivResUpdates);
          m_derivResUpdates *= resFactor;

          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            acc.addValues(iSol+sideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivResUpdates[m_nbrEqs*iSol]);
          }
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);

        // restore the states in the flux points, F_I, the cell average and F*
        restoreFromBackups();
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          pertAvg[iEq] = m_cellAvgBackup[iEq];
        }
        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {
          m_flxPntFluxP0[iFlxPnt] = m_flxPntFluxP0Backup[iFlxPnt];
        }
      }
    }
  }

  // add the values to the jacobian matrix (or direct element blocks)
  getMethodData().assembleJacobBlockFace(acc, m_cells[LEFT]->getID(), m_cells[RIGHT]->getID(), m_nbrSolPnts);

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::computeDivDiscontFlx(vector< RealVector >& residuals)
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

void ConvRHSJacobFluxReconstructionBlending::computePertDivDiscontFlx(RealVector& residuals)
{
  ConvRHSJacobFluxReconstruction::computePertDivDiscontFlx(residuals);

  const CFreal alpha = m_currCellAlpha;
  if (alpha <= 0.0) return;

  residuals *= (1.0 - alpha);
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::setCellData()
{
  ConvRHSFluxReconstruction::setCellData();

  m_currCellAlpha = getCellAlpha(*m_cellStates);
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvRHSJacobFluxReconstructionBlending::computeWaveSpeedP0(const CFuint side)
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

void ConvRHSJacobFluxReconstructionBlending::updateWaveSpeed()
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

void ConvRHSJacobFluxReconstructionBlending::setup()
{
  CFAUTOTRACE;

  ConvRHSJacobFluxReconstruction::setup();

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
  m_cellAvgBackup.resize(m_nbrEqs);
  m_meanTerm.resize(m_nbrEqs);

  // m_nbrFaceFlxPnts is the largest number of flux points of a face here (prisms)
  m_flxPntFluxP0.resize(m_nbrFaceFlxPnts);
  m_flxPntFluxP0Backup.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    m_flxPntFluxP0[iFlxPnt].resize(m_nbrEqs);
    m_flxPntFluxP0Backup[iFlxPnt].resize(m_nbrEqs);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::unsetup()
{
  CFAUTOTRACE;

  for (CFuint side = 0; side < m_cellAvgStates.size(); ++side)
  {
    // the state owns its coordinates (Node built with isOnMesh = false)
    deletePtr(m_cellAvgStates[side]);
  }
  m_cellAvgStates.clear();
  m_p0Data.unsetup();

  ConvRHSJacobFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
