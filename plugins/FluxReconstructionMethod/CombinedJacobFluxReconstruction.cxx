// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/BlockAccumulator.hh"
#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

#include "FluxReconstructionMethod/CombinedJacobFluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/GradientVariables.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CombinedJacobFluxReconstruction::CombinedJacobFluxReconstruction(const std::string& name) :
  DiffRHSJacobFluxReconstruction(name),
  m_updateVarSet(CFNULL),
  m_solEpsilons(),
  m_epsilonLR(),
  m_fluxJacobian(),
  m_gradientFluxJacobian(),
  m_riemannFluxJacobian(),
  m_riemannFluxGradJacobian(),
  m_flxPntRiemannFluxDiff(),
  m_pData(),
  m_unpertContFlx(),
  m_pertContFlx(),
  m_derivContFlx(),
  m_derivContFlxSolPnts(),
  m_derivSolPntRes(),
  m_derivFlxPntFlux(),
  m_derivGradVarsSolPnts(),
  m_derivGradVarsFlxPnt(),
  m_allFaceLocalIdxs(),
  m_cellGradsPtrsBackUp(),
  m_physGradVarsSolPntsBefore(),
  m_avGradDerivs(),
  m_avFaceGradDerivs(),
  m_avFaceGradDerivPtrs(),
  m_avBndResDerivs(),
  m_avBndRes(),
  m_pertAVBndRes()
{
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;

  // setup parent class
  DiffRHSJacobFluxReconstruction::setup();

  // flux at one solution point and its derivatives
  m_unpertContFlx.resize(m_dim+m_ndimplus);
  m_pertContFlx.resize(m_dim+m_ndimplus);
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_unpertContFlx[iDim].resize(m_nbrEqs);
    m_pertContFlx[iDim].resize(m_nbrEqs);
  }
  m_derivContFlx.resize(m_nbrEqs);

  m_derivContFlxSolPnts.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_derivContFlxSolPnts[iSol].resize(m_dim+m_ndimplus);
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_derivContFlxSolPnts[iSol][iDim].resize(m_nbrEqs);
    }
  }
  m_derivSolPntRes.resize(m_nbrEqs);
  m_derivFlxPntFlux.resize(m_nbrEqs);

  // derivatives of the gradient variables
  m_derivGradVarsSolPnts.resize(2);
  m_derivGradVarsFlxPnt.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_derivGradVarsSolPnts[iSide].resize(m_nbrEqs,m_nbrSolPnts);
    m_derivGradVarsFlxPnt[iSide].resize(m_nbrEqs,m_nbFaceFlxPntsMax);
  }
  m_physGradVarsSolPntsBefore.resize(m_nbrEqs,m_nbrSolPnts);

  // local indexes of all the faces of a cell
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  const CFuint nbrFaces = frLocalData[0]->getNbrCellFaces();
  m_allFaceLocalIdxs.resize(nbrFaces);
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    m_allFaceLocalIdxs[iFace] = iFace;
  }

  m_cellGradsPtrsBackUp = m_cellGrads;

  // derivatives of the artificial viscosity gradients and boundary residual
  m_avGradDerivs.resize(2);
  m_avFaceGradDerivs.resize(2);
  m_avFaceGradDerivPtrs.resize(2);
  m_avBndResDerivs.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_avGradDerivs[iSide].resize(m_nbrSolPnts);
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_avGradDerivs[iSide][iSol].resize(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_avGradDerivs[iSide][iSol][iEq].resize(m_dim);
      }
    }

    m_avFaceGradDerivs[iSide].resize(m_nbFaceFlxPntsMax);
    m_avFaceGradDerivPtrs[iSide].resize(m_nbFaceFlxPntsMax);
    for (CFuint iFlx = 0; iFlx < m_nbFaceFlxPntsMax; ++iFlx)
    {
      m_avFaceGradDerivs[iSide][iFlx].resize(m_nbrEqs);
      m_avFaceGradDerivPtrs[iSide][iFlx].resize(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_avFaceGradDerivs[iSide][iFlx][iEq].resize(m_dim);
        m_avFaceGradDerivPtrs[iSide][iFlx][iEq] = &m_avFaceGradDerivs[iSide][iFlx][iEq];
      }
    }

    m_avBndResDerivs[iSide].resize(m_nbrSolPnts*m_nbrEqs);
  }
  m_avBndRes.resize(m_nbrSolPnts*m_nbrEqs);
  m_pertAVBndRes.resize(m_nbrSolPnts*m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;

  // unsetup parent class
  DiffRHSJacobFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::computeBothJacobsDiffFaceTerm()
{
  assembleFaceJacobian(-1);
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::computeOneJacobDiffFaceTerm(const CFuint side)
{
  assembleFaceJacobian(side);
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::computeCellFluxJacobians(const CFuint side)
{
  const CFreal resFactor = getMethodData().getResFactor();

  m_pertSide = side;

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_pertSol = iSol;

    // dereference state
    State& pertState = *(*m_states[side])[iSol];

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      *(m_tempGrad[iEq]) = (*(m_cellGrads[side][iSol]))[iEq];
    }

    // unperturbed diffusive minus convective flux
    m_avgSol = pertState;
    prepareSolPntFluxComputation(pertState.getLocalID());
    m_updateVarSet->computePhysicalData(pertState,m_pData);
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[side][iDim][iSol],0,m_unpertContFlx[iDim]);
      m_unpertContFlx[iDim] -= m_updateVarSet->getFlux()(m_pData,m_neighbCellFluxProjVects[side][iDim][iSol]);
    }

    // derivative with respect to the state
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_pertVar = iVar;

      m_numJacob->perturb(iVar,pertState[iVar]);

      m_avgSol = pertState;
      prepareSolPntFluxComputation(pertState.getLocalID());
      m_updateVarSet->computePhysicalData(pertState,m_pData);
      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      {
        computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[side][iDim][iSol],0,m_pertContFlx[iDim]);
        m_pertContFlx[iDim] -= m_updateVarSet->getFlux()(m_pData,m_neighbCellFluxProjVects[side][iDim][iSol]);
      }

      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      {
        m_numJacob->computeDerivative(m_pertContFlx[iDim],m_unpertContFlx[iDim],m_derivContFlx);
        m_fluxJacobian[side][iSol][iVar][iDim] = resFactor*m_derivContFlx;
      }

      m_numJacob->restore(pertState[iVar]);
    }

    // unperturbed diffusive flux
    m_avgSol = pertState;
    prepareSolPntFluxComputation(pertState.getLocalID());
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[side][iDim][iSol],0,m_unpertContFlx[iDim]);
    }

    // derivative with respect to the gradients
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      for (CFuint iGradDim = 0; iGradDim < m_dim; ++iGradDim)
      {
        m_pertVar = iEq;

        m_numJacob->perturb(iEq,(*(m_tempGrad[iEq]))[iGradDim]);

        m_avgSol = pertState;
        prepareSolPntFluxComputation(pertState.getLocalID());
        for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
        {
          computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[side][iDim][iSol],0,m_pertContFlx[iDim]);
        }

        for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
        {
          m_numJacob->computeDerivative(m_pertContFlx[iDim],m_unpertContFlx[iDim],m_derivContFlx);
          m_gradientFluxJacobian[side][iSol][iEq][iGradDim][iDim] = resFactor*m_derivContFlx;
        }

        m_numJacob->restore((*(m_tempGrad[iEq]))[iGradDim]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::assembleFaceJacobian(const CFint ownedSide)
{
  if (!getMethodData().doComputeJacobian())
  {
    return;
  }

  // number of flux points of the current face
  const CFuint nbrFaceFlxPnts = m_nbrFaceFlxPnts;

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;
  acc.reset();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      acc.setRowColIndex(iSide*m_nbrSolPnts+iSol,(*m_states[iSide])[iSol]->getLocalID());
    }
  }

  // partial derivatives of the flux at the solution points and of the common face flux
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    computeCellFluxJacobians(iSide);
  }
  computeRiemannFluxJacobianNum(resFactor);
  computeRiemannFluxToGradJacobianNum(resFactor);

  // back up the gradients and the states of the cells
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_cellGradsBackUp[iSide][iSol] = *(m_cellGrads[iSide][iSol]);
    }
  }
  std::vector< State* >* cellStatesBackUp = m_cellStates;

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // gradient variables of the perturbed cell before the perturbation
    if (hasPhysicalDiffusionJacobian())
    {
      computeCellGradVars(*(m_states[iSide]),m_gradVarsSolPntsBefore);
    }

    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_pertSide = iSide;
        m_pertSol = iSol;
        m_pertVar = iVar;

        const CFreal invEps = computePertCellGradients();

        computePertCompactFaceGradients(nbrFaceFlxPnts,invEps);

        computeAVGradientDerivatives(iSide,iSol,iVar,false);

        for (CFuint destSide = 0; destSide < 2; ++destSide)
        {
          // only the rows of the owned cell when the other cell belongs to another process
          if (ownedSide >= 0 && destSide != static_cast<CFuint>(ownedSide))
          {
            continue;
          }

          // the volume residual of the perturbed cell is added once, the one of the other cell through this face
          if (destSide != iSide || !m_cellFlags[m_cells[destSide]->getID()])
          {
            addCellVolumeJacobian(acc,destSide,invEps);
          }

          addFaceFluxJacobian(destSide,nbrFaceFlxPnts);
        }
      }
    }
  }

  // restore the gradients, the states and the face data
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      *(m_cellGrads[iSide][iSol]) = m_cellGradsBackUp[iSide][iSol];
    }
  }
  m_cellStates = cellStatesBackUp;
  setFaceData(m_face->getID());
  computeFlxPntStatesAndGrads();
  computeInterfaceFlxCorrection();

  // add the values to the jacobian matrix
  getMethodData().assembleJacobBlockFace(acc,m_cells[LEFT]->getID(),m_cells[RIGHT]->getID(),m_nbrSolPnts);

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

CFreal CombinedJacobFluxReconstruction::computePertCellGradients()
{
  // the gradients start from zero, so the perturbed gradients are the change
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_affectedSolPnts[iSide][iSol] = false;

      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        (*m_cellGrads[iSide][iSol])[iEq] = 0.0;
      }
    }
  }

  // dereference state
  State& pertState = *(*m_states[m_pertSide])[m_pertSol];

  // perturb physical variable in state
  m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);

  const CFreal invEps = 1.0/m_numJacob->getEps();

  // volume term and the liftings of every face of the perturbed cell
  if (hasPhysicalDiffusionJacobian())
  {
    DiffRHSJacobFluxReconstruction::computePerturbedGradientsAnalytical(m_pertSide);
  }

  // restore physical variable in state
  m_numJacob->restore(pertState[m_pertVar]);

  return invEps;
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::computePertCompactFaceGradients(const CFuint nbrFaceFlxPnts, const CFreal invEps)
{
  // derivative of the gradient variables at the solution points
  m_derivGradVarsSolPnts[LEFT] = 0.0;
  m_derivGradVarsSolPnts[RIGHT] = 0.0;

  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_derivGradVarsSolPnts[m_pertSide](iEq,m_pertSol) = hasPhysicalDiffusionJacobian() ? m_pertGradVarsChange[iEq]*invEps : 0.0;
  }

  // extrapolated to the flux points of the current face
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
    {
      const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][iSide][iFlx];

      extrapolateGradVarsToFlxPnt(m_derivGradVarsSolPnts[iSide],(*m_flxSolDep)[flxIdx],(*m_solPolyValsAtFlxPnts)[flxIdx],m_nbrEqs,m_flxPntGradVars);

      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_derivGradVarsFlxPnt[iSide](iEq,iFlx) = m_flxPntGradVars[iEq];
      }
    }
  }

  // the compact face gradient is linear in the gradient variables
  computeCompactBR2FaceGradients(m_derivGradVarsSolPnts[LEFT],m_derivGradVarsSolPnts[RIGHT],m_derivGradVarsFlxPnt[LEFT],m_derivGradVarsFlxPnt[RIGHT]);
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::addCellVolumeJacobian(BlockAccumulator& acc, const CFuint destSide, const CFreal invEps)
{
  // derivative of the flux at the solution points: dF = F_U dU + F_q dq
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      RealVector& derivContFlx = m_derivContFlxSolPnts[iSol][iDim];

      derivContFlx = 0.0;

      if (destSide == m_pertSide && iSol == m_pertSol)
      {
        derivContFlx = m_fluxJacobian[destSide][iSol][m_pertVar][iDim];
      }

      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        for (CFuint iGradDim = 0; iGradDim < m_dim; ++iGradDim)
        {
          derivContFlx += m_gradientFluxJacobian[destSide][iSol][iEq][iGradDim][iDim]*((*m_cellGrads[destSide][iSol])[iEq][iGradDim]*invEps);
        }
      }

      addAVCellFluxDerivative(derivContFlx,destSide,iSol,iDim);
    }
  }

  // FR volume operator: Ddiv dF - sum_f C_f E_f dF
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_derivSolPntRes = 0.0;

    for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
    {
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
        m_derivSolPntRes += (*m_solPolyDerivAtSolPnts)[iSol][iDim][jSol]*m_derivContFlxSolPnts[jSol][iDim];
      }
    }

    for (CFuint iFlx = 0; iFlx < m_nbrTotalFlxPnts; ++iFlx)
    {
      m_derivFlxPntFlux = 0.0;

      const CFuint flxDim = (*m_flxPntFlxDim)[iFlx];

      for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
      {
        m_derivFlxPntFlux += (*m_solPolyValsAtFlxPnts)[iFlx][jSol]*m_derivContFlxSolPnts[jSol][flxDim];
      }

      m_derivSolPntRes -= m_corrFctDiv[iSol][iFlx]*m_derivFlxPntFlux;
    }

    addAVBndResidualDerivative(m_derivSolPntRes,destSide,iSol);

    acc.addValues(destSide*m_nbrSolPnts+iSol,m_pertSide*m_nbrSolPnts+m_pertSol,m_pertVar,&m_derivSolPntRes[0]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::addFaceFluxJacobian(const CFuint destSide, const CFuint nbrFaceFlxPnts)
{
  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint pertFlxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlx];
    const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][destSide][iFlx];

    // derivative of the common face flux: F_U E_f dU + F_q dq_avg
    m_derivFlxPntFlux = m_riemannFluxJacobian[m_pertSide][iFlx][m_pertVar]*(*m_solPolyValsAtFlxPnts)[pertFlxIdx][m_pertSol];

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
        m_derivFlxPntFlux += m_riemannFluxGradJacobian[iFlx][iEq][iDim]*(0.5*((*m_cellGradFlxPnt[LEFT][iFlx][iEq])[iDim]+(*m_cellGradFlxPnt[RIGHT][iFlx][iEq])[iDim]));
      }
    }

    addAVFaceFluxDerivative(m_derivFlxPntFlux,iFlx);

    // face correction of every solution point of the destination cell
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_derivSolPntRes = m_derivFlxPntFlux*(m_faceJacobVecSizeFlxPnts[iFlx][destSide]*m_corrFctDiv[iSol][flxIdx]);

      acc.addValues(destSide*m_nbrSolPnts+iSol,m_pertSide*m_nbrSolPnts+m_pertSol,m_pertVar,&m_derivSolPntRes[0]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::computeCellsWithoutInnerFace()
{
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the cell builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCell = m_cellBuilders[LEFT]->getDataGE();
  geoDataCell.trs = cells;

  const CFuint startIdx = (*elemType)[0].getStartIdx();
  const CFuint endIdx = (*elemType)[0].getEndIdx();

  for (CFuint cellID = startIdx; cellID < endIdx; ++cellID)
  {
    if (m_cellFlags[cellID])
    {
      continue;
    }

    // build the cell
    geoDataCell.idx = cellID;
    m_cells[LEFT] = m_cellBuilders[LEFT]->buildGE();
    m_states[LEFT] = m_cells[LEFT]->getStates();

    prepareIsolatedCellAV(cellID);

    if ((*m_states[LEFT])[0]->isParUpdatable())
    {
      assembleIsolatedCellJacobian(cellID);
    }

    // release the cell
    m_cellBuilders[LEFT]->releaseGE();

    m_cellFlags[cellID] = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::assembleIsolatedCellJacobian(const CFuint cellID)
{
  m_cellStates = m_states[LEFT];

  // metrics of the cell
  m_solJacobDet[LEFT] = m_cells[LEFT]->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_neighbCellFluxProjVects[LEFT][iDim] = m_cells[LEFT]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
  }

  // the gradients corrected with all the faces of the cell
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellGrads[LEFT][iSol] = &gradients[(*m_states[LEFT])[iSol]->getLocalID()];
  }

  // volume residual
  computeUnpertCellDiffResiduals(LEFT);
  updateRHSUnpertCell(LEFT);

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  const CFuint iterFreeze = getMethodData().getFreezeJacobIter();

  if (!getMethodData().doComputeJacobian() ||
      (getMethodData().freezeJacob() && iter >= iterFreeze && (iter - iterFreeze)%getMethodData().getFreezeJacobInterval() != 0))
  {
    return;
  }

  // partial derivatives of the flux at the solution points
  computeCellFluxJacobians(LEFT);

  m_faces[LEFT] = m_cells[LEFT]->getNeighborGeos();

  // dereference the single cell accumulator
  BlockAccumulator& acc = *m_accSC;
  acc.reset();

  // back up the gradients
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellGradsBackUp[LEFT][iSol] = *(m_cellGrads[LEFT][iSol]);
    acc.setRowColIndex(iSol,(*m_states[LEFT])[iSol]->getLocalID());
  }

  // gradient variables of the cell before the perturbation
  if (hasPhysicalDiffusionJacobian())
  {
    computeCellGradVars(*(m_states[LEFT]),m_gradVarsSolPntsBefore);
  }

  cf_assert(m_allFaceLocalIdxs.size() == m_cells[LEFT]->nbNeighborGeos());

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_pertSide = LEFT;
      m_pertSol = iSol;
      m_pertVar = iVar;

      // the gradients start from zero, so the perturbed gradients are the change
      for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_cellGrads[LEFT][jSol])[iEq] = 0.0;
        }
      }

      // dereference state
      State& pertState = *(*m_states[LEFT])[iSol];

      // perturb physical variable in state
      m_numJacob->perturb(iVar,pertState[iVar]);

      const CFreal invEps = 1.0/m_numJacob->getEps();

      // volume term and the liftings of every face of the cell
      if (hasPhysicalDiffusionJacobian())
      {
        addPerturbedVolumeGradient(LEFT);
        addPerturbedFaceLiftings(LEFT,m_allFaceLocalIdxs);
      }

      // restore physical variable in state
      m_numJacob->restore(pertState[iVar]);

      computeAVGradientDerivatives(LEFT,iSol,iVar,true);

      addCellVolumeJacobian(acc,LEFT,invEps);
    }
  }

  // restore the gradients
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    *(m_cellGrads[LEFT][iSol]) = m_cellGradsBackUp[LEFT][iSol];
  }

  // add the values to the jacobian matrix
  getMethodData().assembleJacobBlock(acc,cellID);

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::computeAVGradientDerivatives(const CFuint side, const CFuint iSol, const CFuint iVar, const bool isolated)
{
  if (!hasArtificialViscosityJacobian())
  {
    return;
  }

  // number of flux points of the current face
  const CFuint nbrFaceFlxPnts = m_nbrFaceFlxPnts;

  // the perturbed gradients go to the derivative storage, the physical data is kept
  m_cellGradsPtrsBackUp = m_cellGrads;
  m_physGradVarsSolPntsBefore = m_gradVarsSolPntsBefore;

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
    {
      m_cellGrads[iSide][jSol] = &m_avGradDerivs[iSide][jSol];

      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_avGradDerivs[iSide][jSol][iEq] = 0.0;
      }
    }
  }

  m_pertSide = side;
  m_pertSol = iSol;
  m_pertVar = iVar;

  // artificial viscosity gradient variables before the perturbation
  setAVGradientVars(*m_states[side],m_nbrSolPnts,m_gradVarsSolPntsBefore);

  // dereference state
  State& pertState = *(*m_states[side])[iSol];

  // perturb physical variable in state
  m_numJacob->perturb(iVar,pertState[iVar]);

  const CFreal invEps = 1.0/m_numJacob->getEps();

  // change of the gradient corrected with all the faces
  addPerturbedVolumeGradient(side,true);

  if (isolated)
  {
    addPerturbedFaceLiftings(side,m_allFaceLocalIdxs,true);
  }
  else
  {
    addPerturbedFaceLiftings(side,m_otherFaceLocalIdxs[side],true);
    addPerturbedCurrentFaceGradient(side);
  }

  // restore physical variable in state
  m_numJacob->restore(pertState[iVar]);

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
    {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_avGradDerivs[iSide][jSol][iEq] *= invEps;
      }
    }
  }

  // derivative of the compact gradient of the current face
  if (!isolated)
  {
    m_derivGradVarsSolPnts[LEFT] = 0.0;
    m_derivGradVarsSolPnts[RIGHT] = 0.0;

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_derivGradVarsSolPnts[side](iEq,iSol) = m_pertGradVarsChange[iEq]*invEps;
    }

    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
      {
        const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][iSide][iFlx];

        extrapolateGradVarsToFlxPnt(m_derivGradVarsSolPnts[iSide],(*m_flxSolDep)[flxIdx],(*m_solPolyValsAtFlxPnts)[flxIdx],m_nbrEqs,m_flxPntGradVars);

        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_derivGradVarsFlxPnt[iSide](iEq,iFlx) = m_flxPntGradVars[iEq];
        }
      }
    }

    computeCompactBR2FaceGradients(m_derivGradVarsSolPnts[LEFT],m_derivGradVarsSolPnts[RIGHT],m_derivGradVarsFlxPnt[LEFT],m_derivGradVarsFlxPnt[RIGHT],&m_avFaceGradDerivPtrs);
  }

  // restore the gradient pointers and the physical gradient variables
  m_cellGrads = m_cellGradsPtrsBackUp;
  m_gradVarsSolPntsBefore = m_physGradVarsSolPntsBefore;

  // derivative of the boundary artificial viscosity residual: -resFactor*(R(U+dU) - R(U))/dU
  m_avBndResDerivs[LEFT] = 0.0;
  m_avBndResDerivs[RIGHT] = 0.0;

  if (hasAVBoundaryFlux())
  {
    computeBndFacesAVResidual(*m_cells[side],*m_isFaceOnBoundary[side],*m_faceBCIdx[side],m_solEpsilons[side],m_avBndRes);

    m_numJacob->perturb(iVar,pertState[iVar]);

    const CFreal eps = m_numJacob->getEps();

    computeBndFacesAVResidual(*m_cells[side],*m_isFaceOnBoundary[side],*m_faceBCIdx[side],m_solEpsilons[side],m_pertAVBndRes);

    m_numJacob->restore(pertState[iVar]);

    m_avBndResDerivs[side] = (m_pertAVBndRes-m_avBndRes)*(-getMethodData().getResFactor()/eps);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::addAVCellFluxDerivative(RealVector& derivFlux, const CFuint side, const CFuint iSol, const CFuint iDim)
{
  if (!hasArtificialViscosityJacobian())
  {
    return;
  }

  const CFreal coef = -getMethodData().getResFactor()*m_solEpsilons[side][iSol];

  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iGradDim = 0; iGradDim < m_dim; ++iGradDim)
    {
      derivFlux[iEq] += coef*m_avGradDerivs[side][iSol][iEq][iGradDim]*m_neighbCellFluxProjVects[side][iDim][iSol][iGradDim];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::addAVFaceFluxDerivative(RealVector& derivFlux, const CFuint iFlx)
{
  if (!hasArtificialViscosityJacobian())
  {
    return;
  }

  const CFreal coef = -getMethodData().getResFactor()*0.5*(m_epsilonLR[LEFT][iFlx]+m_epsilonLR[RIGHT][iFlx]);

  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      derivFlux[iEq] += coef*0.5*(m_avFaceGradDerivs[LEFT][iFlx][iEq][iDim]+m_avFaceGradDerivs[RIGHT][iFlx][iEq][iDim])*m_unitNormalFlxPnts[iFlx][iDim];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CombinedJacobFluxReconstruction::addAVBndResidualDerivative(RealVector& derivRes, const CFuint side, const CFuint iSol)
{
  if (!hasArtificialViscosityJacobian())
  {
    return;
  }

  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    derivRes[iEq] += m_avBndResDerivs[side][iSol*m_nbrEqs+iEq];
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
