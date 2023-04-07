#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSJacobFluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {
      
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstruction, 
		       FluxReconstructionSolverData, 
		       FluxReconstructionModule >
DiffBndCorrectionsRHSJacobFluxReconstructionProvider("DiffBndCorrectionsRHSJacob");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstruction::DiffBndCorrectionsRHSJacobFluxReconstruction(const std::string& name) :
  DiffBndCorrectionsRHSFluxReconstruction(name),
  m_cellBuilder(CFNULL),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_faces(),
  m_faceNghbrStates(),
  m_pertResUpdates(),
  m_derivResUpdates(),
  m_gradUpdates(),
  m_pertGrads(),
  m_solJacobDet(),
  m_otherFaceLocalIdxs(),
  m_isFaceOnBoundary(CFNULL),
  m_nghbrCellSide(CFNULL),
  m_currCellSide(CFNULL),
  m_faceOrients(CFNULL),
  m_faceBCIdx(CFNULL),
  m_bcStateComputers(CFNULL),
  m_pertCorrections(),
  m_resUpdates(),
  m_cellGradsBackUp(),
  m_solPolyDerivAtSolPnts(CFNULL),
  m_faceFlxPntConnPerOrient(CFNULL),
  m_pertCellStatesFlxPnt(),
  m_faceMappedCoordDirPO(CFNULL),
  m_pertSol(),
  m_pertVar(),
  m_solFlxDep(CFNULL),
  m_nbrFlxDep(),
  m_cellStatesFlxPntBackup(),
  m_influencedFlxPnt(),
  m_flxPntRiemannFluxBackup(),
  m_cellGradFlxPntBackup(),
  m_dimList(),
  m_ndimplus(),
  m_cellFluxProjVects(),
  m_gradTerm(),
  m_gradTermBefore(),
  m_solSolDep(),
  m_nbrSolSolDep(),
  m_projectedCorr(),
  m_flxLocalCoords(),
  m_faceJacobVecs(),
  m_gradTermFace(),
  m_ghostGradTerm(),
  m_gradTermTemp(),
  m_jacobDet(),
  m_cellStatesFlxPnt2(),
  m_eps()
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstruction::~DiffBndCorrectionsRHSJacobFluxReconstruction()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  DiffBndCorrectionsRHSFluxReconstruction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::executeOnTrs()
{
  CFAUTOTRACE;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current QuadFreeBCFluxReconstruction TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();
  
  CFLog(VERBOSE,"DiffBndCorrectionsRHSJacobFluxReconstruction::executeOnTRS: " << faceTrs->getName() << "\n");

  // get bndFacesStartIdxs from FluxReconstructionMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

  // number of face orientations (should be the same for all TRs)
  cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

  // number of TRs
  const CFuint nbTRs = faceTrs->getNbTRs();
  cf_assert(bndFacesStartIdxs.size() == nbTRs);

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;
  
  // get the geodata of the cell builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCB = m_cellBuilder->getDataGE();
  geoDataCB.trs = cellTrs;
  
  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      CFLog(VERBOSE,"m_orient: " << m_orient << "\n");
      
      // select the correct flx pnts on the face out of all cell flx pnts for the current orient
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        m_flxPntsLocalCoords[iFlx] = (*m_allCellFlxPnts)[(*m_faceFlxPntConn)[m_orient][iFlx]];
      }
      
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // loop over faces with this orientation
      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
        // build the face GeometricEntity
        geoData.idx = faceID;
        m_face = m_faceBuilder->buildGE();
	
	// get the states in the neighbouring cell
        m_cellStates = m_face->getNeighborGeo(0)->getStates();
	
	// compute volume
        m_cellVolume = m_face->getNeighborGeo(0)->computeVolume();
	
	cf_assert(m_cellVolume > 0.0);

        // if cell is parallel updatable, compute the correction flux
        if ((*m_cellStates)[0]->isParUpdatable())
        {
	  // build the neighbouring cell
          const CFuint cellID = m_face->getNeighborGeo(0)->getID();
          geoDataCB.idx = cellID;
          m_intCell = m_cellBuilder->buildGE();
          
          // create a list of the dimensions in which the deriv will be calculated
          for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
          {
            m_cellFluxProjVects[iDim] = m_intCell->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
          }
	  
	  // set the bnd face data
	  setBndFaceData(m_face->getID());//faceID

	  // compute the states, gradients and ghost states, gradients in the flx pnts
	  computeFlxPntStates();

	  // compute FI
          computeInterfaceFlxCorrection();

          // compute the wave speed updates
          computeWaveSpeedUpdates(m_waveSpeedUpd);

          // update the wave speeds
          updateWaveSpeed();

	  // compute the correction -(FI)divh of the bnd face for each sol pnt
          computeCorrection(m_corrections);

	  // update the rhs
          updateRHS();

	  // get all the faces neighbouring the cell
          m_faces = m_intCell->getNeighborGeos();

	  // set the local indexes of the other faces than the current boundary faces
          setOtherFacesLocalIdxs();

          // get the neigbouring states of the other faces
          //setFaceNeighbourStates();

	  // compute solution points Jacobian determinants
          m_solJacobDet = m_intCell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
          
          const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
    
          const CFuint iterFreeze = getMethodData().getFreezeJacobIter();
    
          const CFuint interval = iter - iterFreeze;
      
          if (!getMethodData().freezeJacob() || iter < iterFreeze || interval % getMethodData().getFreezeJacobInterval() == 0)
          {

	    // compute the contribution to the jacobian
            computeJacobDiffBndContribution();
          
          }

	  // release the cell
          m_cellBuilder->releaseGE();
        } 
        
        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computeJacobDiffBndContribution()
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    acc.setRowColIndex(iSol,(*m_cellStates)[iSol]->getLocalID());
    
    // put the perturbed and unperturbed corrections in the correct format
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_resUpdates[m_nbrEqs*iSol+iVar] = m_corrections[iSol][iVar];
    }
  }
  
  // get jacobian determinants at solution points
  m_jacobDet = m_intCell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
  
  /// store values that will be overwritten
  storeBackups();

  // loop over the states in the internal cell to perturb the states
  for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[m_pertSol];
    
    computeCellGradTerm(m_gradTermBefore);

    // Loop over flux points to determine which flx pnts are influenced by the pert
    if (elemShape == CFGeoShape::TRIAG || elemShape == CFGeoShape::TETRA)
    {
      m_influencedFlxPnt = 0;
      m_NbInfluencedFlxPnts = m_nbrFaceFlxPnts;
    }
    else
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
      // get current flx pnt idx
      const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    
      for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
      {
        m_influencedFlxPnt = iFlxPnt;
        m_NbInfluencedFlxPnts= iFlxPnt+1;
        break;
      }
      }
    }

    // loop over the variables in the state
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {
      // perturb physical variable in state
      m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
      
      // compute the perturbed corrected gradients
      computePerturbedGradientsAnalytical();   
//      computePerturbedGradients();

      // compute the perturbed states and ghost states in the flx pnts
      //computeFlxPntStates();
      extrapolatePerturbedState();

      // compute the perturbed interface flx correction
      computeInterfaceFlxCorrection();
      computeCorrection(m_pertCorrections);

      // put the perturbed and unperturbed corrections in the correct format
      for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
	{
          m_pertResUpdates[m_nbrEqs*iState+iVar] = m_pertCorrections[iState][iVar];
        }
      }

      // compute the finite difference derivative of the face term
      m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;

      // add the derivative of the residual updates to the accumulator
      CFuint resUpdIdx = 0;
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
      {
        acc.addValues(iSol,m_pertSol,m_pertVar,&m_derivResUpdates[resUpdIdx]);
      }

      // restore physical variable in state
      m_numJacob->restore(pertState[m_pertVar]);
      
      /// restore overwritten values
      restoreFromBackups();
    }
  }
   //acc.printToScreen();

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::extrapolatePerturbedState()
{ 
  // Loop over flux points to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    // get current flx pnt idx
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
      
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // reset the grads in flx pnts
      *(m_cellGradFlxPnt[iFlxPnt][iEq]) = 0.0;
      
      // extrapolate the grads to current flx pnt
      for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
      {
        const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];
       
        *(m_cellGradFlxPnt[iFlxPnt][iEq]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*((*(m_cellGrads[solIdx]))[iEq]);
      }
    }
  }
   
  // compute ghost gradients
  m_bcStateComputer->computeGhostGradients(m_cellGradFlxPnt,m_flxPntGhostGrads,m_unitNormalFlxPnts,m_flxPntCoords);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::storeBackups()
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    m_cellStatesFlxPntBackup[iFlxPnt] = *(m_cellStatesFlxPnt[iFlxPnt]);
    
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_cellGradFlxPntBackup[iFlxPnt][iVar]) = *(m_cellGradFlxPnt[iFlxPnt][iVar]);
    }
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  { 
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_cellGradsBackUp[iSol][iVar] = (*m_cellGrads[iSol])[iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::restoreFromBackups()
{
  for (CFuint iFlxPnt = m_influencedFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  { 
  (*(m_cellStatesFlxPnt[iFlxPnt]))[m_pertVar] = m_cellStatesFlxPntBackup[iFlxPnt][m_pertVar];  
  *(m_cellGradFlxPnt[iFlxPnt][m_pertVar]) = *(m_cellGradFlxPntBackup[iFlxPnt][m_pertVar]);
  }

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  { 
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      (*m_cellGrads[iSol])[iVar] = m_cellGradsBackUp[iSol][iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::setOtherFacesLocalIdxs()
{
  // get face ID of current boundary face
  const CFuint currFaceID = m_face->getID();

  // loop over the faces of the boundary cell
  const CFuint nbrFaces = m_intCell->nbNeighborGeos();
  cf_assert(m_faces->size() == nbrFaces);
  CFuint iFace = 0;
  for (CFuint faceIdx = 0; faceIdx < nbrFaces; ++faceIdx)
  {
    if ((*m_faces)[faceIdx]->getID() != currFaceID)
    {
      m_otherFaceLocalIdxs[iFace] = faceIdx;
      ++iFace;
    }
  }
  cf_assert(iFace == m_otherFaceLocalIdxs.size());
  cf_assert(iFace == nbrFaces-1);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::setFaceNeighbourStates()
{
  // get neighbour states of other faces
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs.size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[iFace];

    // get neigbouring states
    const CFuint nbrFaceNeighbours = (*m_isFaceOnBoundary)[faceIdx] ? 1 : 2;
    for (CFuint iSide = 0; iSide < nbrFaceNeighbours; ++iSide)
    {
      GeometricEntity* cell = (*m_faces)[faceIdx]->getNeighborGeo(iSide);
      m_faceNghbrStates[iFace][iSide] = cell->getStates();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computePerturbedGradients()
{
  // Add the discontinuous gradient
  computeCellGradTerm(m_gradTerm);
  
  // Loop over solution pnts to reset the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolSolDep; ++iSolPnt)
  {
    const CFuint iSolIdx = (*m_solSolDep)[m_pertSol][iSolPnt];
      
    (*m_cellGrads[iSolIdx])[m_pertVar] = 0.0;
  }
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolSolDep; ++iSolPnt)
  {
    const CFuint iSolIdx = (*m_solSolDep)[m_pertSol][iSolPnt];
      
    // Loop over gradient directions
    for (CFuint iDir = 0; iDir < m_dim; ++iDir)
    {
      // Loop over solution pnts to count factor of all sol pnt polys
      for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
      { 
        const CFuint jSolIdx = (*m_solSolDep)[iSolIdx][jSolPnt];
          
        // project the state on a normal and reuse a RealVector variable of the class to store
        m_projectedCorr = m_gradTerm(m_pertVar,jSolIdx) * m_cellFluxProjVects[iDir+m_ndimplus][jSolIdx];
          
        // compute the grad updates
        (*m_cellGrads[iSolIdx])[m_pertVar] += (*m_solPolyDerivAtSolPnts)[iSolIdx][iDir][jSolIdx]*m_projectedCorr;      
      }
    }
  }   
  
  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
  // Add the contribution of the correction to the gradients for each face
  // compute other face contributions to the gradients
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs.size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[iFace];
    const CFuint orient = (*m_faceOrients)[faceIdx];
    
    // compute face Jacobian vectors
    m_faceJacobVecs = (*m_faces)[faceIdx]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);

    // Loop over flux points to set the normal vectors
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // get face Jacobian vector size
      m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[(*m_faces)[faceIdx]->getID()][iFlxPnt];

      // set unit normal vector
      m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
      
      m_flxPntCoords[iFlxPnt] = (*m_faces)[faceIdx]->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlxPnt]);
    }

    if ((*m_isFaceOnBoundary)[faceIdx])
    {
      // Loop over flux points to extrapolate the states to the flux points
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        *(m_cellStatesFlxPnt2[iFlxPnt]) = 0.0;
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
        for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
        {
            const CFuint iSolIdx = (*m_flxSolDep)[currFlxIdx][iSol];
          *(m_cellStatesFlxPnt2[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSolIdx]*(*((*m_cellStates)[iSolIdx]));
        }
      }
  
      // compute ghost states
      (*m_bcStateComputers)[(*m_faceBCIdx)[faceIdx]]->computeGhostStates(m_cellStatesFlxPnt2,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
      
      computeBndGradTerms2(m_gradTermFace,m_ghostGradTerm);
   
      // compute the boundary face contribution to the gradients
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolSolDep; ++iSolPnt)
      {
        const CFuint iSolIdx = (*m_solSolDep)[m_pertSol][iSolPnt];
        
        for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
        {
          const CFuint flxIdx = (*m_faceFlxPntConn)[faceIdx][iFlx];

          const CFreal avgSol = (m_gradTermFace(m_pertVar,iFlx)+m_ghostGradTerm(m_pertVar,iFlx))/2.0;
	  m_projectedCorr = (avgSol-m_gradTermFace(m_pertVar,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[faceIdx])*m_unitNormalFlxPnts[iFlx];
          
          // compute the bnd correction to the gradients
	  /// @todo Check if this is also OK for triangles!!
	  (*m_cellGrads[iSolIdx])[m_pertVar] += m_projectedCorr*m_corrFctDiv[iSolIdx][flxIdx];
        }
      }
    }
    else
    {
      // cell side with respect to this face
      const CFuint cellSide = (*m_currCellSide)[faceIdx];
      
      // loop over flx pnts to extrapolate the states to the flux points and get the 
      // discontinuous normal flux in the flux points and the Riemann flux
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {     
        // local flux point indices in the left and right cell
        const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[orient][LEFT][iFlxPnt];
        const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[orient][RIGHT][iFlxPnt];
    
        (*(m_pertCellStatesFlxPnt[LEFT][iFlxPnt])) = 0.0;
        (*(m_pertCellStatesFlxPnt[RIGHT][iFlxPnt])) = 0.0;

        // extrapolate the left and right states to the flx pnts
        for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
        {
          const CFuint iSolIdxL = (*m_flxSolDep)[flxPntIdxL][iSol];
          const CFuint iSolIdxR = (*m_flxSolDep)[flxPntIdxR][iSol];
          
          (*(m_pertCellStatesFlxPnt[LEFT][iFlxPnt])) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][iSolIdxL]*(*(*m_faceNghbrStates[iFace][LEFT])[iSolIdxL]);
          (*(m_pertCellStatesFlxPnt[RIGHT][iFlxPnt])) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][iSolIdxR]*(*(*m_faceNghbrStates[iFace][RIGHT])[iSolIdxR]);
        }
      }

      computeFaceGradTerms(m_gradTermFace,m_ghostGradTerm);

      // compute the boundary face contribution to the gradients
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolSolDep; ++iSolPnt)
      {
        const CFuint iSolIdx = (*m_solSolDep)[m_pertSol][iSolPnt];
      
        for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
        {
          const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[orient][cellSide][iFlx];
        
          const CFreal avgSol = (m_gradTermFace(m_pertVar,iFlx)+m_ghostGradTerm(m_pertVar,iFlx))/2.0;

	  if (cellSide == LEFT)
	  {
	    m_projectedCorr = (avgSol-m_gradTermFace(m_pertVar,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDirPO)[orient][cellSide])*m_unitNormalFlxPnts[iFlx];
	  }
	  else
	  {
	    m_projectedCorr = (avgSol-m_ghostGradTerm(m_pertVar,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDirPO)[orient][cellSide])*m_unitNormalFlxPnts[iFlx];
	  }
            
          // compute the face corrections to the gradients
	  /// @todo Check if this is also OK for triangles!!
	  (*m_cellGrads[iSolIdx])[m_pertVar] += m_projectedCorr*m_corrFctDiv[iSolIdx][flxIdx];
        }
      }
    }
  }

  // Add the contribution of the correction of the gradients for this bnd face
  
  // compute face Jacobian vectors
  m_faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
    
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[m_face->getID()][iFlxPnt];

    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
    
    m_flxPntCoords[iFlxPnt] = m_face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlxPnt]);
  }

  for (CFuint iFlxPnt = m_influencedFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  { 
    (*(m_cellStatesFlxPnt[iFlxPnt]))[m_pertVar] = 0.0;
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
        
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint iSolIdx = (*m_flxSolDep)[currFlxIdx][iSol];
      (*(m_cellStatesFlxPnt[iFlxPnt]))[m_pertVar] += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSolIdx]*(*((*m_cellStates)[iSolIdx]))[m_pertVar];
    }
  }
  
  // compute ghost states
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
      
  computeBndGradTerms(m_gradTermFace,m_ghostGradTerm);
   
  // compute the boundary face contribution to the gradients
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolSolDep; ++iSolPnt)
  { 
    const CFuint iSolIdx = (*m_solSolDep)[m_pertSol][iSolPnt];
  
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlx];

      const CFreal avgSol = (m_gradTermFace(m_pertVar,iFlx)+m_ghostGradTerm(m_pertVar,iFlx))/2.0;
      m_projectedCorr = (avgSol-m_gradTermFace(m_pertVar,iFlx))*m_faceJacobVecSizeFlxPnts[iFlx]*m_unitNormalFlxPnts[iFlx];
      
      // compute the bnd correction to the gradients
      /// @todo Check if this is also OK for triangles!!
      (*m_cellGrads[iSolIdx])[m_pertVar] += m_projectedCorr*m_corrFctDiv[iSolIdx][flxIdx];
    }
  }

  for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
  { 
    const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];
    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/m_jacobDet[jSolIdx];

    // update gradients
    (*m_cellGrads[jSolIdx])[m_pertVar] *= invJacobDet;
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computePerturbedGradientsAnalytical()
{
  // Add the discontinuous gradient
    
  computeCellGradTerm(m_gradTerm);

  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_eps[iEq] = m_gradTerm(iEq,m_pertSol) - m_gradTermBefore(iEq,m_pertSol);
  }
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolSolDep; ++iSolPnt)
  {
    const CFuint iSolIdx = (*m_solSolDep)[m_pertSol][iSolPnt];
    
    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/m_jacobDet[iSolIdx];
    
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // Loop over gradient directions
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        m_projectedCorr = m_eps[iEq] * m_cellFluxProjVects[iDir+m_ndimplus][m_pertSol];
	  
        // compute the grad updates
        (*m_cellGrads[iSolIdx])[iEq] += (*m_solPolyDerivAtSolPnts)[iSolIdx][iDir][m_pertSol]*m_projectedCorr*invJacobDet;
      }      
    }
  }
  
  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
  // Perturbed flx pnt idx and cell wide idx
  CFuint pertFlxPnt;
  CFuint pertFlxPntIdx;
  
  // Add the contribution of the correction to the gradients for each face
  // compute other face contributions to the gradients
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs.size();

  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[iFace];
    const CFuint orient = (*m_faceOrients)[faceIdx];

    if ((*m_isFaceOnBoundary)[faceIdx])
    {  
      // compute face Jacobian vectors
      m_faceJacobVecs = (*m_faces)[faceIdx]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
    
      // Loop over flux points to set the normal vectors
      if (elemShape == CFGeoShape::TRIAG || elemShape == CFGeoShape::TETRA)
      {
        pertFlxPnt = 0;      
        pertFlxPntIdx = 0;
        m_NbInfluencedFlxPnts = m_nbrFaceFlxPnts;
      }
      else
      {
        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {
          const CFuint currFlxIdx = (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
        
          for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
          {
            if ((*m_solFlxDep)[m_pertSol][jFlxPnt] == currFlxIdx)
            {
              pertFlxPnt = iFlxPnt;
              pertFlxPntIdx = currFlxIdx;
              m_NbInfluencedFlxPnts = pertFlxPnt+1;
              break;
            }
          }
        }
      }


      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
        // get face Jacobian vector size
        m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[(*m_faces)[faceIdx]->getID()][iFlxPnt];

        // set unit normal vector
        m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
        
        m_flxPntCoords[iFlxPnt] = (*m_faces)[faceIdx]->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlxPnt]);
        
        *(m_cellStatesFlxPnt2[iFlxPnt]) = 0.0;
        
        for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
        {
          const CFuint iSolIdx = (*m_flxSolDep)[currFlxIdx][iSol];
            
          *(m_cellStatesFlxPnt2[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSolIdx]*(*((*m_cellStates)[iSolIdx]));
        }
      }
      
      for (CFuint iFlxPnt = pertFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
      {
        pertFlxPntIdx= (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
        // compute ghost states with pert
        (*m_bcStateComputers)[(*m_faceBCIdx)[faceIdx]]->computeGhostStates(m_cellStatesFlxPnt2,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
        
        computeBndGradTerms2(m_gradTermFace,m_ghostGradTerm);
        
        (*(m_cellStatesFlxPnt2[iFlxPnt]))[m_pertVar] -= m_numJacob->getEps() * (*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol];
        
        // compute ghost states without pert
        (*m_bcStateComputers)[(*m_faceBCIdx)[faceIdx]]->computeGhostStates(m_cellStatesFlxPnt2,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
        
        computeBndGradTerms2(m_gradTermFace,m_gradTermTemp);
        
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          ///@todo check if faceLocalDir is ok & faceFlxPntConn
          m_projectedCorr = 0.5*((m_ghostGradTerm(iEq,iFlxPnt)-m_gradTermTemp(iEq,iFlxPnt))-m_eps[iEq]*(*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol])*(m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[faceIdx])*m_unitNormalFlxPnts[iFlxPnt];
          
          // Loop over solution pnts to calculate the grad updates
          for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
          {
            const CFuint iSolIdx = (*m_flxSolDep)[pertFlxPntIdx][iSolPnt];
          
            // inverse Jacobian determinant
            const CFreal invJacobDet = 1.0/m_jacobDet[iSolIdx];

            /// @todo Check if this is also OK for triangles!!
            (*m_cellGrads[iSolIdx])[iEq] += m_projectedCorr*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
  //          if (m_cells[m_pertSide]->getID() == 1) 
  //	  {
  //          RealVector temp = m_projectedCorrL*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
  //            CFLog(INFO,"Ana Bnd: " << iSolIdx << ", "  << temp << "\n");
  //          }
          }
        }
      }
    }
    else
    {   
      // cell side with respect to this face
      const CFuint cellSide = (*m_currCellSide)[faceIdx];
      
      // compute face Jacobian vectors
      m_faceJacobVecs = (*m_faces)[faceIdx]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
      
      // Loop over flux points to set the normal vectors
      if (elemShape == CFGeoShape::TRIAG || elemShape == CFGeoShape::TETRA)
      {
        pertFlxPnt = 0;      
        pertFlxPntIdx = 0;
        m_NbInfluencedFlxPnts = m_nbrFaceFlxPnts;
      }
      else
      {
        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {
          const CFuint currFlxIdx = (*m_faceFlxPntConnPerOrient)[orient][cellSide][iFlxPnt];
        
          for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
          {
            if ((*m_solFlxDep)[m_pertSol][jFlxPnt] == currFlxIdx)
            {
              pertFlxPnt = iFlxPnt;
              pertFlxPntIdx = currFlxIdx;
              m_NbInfluencedFlxPnts = iFlxPnt+1;
              // get face Jacobian vector size
              m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[(*m_faces)[faceIdx]->getID()][iFlxPnt];

              // set unit normal vector
              m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
          
              break;
            }
          }
        }
      }

      for (CFuint iFlxPnt = pertFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
      {
        pertFlxPntIdx= (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          ///@todo check if faceLocalDir is ok & faceFlxPntConn
          m_projectedCorr = -0.5*m_eps[iEq]*(*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol]*(m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDirPO)[orient][cellSide])*m_unitNormalFlxPnts[iFlxPnt];

          // Loop over solution pnts to calculate the grad updates
          for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
          {
            const CFuint iSolIdx = (*m_flxSolDep)[pertFlxPntIdx][iSolPnt];
          
            // inverse Jacobian determinant
            const CFreal invJacobDet = 1.0/m_jacobDet[iSolIdx];

            /// @todo Check if this is also OK for triangles!!
            (*m_cellGrads[iSolIdx])[iEq] += m_projectedCorr*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
  //          if (m_cells[m_pertSide]->getID() == 5) 
  //	  {
  //            RealVector temp = m_projectedCorrL*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
  //              CFLog(INFO,"Ana otherFace: " << iSolIdx << ", "  << temp << "\n");
  //          }
          }
        }
      }
    }
  }
  
  // Add the contribution of the correction of the gradients for this face
  
  // compute face Jacobian vectors
  m_faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
            
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[m_face->getID()][iFlxPnt];

    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
  
    m_flxPntCoords[iFlxPnt] = m_face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlxPnt]);
  }
  
  for (CFuint iFlxPnt = m_influencedFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  { 
  const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];

  *(m_cellStatesFlxPnt[iFlxPnt]) = 0.0;
        
  for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
  {
    const CFuint iSolIdx = (*m_flxSolDep)[currFlxIdx][iSol];
            
    *(m_cellStatesFlxPnt[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSolIdx]*(*((*m_cellStates)[iSolIdx]));
  }
  
  // compute ghost states with pert
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
      
  computeBndGradTerms(m_gradTermFace,m_ghostGradTerm);
      
  (*(m_cellStatesFlxPnt[iFlxPnt]))[m_pertVar] -= m_numJacob->getEps() * (*m_solPolyValsAtFlxPnts)[currFlxIdx][m_pertSol];
      
  // compute ghost states without pert
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
      
  computeBndGradTerms(m_gradTermFace,m_gradTermTemp);
      
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    ///@todo check if faceLocalDir is ok & faceFlxPntConn
    m_projectedCorr = 0.5*((m_ghostGradTerm(iEq,iFlxPnt)-m_gradTermTemp(iEq,iFlxPnt))-m_eps[iEq])*(*m_solPolyValsAtFlxPnts)[currFlxIdx][m_pertSol]*m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*m_unitNormalFlxPnts[iFlxPnt];

    // Loop over solution pnts to calculate the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
    {
      const CFuint iSolIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
    
      // inverse Jacobian determinant
      const CFreal invJacobDet = 1.0/m_jacobDet[iSolIdx];

      /// @todo Check if this is also OK for triangles!!
      (*m_cellGrads[iSolIdx])[iEq] += m_projectedCorr*m_corrFctDiv[iSolIdx][currFlxIdx]*invJacobDet;
    }
  }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      gradTerm(iEq,iFlx) = (*(m_cellStatesFlxPnt[iFlx]->getData()))[iEq];
      ghostGradTerm(iEq,iFlx) = (*(m_flxPntGhostSol[iFlx]->getData()))[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computeBndGradTerms2(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      gradTerm(iEq,iFlx) = (*(m_cellStatesFlxPnt2[iFlx]->getData()))[iEq];
      ghostGradTerm(iEq,iFlx) = (*(m_flxPntGhostSol[iFlx]->getData()))[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computeCellGradTerm(RealMatrix& gradTerm)
{
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iSol = 0; iSol < m_nbrFaceFlxPnts; ++iSol)
    {
      gradTerm(iEq,iSol) = (*((*m_cellStates)[iSol]->getData()))[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      gradTermL(iEq,iFlx) = (*(m_pertCellStatesFlxPnt[LEFT][iFlx]->getData()))[iEq];
      gradTermR(iEq,iFlx) = (*(m_pertCellStatesFlxPnt[RIGHT][iFlx]->getData()))[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;

  // setup parent class
  DiffBndCorrectionsRHSFluxReconstruction::setup();
  
  // get CellToFaceGeBuilder
  m_cellBuilder      = getMethodData().getCellBuilder();
  m_isFaceOnBoundary = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
  m_nghbrCellSide    = m_cellBuilder->getGeoBuilder()->getNeighbrCellSide ();
  m_currCellSide     = m_cellBuilder->getGeoBuilder()->getCurrentCellSide ();
  m_faceOrients      = m_cellBuilder->getGeoBuilder()->getFaceOrient      ();
  m_faceBCIdx        = m_cellBuilder->getGeoBuilder()->getFaceBCIdx       ();

  // get BC state computers
  m_bcStateComputers = getMethodData().getBCStateComputers();

  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);

  //get element shape
  elemShape = frLocalData[0]->getShape();
  
  //Setting ndimplus, needed for Triag (and tetra, prism)
  if (elemShape == CFGeoShape::TRIAG)
    {
      m_ndimplus=3;
    }
  else if (elemShape == CFGeoShape::TETRA)
    {
      m_ndimplus=4;
    }
  else
    {
      m_ndimplus=0;
    }
  // get the coefs for derivation of the states in the sol pnts
  m_solPolyDerivAtSolPnts = frLocalData[0]->getCoefSolPolyDerivInSolPnts();
  
  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConnPerOrient = frLocalData[0]->getFaceFlxPntConnPerOrient();
  
  // get flux point mapped coordinate directions per orient
  m_faceMappedCoordDirPO = frLocalData[0]->getFaceMappedCoordDirPerOrient();
  
  // get the face local coords of the flux points on one face
  m_flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();

  m_solSolDep = frLocalData[0]->getSolPntSolDependency();

  m_nbrSolSolDep = ((*m_solSolDep)[0]).size();
  
  m_solFlxDep = frLocalData[0]->getSolPntFlxDependency();

  m_nbrFlxDep = ((*m_solFlxDep)[0]).size();

  // get the number of faces in a cell
  const CFuint nbrFaces = frLocalData[0]->getNbrCellFaces();
  const CFuint nbrFacesM1 = nbrFaces - 1;

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(m_nbrSolPnts,m_nbrSolPnts,m_nbrEqs));

  // resize m_otherFaceLocalIdxs
  m_otherFaceLocalIdxs.resize(nbrFacesM1);

  // resize m_faceNghbrStates
  m_faceNghbrStates.resize(nbrFacesM1);
  for (CFuint iFace = 0; iFace < nbrFacesM1; ++iFace)
  {
    m_faceNghbrStates[iFace].resize(2);
  }

  // resize variables
  const CFuint nbrCellResiduals = m_nbrSolPnts*m_nbrEqs;
  m_pertResUpdates .resize(nbrCellResiduals);
  m_derivResUpdates.resize(nbrCellResiduals);
  m_resUpdates.resize(nbrCellResiduals);
  m_pertCorrections.resize(m_nbrSolPnts);
  m_cellGradsBackUp.resize(m_nbrSolPnts);
  m_cellGradFlxPntBackup.resize(m_nbrFaceFlxPnts);
  m_cellFluxProjVects.resize(m_dim+m_ndimplus);
  m_gradTerm.resize(m_nbrEqs,m_nbrSolPnts);
  m_gradTermBefore.resize(m_nbrEqs,m_nbrSolPnts);
  m_gradTermFace.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_ghostGradTerm.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_gradTermTemp.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_projectedCorr.resize(m_dim);
  m_faceJacobVecs.resize(m_nbrFaceFlxPnts);
  m_jacobDet.resize(m_nbrSolPnts);
  m_eps.resize(m_nbrEqs);

  // allocate memory for perturbed gradients
  m_pertGrads.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellGradsBackUp[iSol].resize(m_nbrEqs);
    m_pertCorrections[iSol].resize(m_nbrEqs);
    m_pertGrads[iSol] = new vector<RealVector>(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_cellGradsBackUp[iSol][iEq].resize(m_dim);
      (*m_pertGrads[iSol])[iEq].resize(m_dim);
    }
  }

  // resize gradient updates
  m_gradUpdates.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_gradUpdates[iSide].resize(m_nbrSolPnts);
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_gradUpdates[iSide][iSol].resize(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_gradUpdates[iSide][iSol][iEq].resize(m_dim);
      }
    }
  }

  // resize neighbouring cell solution point Jacobian determinants
  m_solJacobDet.resize(m_nbrSolPnts);
  
  // create internal and ghost states
  m_pertCellStatesFlxPnt.resize(2);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_pertCellStatesFlxPnt[LEFT].push_back(new State());
    m_pertCellStatesFlxPnt[RIGHT].push_back(new State());
    m_faceJacobVecs[iFlx].resize(m_dim);
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_pertCellStatesFlxPnt[LEFT][iFlx]->setLocalID(iFlx);
    m_pertCellStatesFlxPnt[RIGHT][iFlx]->setLocalID(iFlx);
  }
  
  m_cellStatesFlxPntBackup.resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFluxBackup.resize(m_nbrFaceFlxPnts);
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPntBackup[iFlx].resize(m_nbrEqs); 
    m_flxPntRiemannFluxBackup[iFlx].resize(m_nbrEqs);
    
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_cellGradFlxPntBackup[iFlx].push_back(new RealVector(m_dim));
    }
  }
  
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_cellFluxProjVects[iDim].resize(m_nbrSolPnts);
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      m_cellFluxProjVects[iDim][iSolPnt].resize(m_dim);
    }
  }
  
  // create a list of the dimensions in which the deriv will be calculated
  m_dimList.resize(m_dim+m_ndimplus);
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_dimList[iDim].resize(m_nbrSolPnts);
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      m_dimList[iDim][iSolPnt] = iDim;
    }
  }
  
  // create internal states
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPnt2.push_back(new State());
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPnt2[iFlx]->setLocalID(iFlx);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  for (CFuint iSol = 0; iSol < m_pertGrads.size(); ++iSol)
  {
    deletePtr(m_pertGrads[iSol]);
  }
  m_pertGrads.clear();
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    deletePtr(m_cellStatesFlxPnt2[iFlx]);
    deletePtr(m_pertCellStatesFlxPnt[LEFT][iFlx]);
    deletePtr(m_pertCellStatesFlxPnt[RIGHT][iFlx]);
  }
  m_pertCellStatesFlxPnt.clear();
  m_cellStatesFlxPnt2.clear();
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      deletePtr(m_cellGradFlxPntBackup[iFlx][iGrad]); 
    }
    m_cellGradFlxPntBackup[iFlx].clear();
  }

  // unsetup parent class
  DiffBndCorrectionsRHSFluxReconstruction::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
