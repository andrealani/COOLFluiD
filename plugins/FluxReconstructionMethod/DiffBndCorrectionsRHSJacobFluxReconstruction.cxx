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
  m_pertCorrections(),
  m_resUpdates(),
  m_gradUpdates(),
  m_pertGrads(),
  m_solJacobDet(),
  m_cellGradsBackUp(),
  m_solPntsLocalCoords(),
  m_otherFaceLocalIdxs(),
  m_pertCellStatesFlxPnt(),
  m_faceMappedCoordDirPO(CFNULL),
  m_faceFlxPntConnPerOrient(CFNULL),
  m_isFaceOnBoundary(CFNULL),
  m_nghbrCellSide(CFNULL),
  m_currCellSide(CFNULL),
  m_faceOrients(CFNULL),
  m_faceBCIdx(CFNULL),
  m_solPolyDerivAtSolPnts(CFNULL),
  m_bcStateComputers(CFNULL)
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
	
	//CFLog(VERBOSE,"faceID: " << faceID << ", real face ID: " << m_face->getID() << "\n");
	//CFLog(VERBOSE,"cellID: " << m_intCell->getID() << "\n");
	//CFLog(VERBOSE,"coord state 0: " << (((*m_cellStates)[0])->getCoordinates()) << "\n");
	//CFLog(VERBOSE,"state 0: " << *(((*m_cellStates)[0])->getData()) << "\n");

        // if cell is parallel updatable, compute the correction flux
        if ((*m_cellStates)[0]->isParUpdatable())
        {
	  // build the neighbouring cell
          const CFuint cellID = m_face->getNeighborGeo(0)->getID();
          geoDataCB.idx = cellID;
          m_intCell = m_cellBuilder->buildGE();
	  
	  // set the face ID in the BCStateComputer
	  m_bcStateComputer->setFaceID(m_face->getID());

	  // set the bnd face data
	  setBndFaceData(m_face->getID());//faceID

	  // compute the states, gradients and ghost states, gradients in the flx pnts
	  computeFlxPntStates();

	  // compute FI-FD
          computeInterfaceFlxCorrection();

          // compute the wave speed updates
          computeWaveSpeedUpdates(m_waveSpeedUpd);
      
          // update the wave speeds
          updateWaveSpeed();

	  // compute the correction -(FI-FD)divh of the bnd face for each sol pnt
          computeCorrection(m_corrections);

	  // update the rhs
          updateRHS();
	  
	  // get all the faces neighbouring the cell
          m_faces = m_intCell->getNeighborGeos();
	  
	  // set the local indexes of the other faces than the current boundary faces
          setOtherFacesLocalIdxs();

          // get the neigbouring states of the other faces
          setFaceNeighbourStates();
	  
	  // compute solution points Jacobian determinants
          m_solJacobDet = m_intCell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
	  
	  // compute the contribution to the jacobian
          computeJacobDiffBndContribution();
	  
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
  }
  
  // put the perturbed and unperturbed corrections in the correct format and backup the gradients in the sol pnts
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_resUpdates[m_nbrEqs*iState+iVar] = m_corrections[iState][iVar];
      m_cellGradsBackUp[iState][iVar] = (*m_cellGrads[iState])[iVar];
    }
  }

  // loop over the states in the internal cell to perturb the states
  for (CFuint iSolPert = 0; iSolPert < m_nbrSolPnts; ++iSolPert)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[iSolPert];

    // loop over the variables in the state
    for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
    {
      // perturb physical variable in state
      m_numJacob->perturb(iEqPert,pertState[iEqPert]);
      
      // compute the perturbed corrected gradients
      computePerturbedGradients();
      
      // compute the perturbed states and ghost states in the flx pnts
      computeFlxPntStates();

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
        acc.addValues(iSol,iSolPert,iEqPert,&m_derivResUpdates[resUpdIdx]);
      }

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);
      
      // restore the gradients in the sol pnts
      for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          (*m_cellGrads[iState])[iVar] = m_cellGradsBackUp[iState][iVar];
        }
      }
    }
  }
//   acc.printToScreen();

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
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
  vector< vector< RealVector > > cellFluxProjVects;
  cellFluxProjVects.resize(m_dim);
  // create a list of the dimensions in which the deriv will be calculated
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    vector<CFuint> dimList;
    dimList.resize(m_nbrSolPnts);
    cellFluxProjVects[iDim].resize(m_nbrSolPnts);
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      dimList[iSolPnt] = iDim;
      cellFluxProjVects[iDim][iSolPnt].resize(m_dim);
    }
    cellFluxProjVects[iDim] = m_intCell->computeMappedCoordPlaneNormalAtMappedCoords(dimList,*m_solPntsLocalCoords);
  }
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[0][iSolPnt][iEq] = 0.0;

      // Loop over gradient directions
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        // Loop over solution pnts to count factor of all sol pnt polys
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolPnts; ++jSolPnt)
        {
	  const RealVector projectedState = ((*(*m_cellStates)[jSolPnt])[iEq]) * cellFluxProjVects[iDir][jSolPnt];
	  
          // compute the grad updates
          m_gradUpdates[0][iSolPnt][iEq] += (*m_solPolyDerivAtSolPnts)[iSolPnt][iDir][jSolPnt]*projectedState;
       
	  if (fabs(m_gradUpdates[0][iSolPnt][iEq][iDir]) < MathTools::MathConsts::CFrealEps())
          {
            m_gradUpdates[0][iSolPnt][iEq][iDir] = 0.0;
	  }
	}
      }
    }
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      (*m_cellGrads[iSol])[iGrad] = m_gradUpdates[0][iSol][iGrad];
    }
  }
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  // compute flux point coordinates
  SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
    
  // compute face Jacobian vectors
  vector< RealVector > faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*flxLocalCoords);
  
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

    // Loop over flux points to set the normal vectors
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // get face Jacobian vector size
      m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[(*m_faces)[faceIdx]->getID()][iFlxPnt];

      // set unit normal vector
      m_unitNormalFlxPnts[iFlxPnt] = faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
    }

    if ((*m_isFaceOnBoundary)[faceIdx])
    {
      // Loop over flux points to extrapolate the states to the flux points
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        *(m_cellStatesFlxPnt[iFlxPnt]) = 0.0;
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          *(m_cellStatesFlxPnt[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSol]*(*((*m_cellStates)[iSol]));
        }
      }
  
      // compute ghost states
      (*m_bcStateComputers)[(*m_faceBCIdx)[faceIdx]]->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
   
      // compute the boundary face contribution to the gradients
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
      {
        // Loop over  variables
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          // set the grad updates to 0 
          m_gradUpdates[0][iSolPnt][iEq] = 0.0;
      
          // compute the bnd correction to the gradients
          for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
          {
            const CFreal avgSol = ((*m_cellStatesFlxPnt[iFlx])[iEq]+(*(m_flxPntGhostSol[iFlx]))[iEq])/2.0;
	    const RealVector projectedCorr = (avgSol-(*m_cellStatesFlxPnt[iFlx])[iEq])*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[orient])*m_unitNormalFlxPnts[iFlx];
	    /// @todo Check if this is also OK for triangles!!
	    m_gradUpdates[0][iSolPnt][iEq] += projectedCorr*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConn)[faceIdx][iFlx]];
          }
        }
      }

      // add the contribution to the gradients
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_cellGrads[iSol])[iEq] += m_gradUpdates[0][iSol][iEq];
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
    
        *(m_pertCellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
        *(m_pertCellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;

        // extrapolate the left and right states to the flx pnts
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          *(m_pertCellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][iSol]*(*(*m_faceNghbrStates[iFace][LEFT])[iSol]);
          *(m_pertCellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][iSol]*(*(*m_faceNghbrStates[iFace][RIGHT])[iSol]);
        }
      }
      
      // compute the boundary face contribution to the gradients
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
      {
        // Loop over  variables
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          //set the grad updates to 0 
          m_gradUpdates[cellSide][iSolPnt][iEq] = 0.0;
      
          // compute the face corrections to the gradients
          for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
          {
            const CFreal avgSol = ((*m_pertCellStatesFlxPnt[LEFT][iFlx])[iEq]+(*m_pertCellStatesFlxPnt[RIGHT][iFlx])[iEq])/2.0;
	    const RealVector projectedCorr = (avgSol-(*m_pertCellStatesFlxPnt[cellSide][iFlx])[iEq])*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDirPO)[orient][cellSide])*m_unitNormalFlxPnts[iFlx];
	    /// @todo Check if this is also OK for triangles!!
	    m_gradUpdates[cellSide][iSolPnt][iEq] += projectedCorr*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConnPerOrient)[orient][cellSide][iFlx]];
          }
        }
      }

      // add the contribution to the gradients
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_cellGrads[iSol])[iEq] += m_gradUpdates[cellSide][iSol][iEq];
        }
      }
    }
  }
  
  // Add the contribution of the correction of the gradients for this bnd face
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[m_face->getID()][iFlxPnt];

    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
  }
    
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    *(m_cellStatesFlxPnt[iFlxPnt]) = 0.0;
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      *(m_cellStatesFlxPnt[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSol]*(*((*m_cellStates)[iSol]));
    }
  }
  
  // compute ghost states
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
   
  // compute the boundary face contribution to the gradients
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // set the grad updates to 0 
      m_gradUpdates[0][iSolPnt][iEq] = 0.0;
      
      // compute the bnd correction to the gradients
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        const CFreal avgSol = ((*m_cellStatesFlxPnt[iFlx])[iEq]+(*(m_flxPntGhostSol[iFlx]))[iEq])/2.0;
        const RealVector projectedCorr = (avgSol-(*m_cellStatesFlxPnt[iFlx])[iEq])*m_faceJacobVecSizeFlxPnts[iFlx]*m_unitNormalFlxPnts[iFlx];
        /// @todo Check if this is also OK for triangles!!
        m_gradUpdates[0][iSolPnt][iEq] += projectedCorr*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConn)[m_orient][iFlx]];
      }
    }
  }

  // get jacobian determinants at solution points
  const std::valarray<CFreal> jacobDet =
      m_intCell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/jacobDet[iSol];

    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      (*m_cellGrads[iSol])[iGrad] += m_gradUpdates[0][iSol][iGrad];
      (*m_cellGrads[iSol])[iGrad] *= invJacobDet;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;

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
  
  // get the coefs for derivation of the states in the sol pnts
  m_solPolyDerivAtSolPnts = frLocalData[0]->getCoefSolPolyDerivInSolPnts();
  
  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConnPerOrient = frLocalData[0]->getFaceFlxPntConnPerOrient();
  
  // get flux point mapped coordinate directions per orient
  m_faceMappedCoordDirPO = frLocalData[0]->getFaceMappedCoordDirPerOrient();

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
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_pertCellStatesFlxPnt[LEFT][iFlx]->setLocalID(iFlx);
    m_pertCellStatesFlxPnt[RIGHT][iFlx]->setLocalID(iFlx);
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
    deletePtr(m_pertCellStatesFlxPnt[LEFT][iFlx]);
    deletePtr(m_pertCellStatesFlxPnt[RIGHT][iFlx]);
  }
  m_pertCellStatesFlxPnt.clear();

  DiffBndCorrectionsRHSFluxReconstruction::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
