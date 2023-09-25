// Copyright (C) 2022 KU Leuven CmPA, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionPoisson/PoissonJacobGradientComputer.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Poisson;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< PoissonJacobGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
PoissonJacobGradientComputerProvider("ConvRHSJacobPoisson");
  
//////////////////////////////////////////////////////////////////////////////
  
PoissonJacobGradientComputer::PoissonJacobGradientComputer(const std::string& name) :
  ConvRHSJacobFluxReconstruction(name),
  m_diffVarSetPoisson(CFNULL)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

void PoissonJacobGradientComputer::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ConvRHSJacobFluxReconstructionPoisson::execute()\n");
  
  // boolean telling whether there is a diffusive term
  //const bool hasDiffTerm = getMethodData().hasDiffTerm() || getMethodData().hasArtificialViscosity();
  
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
  
  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    CFLog(VERBOSE, "Orient = " << m_orient << "\n");
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
    m_nbrFaceFlxPnts = (*m_faceFlxPntConnPerOrient)[m_orient][0].size();

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoDataFace.idx = faceID;
      m_face = m_faceBuilder->buildGE();
      
      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();
      
      // if one of the neighbouring cells is parallel updatable or if the gradients have to be computed, set the bnd face data and compute the discontinuous flx
      //if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable() || hasDiffTerm)
      //{
	// set the bnd face data
        setFaceData(m_face->getID());//faceID

	// compute the left and right states in the flx pnts
        computeFlxPntStates();
      //}

      // if one of the neighbouring cells is parallel updatable, compute the correction flux
      //if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
     // {
	// compute the interface flux
	//computeInterfaceFlxCorrection();
	
	// compute the wave speed updates
        //computeWaveSpeedUpdates(m_waveSpeedUpd);
	
        // update the wave speed
        //updateWaveSpeed();
	
	// compute the correction for the left neighbour
	//computeCorrection(LEFT, m_divContFlxL);

	// compute the correction for the right neighbour
	//computeCorrection(RIGHT, m_divContFlxR);

	// update RHS
	//updateRHSBothSides();
      //}
      
      // if there is a diffusive term, compute the gradients
      //if (hasDiffTerm)
      //{
        computeGradientFaceCorrections();
      //}

//      // compute the contribution to the numerical jacobian
//      if ((*m_states[LEFT])[0]->isParUpdatable() && (*m_states[RIGHT])[0]->isParUpdatable())
//      {
//        computeBothJacobs();
//      }
//      else if ((*m_states[LEFT])[0]->isParUpdatable())
//      {
//        computeOneJacob(LEFT);
//      }
//      else if ((*m_states[RIGHT])[0]->isParUpdatable())
//      {
//        computeOneJacob(RIGHT);
//      }

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

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoDataCell.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // if the states in the cell are parallel updatable or the gradients need to be computed, set the cell data
      //if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
      //{
	// set the cell data
	setCellData();
      //}
      
//      // if the states in the cell are parallel updatable, compute the divergence of the discontinuous flx (-divFD+divhFD)
//      if ((*m_cellStates)[0]->isParUpdatable())
//      {
//	// compute the residual updates (-divFC)
//	computeDivDiscontFlx(m_divContFlx);
//
//	// update RHS
//        updateRHS();
//      }
      
      // if there is a diffusive term, compute the gradients
      //if (hasDiffTerm)
      //{
	computeGradients();
      //}

      // if the states in the cell are parallel updatable, compute the contribution to the numerical jacobian
//      if ((*m_cellStates)[0]->isParUpdatable())
//      {
//	// add the contributions to the Jacobian
//	computeJacobConvCorrection();
//      }

      // divide by the Jacobian to transform the residuals back to the physical domain
      //divideByJacobDet();
      
      // print out the residual updates for debugging
      if(m_cell->getID() == 1944)
      {
	CFLog(VERBOSE, "ID  = " << m_cell->getID() << "\n");
        CFLog(VERBOSE, "ConvUpdate = \n");
        // get the datahandle of the rhs
        DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            CFLog(VERBOSE, "" << rhs[resID+iVar] << " ");
          }
          CFLog(VERBOSE,"\n");
          DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
          CFLog(VERBOSE, "UpdateCoeff: " << updateCoeff[(*m_cellStates)[iState]->getLocalID()] << "\n");
        }
      }
      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
}
  
//////////////////////////////////////////////////////////////////////////////

void PoissonJacobGradientComputer::computeGradients()
{
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_tempStates[iSol] = (*m_cellStates)[iSol]->getData();
    
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[0][iSol][iEq] = 0.0;
    }
  }
  
  m_diffVarSetPoisson->setGradientVars(m_tempStates,m_tempGradTerm,m_nbrSolPnts);
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // Loop over gradient directions
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
	m_projectedCorrL = m_tempGradTerm(iEq,iSolPnt) * m_cellFluxProjVects[iDir][iSolPnt];
	
        // Loop over solution pnts to count factor of all sol pnt polys
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
        { 
          const CFuint jSolIdx = (*m_solSolDep)[iSolPnt][jSolPnt];
          // compute the grad updates
          m_gradUpdates[0][jSolIdx][iEq] += (*m_solPolyDerivAtSolPnts)[jSolIdx][iDir][iSolPnt]*m_projectedCorrL;
	}
      }
    }
  }
  
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  
  // get the volumes
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // get state ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/volumes[solID];

    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradients[solID][iGrad] += m_gradUpdates[0][iSol][iGrad];
      gradients[solID][iGrad] *= invJacobDet;
    }
  }
  
  // if needed, compute the gradients for the artificial viscosity
  if (getMethodData().getUpdateVarStr() == "Cons" && getMethodData().hasArtificialViscosity())
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        //set the grad updates to 0 
        m_gradUpdates[0][iSol][iEq] = 0.0;
      }
    }

    // Loop over solution pnts to calculate the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        // Loop over gradient directions
        for (CFuint iDir = 0; iDir < m_dim; ++iDir)
        {
          m_projectedCorrL = (*((*m_cellStates)[iSolPnt]))[iEq] * m_cellFluxProjVects[iDir][iSolPnt];

          // Loop over solution pnts to count factor of all sol pnt polys
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
          {
            const CFuint jSolIdx = (*m_solSolDep)[iSolPnt][jSolPnt];
	  
            // compute the grad updates
            m_gradUpdates[0][jSolIdx][iEq] += (*m_solPolyDerivAtSolPnts)[jSolIdx][iDir][iSolPnt]*m_projectedCorrL;
	  }
        }
      }
    }
  
    // get the gradientsAV
    DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      // get state ID
      const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

      // inverse Jacobian determinant
      const CFreal invJacobDet = 1.0/volumes[solID];

      // update gradientsAV
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        gradientsAV[solID][iGrad] += m_gradUpdates[0][iSol][iGrad];
        gradientsAV[solID][iGrad] *= invJacobDet;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PoissonJacobGradientComputer::computeGradientFaceCorrections()
{
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
    m_tempStatesR[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
  }
  
  m_diffVarSetPoisson->setGradientVars(m_tempStatesL,m_tempGradTermL,m_nbrFaceFlxPnts);
  m_diffVarSetPoisson->setGradientVars(m_tempStatesR,m_tempGradTermR,m_nbrFaceFlxPnts);
  
  // Loop over solution pnts to reset the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[LEFT][iSolPnt][iEq] = 0.0;
      m_gradUpdates[RIGHT][iSolPnt][iEq] = 0.0;
    }
  }
  
  // compute the face corrections to the gradients
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlx];
    const CFuint flxIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx];
    
    m_nbrSolDep = ((*m_flxSolDep)[flxIdxL]).size();

    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFreal avgSol = (m_tempGradTermL(iEq,iFlx)+m_tempGradTermR(iEq,iFlx))/2.0;
      m_projectedCorrL = (avgSol-m_tempGradTermL(iEq,iFlx))*m_faceJacobVecSizeFlxPnts[iFlx][LEFT]*m_unitNormalFlxPnts[iFlx];
      m_projectedCorrR = (avgSol-m_tempGradTermR(iEq,iFlx))*m_faceJacobVecSizeFlxPnts[iFlx][RIGHT]*m_unitNormalFlxPnts[iFlx];

      // Loop over solution pnts to calculate the grad updates
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
      {
        const CFuint iSolIdxL = (*m_flxSolDep)[flxIdxL][iSolPnt];
        const CFuint iSolIdxR = (*m_flxSolDep)[flxIdxR][iSolPnt];

	/// @todo Check if this is also OK for triangles!!
	m_gradUpdates[LEFT][iSolIdxL][iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdxL][flxIdxL];
	m_gradUpdates[RIGHT][iSolIdxR][iEq] += m_projectedCorrR*m_corrFctDiv[iSolIdxR][flxIdxR];
      }
    }
  }
  
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      // get state ID
      const CFuint solID = (*m_states[iSide])[iSol]->getLocalID();

      // update gradients
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        gradients[solID][iGrad] += m_gradUpdates[iSide][iSol][iGrad];
      }
    }
  }
  
  // if needed, compute the gradients for the artificial viscosity
  if (getMethodData().getUpdateVarStr() == "Cons" && getMethodData().hasArtificialViscosity())
  {
    // Loop over solution pnts to reset the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      // Loop over  variables
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        //set the grad updates to 0 
        m_gradUpdates[LEFT][iSolPnt][iEq] = 0.0;
        m_gradUpdates[RIGHT][iSolPnt][iEq] = 0.0;
      }
    }

    // Loop over  variables
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFuint flxIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlx];
      const CFuint flxIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx];

      m_nbrSolDep = ((*m_flxSolDep)[flxIdxL]).size();

      // compute the face corrections to the gradients
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        const CFreal avgSol = ((*(m_cellStatesFlxPnt[LEFT][iFlx]))[iEq]+(*(m_cellStatesFlxPnt[RIGHT][iFlx]))[iEq])/2.0;
	m_projectedCorrL = (avgSol-(*(m_cellStatesFlxPnt[LEFT][iFlx]))[iEq])*m_faceJacobVecSizeFlxPnts[iFlx][LEFT]*m_unitNormalFlxPnts[iFlx];
	m_projectedCorrR = (avgSol-(*(m_cellStatesFlxPnt[RIGHT][iFlx]))[iEq])*m_faceJacobVecSizeFlxPnts[iFlx][RIGHT]*m_unitNormalFlxPnts[iFlx];

        // Loop over solution pnts to calculate the grad updates
        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint iSolIdxL = (*m_flxSolDep)[flxIdxL][iSolPnt];
          const CFuint iSolIdxR = (*m_flxSolDep)[flxIdxR][iSolPnt];

	  /// @todo Check if this is also OK for triangles!!
	  m_gradUpdates[LEFT][iSolIdxL][iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdxL][flxIdxL];
	  m_gradUpdates[RIGHT][iSolIdxR][iEq] += m_projectedCorrR*m_corrFctDiv[iSolIdxR][flxIdxR];
        }
      }
    }
  
    // get the gradientsAV
    DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        // get state ID
        const CFuint solID = (*m_states[iSide])[iSol]->getLocalID();

        // update gradientsAV
        for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
        {
          gradientsAV[solID][iGrad] += m_gradUpdates[iSide][iSol][iGrad];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PoissonJacobGradientComputer::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvRHSJacobFluxReconstruction::setup();
  
  // get the diffusive varset
  m_diffVarSetPoisson = (getMethodData().getDiffusiveVar()).d_castTo< Physics::Poisson::PoissonDiffVarSet >();
  cf_assert(m_diffVarSetPoisson.isNotNull());
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  
  m_tempGradTerm.resize(m_nbrEqs,m_nbrSolPnts);
  m_tempGradTermL.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermR.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStates.resize(m_nbrSolPnts);
  m_tempStatesL.resize(m_nbrFaceFlxPnts);
  m_tempStatesR.resize(m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

