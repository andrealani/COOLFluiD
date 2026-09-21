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
#include "FluxReconstructionPoisson/FluxReconstructionPoisson.hh"
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
		       FluxReconstructionPoissonModule >
PoissonJacobGradientComputerProvider("ConvRHSJacobPoisson");
  
//////////////////////////////////////////////////////////////////////////////
  
PoissonJacobGradientComputer::PoissonJacobGradientComputer(const std::string& name) :
  ConvRHSJacobFluxReconstruction(name)
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

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

