// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/DiffRHSFluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {
    
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffRHSFluxReconstruction,
               FluxReconstructionSolverData,
               FluxReconstructionModule >
diffRHSFluxReconstructionProvider("DiffRHS");

//////////////////////////////////////////////////////////////////////////////
  
DiffRHSFluxReconstruction::DiffRHSFluxReconstruction(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_gradients("gradients"),
  socket_gradientsAV("gradientsAV"),
  socket_posPrev("posPrev"),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  m_diffusiveVarSet(CFNULL),
  m_cellBuilder(CFNULL),
  m_iElemType(),
  m_cell(),
  m_cellStates(),
  m_cellStatesFlxPnt(),
  m_cellFlx(),
  m_solPntsLocalCoords(CFNULL),
  m_flxPntsLocalCoords(CFNULL),
  m_faceFlxPntConnPerOrient(CFNULL),
  m_faceFlxPntConn(CFNULL),
  m_faceConnPerOrient(CFNULL),
  m_faceBuilder(CFNULL),
  m_nbrEqs(),
  m_dim(),
  m_ndimplus(),
  m_orient(),
  m_nbrSolPnts(),
  m_nbrFaceFlxPnts(),
  m_face(),
  m_cells(),
  m_riemannFluxComputer(CFNULL),
  m_corrFctComputer(CFNULL),
  m_corrFct(),
  m_corrFctDiv(),
  m_states(),
  m_flxPntRiemannFlux(),
  m_contFlx(),
  m_divContFlx(),
  m_waveSpeedUpd(),
  m_faceJacobVecAbsSizeFlxPnts(),
  m_faceIntegrationCoefs(CFNULL),
  m_faceMappedCoordDir(CFNULL),
  m_faceLocalDir(CFNULL),
  m_unitNormalFlxPnts(),
  m_faceJacobVecSizeFlxPnts(),
  m_flxPntCoords(),
  m_cellFluxProjVects(),
  m_cellGrads(),
  m_cellGradFlxPnt(),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_solPolyDerivAtSolPnts(CFNULL),
  m_faceInvCharLengths(),
  m_cellVolume(),
  m_cflConvDiffRatio(),
  m_faceFlxPntCellMappedCoords(CFNULL),
  m_extrapolatedFluxes(),
  m_flxPntFlxDim(CFNULL),
  m_avgSol(),
  m_avgGrad(),
  m_flxLocalCoords(CFNULL),
  m_flxSolDep(CFNULL),
  m_solSolDep(CFNULL),
  m_solFlxDep(CFNULL),
  m_nbrSolDep(),
  m_nbrFlxDep(),
  m_nbrSolSolDep(),
  m_dimList(),
  m_tempGrad(),
  m_nbrTotalFlxPnts(),
  m_faceJacobVecs(),
  m_jacobDets(),
  m_order()
  {
    addConfigOptionsTo(this);
    
    m_addRiemannToGradJacob = true;
    //setParameter( "AddRiemannToGradJacob", &m_addRiemannToGradJacob);
    
    m_addRiemannToGradCrossCellJacob = true;
    //setParameter( "AddRiemannToGradCrossCellJacob", &m_addRiemannToGradCrossCellJacob);
    
    m_addFluxToGradCrossCellJacob = true;
    //setParameter( "AddFluxToGradCrossCellJacob", &m_addFluxToGradCrossCellJacob);
  }
  
  
//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
  //options.addConfigOption< bool >("AddRiemannToGradJacob","Boolean telling whether to add Riemann Flux Jacobian to the Gradient.");
  
  //options.addConfigOption< bool >("AddRiemannToGradCrossCellJacob","Boolean telling whether to add the Riemann Flux Jacobian to the cross-cell gradient.");
  
  //options.addConfigOption< bool >("AddFluxToGradCrossCellJacob","Boolean telling whether to add the discontinuous flux Jacobian to the cross-cell gradient.");
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}
  
//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
DiffRHSFluxReconstruction::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_gradients);
  result.push_back(&socket_gradientsAV);
  result.push_back(&socket_posPrev);
  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "DiffRHSFluxReconstruction::execute()\n");
  
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
      
      // compute volume
      m_cellVolume[LEFT] = m_cells[LEFT]->computeVolume();
      m_cellVolume[RIGHT] = m_cells[RIGHT]->computeVolume();
      
      cf_assert(m_cellVolume[LEFT] > 0.0);
      cf_assert(m_cellVolume[RIGHT] > 0.0);
      
      // set the face data
      setFaceData(m_face->getID());//faceID
    
      // compute the left and right states and gradients in the flx pnts
      computeFlxPntStatesAndGrads();
    
      // compute the common interface flux
      computeInterfaceFlxCorrection();

      // compute the wave speed updates
      computeWaveSpeedUpdates(m_waveSpeedUpd);

      // update the wave speed
      updateWaveSpeed();
    
    // if one of the neighbouring cells is parallel updatable, compute the correction flux
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {

    // compute the correction for the left neighbour
    computeCorrection(LEFT, m_divContFlx);

    // update RHS
    updateRHS();
    
    // compute the correction for the right neighbour
    computeCorrection(RIGHT, m_divContFlx);
    
    // update RHS
    updateRHS();
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

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoDataCell.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();
      
      // if the states in the cell are parallel updatable, compute the resUpdates (-divFC)
      if ((*m_cellStates)[0]->isParUpdatable())
      {
    // set the cell data
    setCellData();

    // compute the divergence of the discontinuous flux (-divFD+divhFD)
    computeDivDiscontFlx(m_divContFlx);
      
    // update RHS
        updateRHS();
      }
      
      // divide by the Jacobian to transform the residuals back to the physical domain
      //divideByJacobDet();
      
      // print out the residual updates for debugging
      if(m_cell->getID() == 35)
      {
    CFLog(VERBOSE, "ID  = " << (*m_cellStates)[0]->getLocalID() << "\n");
        CFLog(VERBOSE, "TotalUpdate = \n");
        // get the datahandle of the rhs
        DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            CFLog(VERBOSE, "" << rhs[resID+iVar] << " ");
          }
          //if (abs(rhs[resID+1])>1e-3) CFLog(INFO, "state: " << (*m_cellStates)[iState]->getLocalID() << ", resU: " << rhs[resID+1] << "\n");
          CFLog(VERBOSE,"\n");
          DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
          CFLog(VERBOSE, "UpdateCoeff: " << updateCoeff[(*m_cellStates)[iState]->getLocalID()] << "\n");
        }
      }
      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  
//  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
//
//    for (CFuint i = 0; i < rhs.size()/4; ++i)
//    {
//    CFLog(INFO, "stateID: " << i << ", resV: " << rhs[i*4+2] << "\n");
//    }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::computeInterfaceFlxCorrection()
{
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // compute the average sol and grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
             
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0;
    }
    
    prepareFluxComputation();
     
    // compute the Riemann flux
//     m_flxPntRiemannFlux[iFlxPnt] = m_diffusiveVarSet->getFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0);
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFlux[iFlxPnt]);
     
    // compute FI in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::setFaceData(CFuint faceID)
{
  // compute face Jacobian vectors
  m_faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);

  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector sizes in the flux points
    DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[faceID][iFlxPnt];

    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient][LEFT];
    m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient][RIGHT];

    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
  }
  
  // compute Jacobian determinants
  // Find a way to avoid havin to compute this!
  m_jacobDets[LEFT] = m_cells[LEFT]->computeGeometricShapeFunctionJacobianDeterminant((*m_faceFlxPntCellMappedCoords)[m_orient][LEFT]);
  m_jacobDets[RIGHT] = m_cells[RIGHT]->computeGeometricShapeFunctionJacobianDeterminant((*m_faceFlxPntCellMappedCoords)[m_orient][RIGHT]);

  // compute inverse characteristic lengths
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_faceInvCharLengths[iFlx] = m_faceJacobVecAbsSizeFlxPnts[iFlx]/(m_jacobDets[LEFT][iFlx] + m_jacobDets[RIGHT][iFlx]);
  }
  
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      const CFuint stateID = (*(m_states[iSide]))[iState]->getLocalID();
      m_cellGrads[iSide][iState] = &gradients[stateID];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::computeFlxPntStatesAndGrads()
{
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];
    
    // reset states in flx pnt
    *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;

    // reset the grads in the flx pnts
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) = 0.0;
      *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]) = 0.0;
    }

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdxL = (*m_flxSolDep)[flxPntIdxL][iSol];
      const CFuint solIdxR = (*m_flxSolDep)[flxPntIdxR][iSol];
 
      // add the contributions of the current sol pnt
      *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*(*((*(m_states[LEFT]))[solIdxL]));
      *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*(*((*(m_states[RIGHT]))[solIdxR]));

      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*((*(m_cellGrads[LEFT][solIdxL]))[iVar]);
    *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*((*(m_cellGrads[RIGHT][solIdxR]))[iVar]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::computeDivDiscontFlx(vector< RealVector >& residuals)
{
  // reset the extrapolated fluxes
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrTotalFlxPnts; ++iFlxPnt)
  {
    m_extrapolatedFluxes[iFlxPnt] = 0.0;
  }

  // Loop over solution points to calculate the discontinuous flux.
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_tempGrad[iVar]) = (*(m_cellGrads[0][iSolPnt]))[iVar];
    }

    m_avgSol = *((*m_cellStates)[iSolPnt]->getData());

    prepareFluxComputation();

    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
//       m_contFlx[iSolPnt][iDim] = m_diffusiveVarSet->getFlux(m_avgSol,grad,m_cellFluxProjVects[iDim][iSolPnt],0);
       computeFlux(m_avgSol,m_tempGrad,m_cellFluxProjVects[iDim][iSolPnt],0,m_contFlx[iSolPnt][iDim]);
//       for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//        {
//       if (m_cell->getID() == 11) CFLog(INFO,"iSol: " << iSolPnt << ", iDir: " << iDim << ", iEq: " << iEq << ", state: " << m_avgSol[iEq] <<
//               ", grad: " << (*(m_tempGrad[iEq]))[iDim] << ", n: " << m_cellFluxProjVects[iDim][iSolPnt] << ", flx: " << m_contFlx[iSolPnt][iDim][iEq] << "\n");
//       }

    }

    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];
      const CFuint dim = (*m_flxPntFlxDim)[flxIdx];

      m_extrapolatedFluxes[flxIdx] += (*m_solPolyValsAtFlxPnts)[flxIdx][iSolPnt]*(m_contFlx[iSolPnt][dim]);
    }
  }

  // Loop over solution pnts to calculate the divergence of the discontinuous flux
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // reset the divergence of FC
    residuals[iSolPnt] = 0.0;

    // Loop over solution pnt to count factor of all sol pnt polys
    for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
    {
      const CFuint jSolIdx = (*m_solSolDep)[iSolPnt][jSolPnt];

      // Loop over deriv directions and sum them to compute divergence
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        const CFreal polyCoef = (*m_solPolyDerivAtSolPnts)[iSolPnt][iDir][jSolIdx];

        // Loop over conservative fluxes
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          // Store divFD in the vector that will be divFC
          residuals[iSolPnt][iEq] += polyCoef*(m_contFlx[jSolIdx][iDir+m_ndimplus][iEq]);
          //if (m_cell->getID() == 11) CFLog(INFO,"State: " << iSolPnt << ", jSol: " << jSolIdx << ", iDir: " << iDir << ", var: " << iEq << ", flx: " << m_contFlx[jSolIdx][iDir][iEq] << "\n");
    }
      }
    }

    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];

      // get the divergence of the correction function
      const CFreal divh = m_corrFctDiv[iSolPnt][flxIdx];
  
      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        residuals[iSolPnt][iVar] += -m_extrapolatedFluxes[flxIdx][iVar] * divh;
//if (m_cell->getID() == 11 && iVar == 2) CFLog(INFO,"State: " << iSolPnt << ", iFlx: " << flxIdx << ", var: "
//        << iVar << ", up: " << -m_extrapolatedFluxes[flxIdx][iVar] * divh << ", flux: " << -m_extrapolatedFluxes[flxIdx][iVar] << ", divh: " << divh << "\n");
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::setCellData()
{
  // create a list of the dimensions in which the deriv will be calculated
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_cellFluxProjVects[iDim] = m_cell->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
  }
  
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // get the grads in the current cell
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    const CFuint stateID = (*(m_cellStates))[iState]->getLocalID();
    m_cellGrads[0][iState] = &gradients[stateID];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::updateRHS()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // update rhs
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      rhs[resID+iVar] += resFactor*m_divContFlx[iState][iVar];
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::updateWaveSpeed()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // loop over the sol pnts of both sides to update the wave speeds
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrSolPnts = m_states[iSide]->size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      const CFuint solID = (*m_states[iSide])[iSol]->getLocalID();
      updateCoeff[solID] += m_waveSpeedUpd[iSide]*(2.0*m_order+1);
      CFLog(DEBUG_MIN, "updateCoeff = " << updateCoeff[solID] << "\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::computeCorrection(CFuint side, vector< RealVector >& corrections)
{
  cf_assert(corrections.size() == m_nbrSolPnts);

  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // reset the corrections which will be stored in divContFlx in order to be able to reuse updateRHS()
    corrections[iSolPnt] = 0.0;
  }
  
  // loop over sol pnts to compute the corrections
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt];

    // the current correction factor (stored in cellFlx)
    const RealVector& currentCorrFactor = m_cellFlx[side][iFlxPnt];

    // compute the term due to each flx pnt
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

      // divergence of the correction function
      const CFreal divh = m_corrFctDiv[solIdx][flxIdx];
    
      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        corrections[solIdx][iVar] += currentCorrFactor[iVar] * divh;
//if (m_cells[side]->getID() == 11 && iVar == 2 && flxIdx == 1) CFLog(INFO, "State: " << solIdx << ", flx: " << flxIdx << ", var: " << iVar << ", up: " << currentCorrFactor[iVar] * divh << ", divh: " << divh << "\n");
      }
    }
  }
  
  // in order to use updateRHS, m_cellStates should have the correct states
  m_cellStates = m_cells[side]->getStates();
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::divideByJacobDet()
{
  // this is achieved by multiplying the update coefs with the Jacobian determinant
  // (and dividing by the cell volume)

  // get the updateCoeff
  DataHandle< CFreal > updateCoeff = socket_updateCoeff.getDataHandle();

  // get the cell volume
  const CFreal invCellVolume = 1.0/m_cell->computeVolume();

  // get jacobian determinants at solution points
  const std::valarray<CFreal> jacobDet =
      m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  // get number of solution points
  const CFuint nbrSolPnts = m_cellStates->size();

  // loop over residuals
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // get solution point ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

    // divide update coeff by volume
    updateCoeff[solID] *= jacobDet[iSol]*invCellVolume;
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  CFreal visc = 1.0;
  
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    waveSpeedUpd[iSide] = 0.0;
    for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
    //if (waveSpeedUpd[iSide] > 10.0) CFLog(INFO, "wvspDiff: " << waveSpeedUpd[iSide] << "\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::computeFlux(const RealVector& values, const std::vector< RealVector* >& gradients, const RealVector& normal, const CFreal& radius, RealVector& flux)
{
  flux = m_diffusiveVarSet->getFlux(values,gradients,normal,radius);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::setup()
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "Diff setup\n");
  // setup parent class
  FluxReconstructionSolverCom::setup();
  
  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm();

  if (hasDiffTerm)
  {
    // get the diffusive varset
    m_diffusiveVarSet = getMethodData().getDiffusiveVar();
  }
  
  // get face builder
  m_faceBuilder = getMethodData().getFaceBuilder();
  
  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();
  
  // get the Riemann flux
  m_riemannFluxComputer = getMethodData().getRiemannFlux();
  
  // get the correction function computer
  m_corrFctComputer = getMethodData().getCorrectionFunction();
  
  m_waveSpeedUpd.resize(2);
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_order = static_cast<CFuint>(order);

  const CFGeoShape::Type elemShape = frLocalData[0]->getShape();
  
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
  
  // compute flux point coordinates
  SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  m_nbrFaceFlxPnts = flxLocalCoords->size();
  
  // number of sol points
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();

  cf_assert(m_nbrSolPnts == (frLocalData[0]->getSolPntsLocalCoords())->size());

  // dimensionality and number of equations
  m_dim   = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // get solution point local coordinates
  m_solPntsLocalCoords = frLocalData[0]->getSolPntsLocalCoords();
    
  // get flux point local coordinates
  m_flxPntsLocalCoords = frLocalData[0]->getFlxPntsLocalCoords();
   
  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConnPerOrient = frLocalData[0]->getFaceFlxPntConnPerOrient();
    
  // get the face connectivity per orientation
  m_faceConnPerOrient = frLocalData[0]->getFaceConnPerOrient();
  
  // get the face integration coefficient
  m_faceIntegrationCoefs = frLocalData[0]->getFaceIntegrationCoefs();
  
  // get flux point mapped coordinate directions per orient
  m_faceMappedCoordDir = frLocalData[0]->getFaceMappedCoordDirPerOrient();
  
  // get flux point mapped coordinate directions
  m_faceLocalDir = frLocalData[0]->getFaceMappedCoordDir();
    
  // get the face - flx pnt connectivity
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConn();
  
  // get convective/diffusive CFL ratio
  m_cflConvDiffRatio = frLocalData[0]->getCFLConvDiffRatio();
  
  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();
  
  // get the coefs for derivation of the states in the sol pnts
  m_solPolyDerivAtSolPnts = frLocalData[0]->getCoefSolPolyDerivInSolPnts();
  
  // get face flux point cell mapped coordinates
  m_faceFlxPntCellMappedCoords = frLocalData[0]->getFaceFlxPntCellMappedCoordsPerOrient();
  
  // get the dimension on which to project the flux in a flux point
  m_flxPntFlxDim = frLocalData[0]->getFluxPntFluxDim();

  // get the face local coords of the flux points on one face
  m_flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();

  m_flxSolDep = frLocalData[0]->getFlxPntSolDependency();
  
  m_solSolDep = frLocalData[0]->getSolPntSolDependency();

  m_solFlxDep = frLocalData[0]->getSolPntFlxDependency();

  m_nbrSolDep = ((*m_flxSolDep)[0]).size();
  m_nbrFlxDep = ((*m_solFlxDep)[0]).size();
  m_nbrSolSolDep = ((*m_solSolDep)[0]).size();
  
  // create internal and ghost states
  m_cellStatesFlxPnt.resize(2);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPnt[LEFT].push_back(new State());
    m_cellStatesFlxPnt[RIGHT].push_back(new State());
  }
  
  RealVector dummyCoord;
  dummyCoord.resize(m_dim);
  dummyCoord = 0.0;

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPnt[LEFT][iFlx]->setLocalID(iFlx);
    m_cellStatesFlxPnt[RIGHT][iFlx]->setLocalID(iFlx);
    
    m_cellStatesFlxPnt[LEFT][iFlx]->setSpaceCoordinates(new Node(dummyCoord,false));
    m_cellStatesFlxPnt[RIGHT][iFlx]->setSpaceCoordinates(new Node(dummyCoord,false));
  }
  
  for (CFuint iFlxPnt = 0; iFlxPnt < m_flxPntsLocalCoords->size(); ++iFlxPnt)
  {
    RealVector temp(m_nbrEqs);
    m_extrapolatedFluxes.push_back(temp);
  }
  
  // Resize vectors
  m_cells.resize(2);
  m_states.resize(2);
  m_cellFlx.resize(2);
  m_cellGrads.resize(2);
  m_cellGradFlxPnt.resize(2);
  m_cellVolume.resize(2);
  m_jacobDets.resize(2);
  m_faceJacobVecAbsSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_cellGrads[LEFT].resize(m_nbrSolPnts);
  m_cellGrads[RIGHT].resize(m_nbrSolPnts);
  m_cellGradFlxPnt[LEFT].resize(m_nbrFaceFlxPnts);
  m_cellGradFlxPnt[RIGHT].resize(m_nbrFaceFlxPnts);
  m_cellFlx[LEFT].resize(m_nbrFaceFlxPnts);
  m_cellFlx[RIGHT].resize(m_nbrFaceFlxPnts);
  m_jacobDets[LEFT].resize(m_nbrFaceFlxPnts);
  m_jacobDets[RIGHT].resize(m_nbrFaceFlxPnts);
  m_faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_unitNormalFlxPnts.resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFlux.resize(m_nbrFaceFlxPnts);
  m_contFlx.resize(m_nbrSolPnts);
  m_divContFlx.resize(m_nbrSolPnts);
  m_corrFctDiv.resize(m_nbrSolPnts);
  m_cellFluxProjVects.resize(m_dim+m_ndimplus);
  m_flxPntCoords.resize(m_nbrFaceFlxPnts);
  m_faceInvCharLengths.resize(m_nbrFaceFlxPnts);
  m_avgSol.resize(m_nbrEqs);
  m_avgGrad.resize(m_nbrEqs);
  m_tempGrad.resize(m_nbrEqs);
  m_faceJacobVecs.resize(m_nbrFaceFlxPnts);
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx].resize(m_dim);
    m_faceJacobVecSizeFlxPnts[iFlx].resize(2);
    m_unitNormalFlxPnts[iFlx].resize(m_dim);
    m_cellFlx[LEFT][iFlx].resize(m_nbrEqs);
    m_cellFlx[RIGHT][iFlx].resize(m_nbrEqs);
    m_flxPntRiemannFlux[iFlx].resize(m_nbrEqs);
    m_faceJacobVecs[iFlx].resize(m_dim);
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_cellGradFlxPnt[LEFT][iFlx].push_back(new RealVector(m_dim));
      m_cellGradFlxPnt[RIGHT][iFlx].push_back(new RealVector(m_dim));
    }
  }

  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    m_contFlx[iSolPnt].resize(m_dim+m_ndimplus);
    m_divContFlx[iSolPnt].resize(m_nbrEqs);
    m_corrFctDiv[iSolPnt].resize(m_flxPntsLocalCoords->size());
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_contFlx[iSolPnt][iDim].resize(m_nbrEqs);
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
  
  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
  {
    m_avgGrad[iVar] = new RealVector(m_dim);
    m_tempGrad[iVar] = new RealVector(m_dim);
  }
  
  // compute the divergence of the correction function
  m_corrFctComputer->computeDivCorrectionFunction(frLocalData[0],m_corrFctDiv);

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

  m_nbrTotalFlxPnts = m_flxPntsLocalCoords->size();
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    deletePtr(m_cellStatesFlxPnt[LEFT][iFlx]);
    deletePtr(m_cellStatesFlxPnt[RIGHT][iFlx]);
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      deletePtr(m_cellGradFlxPnt[LEFT][iFlx][iGrad]);
      deletePtr(m_cellGradFlxPnt[RIGHT][iFlx][iGrad]);
    }
    m_cellGradFlxPnt[LEFT][iFlx].clear();
    m_cellGradFlxPnt[RIGHT][iFlx].clear();
  }
  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
  {
    deletePtr(m_avgGrad[iVar]);
  }
  m_avgGrad.clear();
  m_tempGrad.clear();
  m_cellStatesFlxPnt[LEFT].clear();
  m_cellStatesFlxPnt[RIGHT].clear();
  m_cellStatesFlxPnt.clear();
  m_cellGradFlxPnt[LEFT].clear();
  m_cellGradFlxPnt[RIGHT].clear();
  m_cellGradFlxPnt.clear();
  
  // unsetup parent class
  FluxReconstructionSolverCom::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

