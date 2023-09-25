// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/ConvDiffCLLAVJacobFluxReconstruction.hh"
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
  
ConvDiffCLLAVJacobFluxReconstruction::ConvDiffCLLAVJacobFluxReconstruction(const std::string& name) :
  ConvDiffLLAVJacobFluxReconstruction(name),
  m_corrFctLLAV(),
  m_alphaValues()
  {
    addConfigOptionsTo(this);
    
    m_LLAVCorrFctOrder = 0;
    setParameter( "LLAVCorrFctOrder", &m_LLAVCorrFctOrder);
    
    m_LLAVCorrFctFactor = 0.0;
    setParameter( "LLAVCorrFctFactor", &m_LLAVCorrFctFactor);
  }

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  ConvDiffLLAVJacobFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("LLAVCorrFctOrder","VCJH corr fct used for LLAV that will be used is tied to an FR scheme of this order (so corr fct will be this order + 1).");
  
  options.addConfigOption< CFreal >("LLAVCorrFctFactor","Factor for the VCJH corr fct used for LLAV.");
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstruction::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ConvDiffCLLAVJacobFluxReconstruction::execute()\n");
  
  ////////////////////INITIALIZATION/////////////////////////
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
  geoDataCell.trs = cells;
  
  // reset cell flags, these tell whether the cell flux part of the cell has been added yet
  for (CFuint iCell = 0; iCell < m_cellFlags.size(); ++iCell)
  {
    m_cellFlags[iCell] = false;
  }
  
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
  
  // get the geodata of the cell builders and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCBL = m_cellBuilders[LEFT]->getDataGE();
  geoDataCBL.trs = cells;
  CellToFaceGEBuilder::GeoData& geoDataCBR = m_cellBuilders[RIGHT]->getDataGE();
  geoDataCBR.trs = cells;
  
  // reset epsilon in vertices
  m_nodeEpsilons = 0.0;
  
  // get current residual
  const CFreal residual = SubSystemStatusStack::getActive()->getResidual();
  
  // get current iteration
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  // check if LLAV should be frozen
  m_useMax = residual < m_freezeLimiterRes || iter > m_freezeLimiterIter;
  
  // initialize Smax and eps_total
  m_Smax = -100.0;
  m_totalEps = 0.0;
  
  ////////////////////COMPUTE EPSILON AND GRADIENTS/////////////////////////
  
  //// Loop over the elements to compute the artificial viscosities and cell
  //// part of gradients
  
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
      m_elemIdx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();
      
//       // if the states in the cell are parallel updatable, compute the resUpdates (-divFC)
//       if ((*m_cellStates)[0]->isParUpdatable())
//       {
	// compute the states projected on order P-1
	computeProjStates(m_statesPMinOne);
	
	// compute the artificial viscosity
	computeEpsilon();
	
	// store epsilon
	storeEpsilon();
//       } 
      
      // add the cell part to the gradients
      computeGradients();
      
      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  
  //// Loop over faces to add the face part to the gradients and AV
  
  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoDataFace.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
      m_nbrFaceFlxPnts = (*m_faceFlxPntConnPerOrient)[m_orient][0].size();

      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();
 
      // set the needed data
      setFaceDataForGradients(m_face->getID());
          
      // compute the face correction term of the corrected gradients
      computeGradientFaceCorrections();
      
      // compute solution points epsilons
      for (CFuint iSide = 0; iSide < 2; ++iSide)
      {   
        const CFuint iOtherSide = iSide == LEFT ? RIGHT : LEFT;
            
        // loop over flx pnts to correct the AV
        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {     
          // local flux point indices
          const CFuint flxPntIdx = (*m_faceFlxPntConnPerOrient)[m_orient][iSide][iFlxPnt];

          // correct sol AV
          for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
          {
            const CFuint solIdx = (*m_flxSolDep)[flxPntIdx][iSol];
            
            const CFuint solID = (*(m_states[iSide]))[solIdx]->getLocalID();
            
//            const CFreal coord = ((*(m_states[iSide]))[solIdx]->getCoordinates())[YY];
//
//            if (coord < 0.32 && coord > 0.0)
//            {
//                CFLog(INFO, "coord: " << (*(m_states[iSide]))[solIdx]->getCoordinates()[XX] << " , " << coord << ", eps: " << m_nodeEpsilons[solID] << ", eps after: " << m_nodeEpsilons[solID]+m_corrFctLLAV[solIdx][flxPntIdx]*0.5*(m_cellEpsilons[m_cells[iOtherSide]->getID()]-m_cellEpsilons[m_cells[iSide]->getID()]) << ", corr: " << m_corrFctLLAV[solIdx][flxPntIdx] << ", otherE: " << m_cellEpsilons[m_cells[iOtherSide]->getID()] << ", thisE: " << m_cellEpsilons[m_cells[iSide]->getID()] << "\n");
//            }
            
            //m_nodeEpsilons[solID] += m_LLAVRelax*m_corrFctLLAV[solIdx][flxPntIdx]*0.5*(m_cellEpsilons[m_cells[iOtherSide]->getID()]+m_cellEpsilons[m_cells[iSide]->getID()]);
            //m_nodeEpsilons[solID] += m_LLAVRelax*m_corrFctLLAV[solIdx][flxPntIdx]*0.5*(m_cellEpsilons[m_cells[iOtherSide]->getID()]-m_cellEpsilons[m_cells[iSide]->getID()]);
            m_nodeEpsilons[solID] += m_LLAVRelax*m_corrFctLLAV[solIdx][flxPntIdx]/3.0*(m_cellEpsilons[m_cells[iOtherSide]->getID()]);
            
            m_nodeEpsilons[solID] = max(0.0,m_nodeEpsilons[solID]);
            
            cf_assert(m_nodeEpsilons[solID] > -0.01);
          }
        }
      }

      // release the cells
      m_cellBuilders[LEFT ]->releaseGE();
      m_cellBuilders[RIGHT]->releaseGE();
      
      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
  
  //// print outputs of LLAV
  if (m_printLLAV && iter%m_showrate == 0)
  {
    const std::string nsp = this->getMethodData().getNamespace();
  
#ifdef CF_HAVE_MPI
    MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
    PE::GetPE().setBarrier(nsp);
    const CFuint count = 1;
    MPI_Allreduce(&m_totalEps, &m_totalEpsGlobal, count, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&m_Smax, &m_SmaxGlobal, count, MPI_DOUBLE, MPI_MAX, comm);
#endif
    
    if (PE::GetPE().GetRank(nsp) == 0) 
    {
      // print total artificial viscosity
      CFLog(INFO, "total eps: " << m_totalEpsGlobal << ", Smax: " << m_SmaxGlobal << "\n");
    }

    PE::GetPE().setBarrier(nsp);
  
    m_Smax = m_SmaxGlobal;
  }
  
  m_flagComputeNbNghb = false;
  
  ////////////////////COMPUTE RHS////////////////////////
  
  // get the cell volumes
  DataHandle< CFreal > cellVolumes = socket_cellVolumes.getDataHandle();
  
  //// Loop over faces to compute the RHS
  
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
      m_faceID = faceID; 
      
      // build the face GeometricEntity
      geoDataFace.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
      m_nbrFaceFlxPnts = (*m_faceFlxPntConnPerOrient)[m_orient][0].size();

      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();
      
      // compute volume
      m_cellVolume[LEFT] = cellVolumes[m_cells[LEFT]->getID()];
      m_cellVolume[RIGHT] = cellVolumes[m_cells[RIGHT]->getID()];
      
      cf_assert(m_cellVolume[LEFT] > 0.0);
      cf_assert(m_cellVolume[RIGHT] > 0.0);
      
//      // if one of the neighbouring cells is parallel updatable, compute the correction flux
//      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
//      {
	// build the neighbouring cells
        const CFuint cellIDL = m_face->getNeighborGeo(LEFT)->getID();
        geoDataCBL.idx = cellIDL;
        m_cells[LEFT] = m_cellBuilders[LEFT ]->buildGE();
        const CFuint cellIDR = m_face->getNeighborGeo(RIGHT)->getID();
        geoDataCBR.idx = cellIDR;
        m_cells[RIGHT] = m_cellBuilders[RIGHT]->buildGE();
        
        // get the sol pnt jacobians
        DataHandle< CFreal > volumes = socket_volumes.getDataHandle();
        
        // get the sol pnt normals
        DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
        
        // compute solution points Jacobian determinants and epsilons
	for (CFuint iSide = 0; iSide < 2; ++iSide)
        {   
          //const CFuint iOtherSide = iSide == LEFT ? RIGHT : LEFT;
            
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            m_solJacobDet[iSide][iSol] = volumes[(*(m_states[iSide]))[iSol]->getLocalID()];
            
            const CFuint solID = (*(m_states[iSide]))[iSol]->getLocalID();
            
            m_solEpsilons[iSide][iSol] = m_nodeEpsilons[solID];
      
            for (CFuint iDim = 0; iDim < m_dim; ++iDim)
            {
              for (CFuint jDim = 0; jDim < m_dim; ++jDim)
              {
                m_neighbCellFluxProjVects[iSide][iDim][iSol][jDim] = solPntNormals[solID*m_dim*m_dim+iDim*m_dim+jDim];
              }
            }
          }
	}

	// set the face data
	setFaceData(m_face->getID());//faceID

	// compute the left and right states and gradients in the flx pnts
	computeFlxPntStatesAndGrads();

	// compute FI
	computeInterfaceFlxCorrection();

	// compute the wave speed updates
        computeWaveSpeedUpdates(m_waveSpeedUpd);

        // update the wave speed
        updateWaveSpeed();

	// compute the correction for the left neighbour
	computeCorrection(LEFT, m_divContFlxL);
	m_divContFlx = m_divContFlxL;

	// update RHS
	updateRHS();
	
	// compute the correction for the right neighbour
	computeCorrection(RIGHT, m_divContFlxR);
	m_divContFlx = m_divContFlxR;
	
	// update RHS
	updateRHS();
        
        CFreal leftAvAlpha = 0.0;
        CFreal rightAvAlpha = 0.0;
        
        if (m_LLAVSubCellRedistribution)
        {       
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            m_alphaValues[LEFT][iSol] = sqrt(((*(m_cellGrads[LEFT][iSol]))[0])[XX]*((*(m_cellGrads[LEFT][iSol]))[0])[XX]+((*(m_cellGrads[LEFT][iSol]))[0])[YY]*((*(m_cellGrads[LEFT][iSol]))[0])[YY]);
            m_alphaValues[RIGHT][iSol] = sqrt(((*(m_cellGrads[RIGHT][iSol]))[0])[XX]*((*(m_cellGrads[RIGHT][iSol]))[0])[XX]+((*(m_cellGrads[RIGHT][iSol]))[0])[YY]*((*(m_cellGrads[RIGHT][iSol]))[0])[YY]);
            leftAvAlpha += m_alphaValues[LEFT][iSol];
            rightAvAlpha += m_alphaValues[RIGHT][iSol];
          }
          
          leftAvAlpha /= m_nbrSolPnts;
          rightAvAlpha /= m_nbrSolPnts;
          
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            m_alphaValues[LEFT][iSol] /= leftAvAlpha;
            m_alphaValues[RIGHT][iSol] /= rightAvAlpha;
          }
       
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            m_solEpsilons[LEFT][iSol] *= m_alphaValues[LEFT][iSol];
            m_solEpsilons[RIGHT][iSol] *= m_alphaValues[RIGHT][iSol];
          }
        }
        
        // compute needed cell contributions: what used to be cell loop is incorporated here!!
        if (!m_cellFlags[cellIDL])// && (*m_states[LEFT ])[0]->isParUpdatable())
        {
          // compute cell contribution 
          computeUnpertCellDiffResiduals(LEFT);
          
	  // update RHS
	  updateRHSUnpertCell(LEFT);
        }
        if (!m_cellFlags[cellIDR])// && (*m_states[RIGHT])[0]->isParUpdatable())
        {
          // compute cell contribution 
          computeUnpertCellDiffResiduals(RIGHT);

	  // update RHS
	  updateRHSUnpertCell(RIGHT);
        }

	// get all the faces neighbouring the cells
        m_faces[LEFT ] = m_cells[LEFT ]->getNeighborGeos();
        m_faces[RIGHT] = m_cells[RIGHT]->getNeighborGeos();

        // set the local indexes of the other faces than the current faces
        setOtherFacesLocalIdxs();
        
        const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
    
        const CFuint iterFreeze = getMethodData().getFreezeJacobIter();
    
        const CFuint interval = iter - iterFreeze;
      
        if (!getMethodData().freezeJacob() || iter < iterFreeze || interval % getMethodData().getFreezeJacobInterval() == 0)
        {
        
          // compute the diffusive face term contribution to the jacobian
          if ((*m_states[LEFT])[0]->isParUpdatable() && (*m_states[RIGHT])[0]->isParUpdatable())
          {
            computeBothJacobsDiffFaceTerm();
          }
          else if ((*m_states[LEFT])[0]->isParUpdatable())
          {
            computeOneJacobDiffFaceTerm(LEFT );
          }
          else if ((*m_states[RIGHT])[0]->isParUpdatable())
          {
            computeOneJacobDiffFaceTerm(RIGHT);
          }
        
        }

        // release the cells
        m_cellBuilders[LEFT ]->releaseGE();
        m_cellBuilders[RIGHT]->releaseGE();
        
        m_cellFlags[cellIDL] = true;
        m_cellFlags[cellIDR] = true;
//      }
      
      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstruction::storeEpsilon()
{
  if (!m_useMax) 
  {
    m_cellEpsilons[m_cell->getID()] = m_epsilon;
    //if(m_cell->getID()==1) CFLog(INFO, "eps: " << newEps << ", eps0: " << m_epsilon0 << ", S: " << m_s << "\n");
    m_totalEps += m_epsilon;
  }
  else
  {
    const CFreal maxEps = max(m_epsilon, m_cellEpsilons[m_cell->getID()]);
    m_cellEpsilons[m_cell->getID()] = maxEps;
    m_totalEps += maxEps;
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
          
    m_nodeEpsilons[stateID] += m_cellEpsilons[m_cell->getID()]/3.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstruction::computeEpsilon()
{
  computeEpsilon0();
  
  computeSmoothness();
  
  if (m_s < m_s0 - m_kappa)
  {
    m_epsilon = 0.0;
  }
  else if (m_s > m_s0 + m_kappa)
  {
    m_epsilon = m_epsilon0;
  }
  else
  {
    m_epsilon = m_epsilon0*0.5*(1.0 + sin(0.5*MathTools::MathConsts::CFrealPi()*(m_s-m_s0)/m_kappa));
  }

  if (m_useWallCutOff)
  {
    // Get the wall distance
    DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
    CFreal centroidDistance = 0.0;
      
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
      centroidDistance += wallDist[stateID];
    }
    
    centroidDistance /= m_nbrSolPnts;
    
    if (centroidDistance < m_wallCutOff) 
    {
      if (centroidDistance < 0.5*m_wallCutOff)
      {
        m_epsilon = 0.0; 
      }
      else
      {
        m_epsilon *= 0.5*(1.0 + sin(0.5*MathTools::MathConsts::CFrealPi()*(centroidDistance-0.75*m_wallCutOff)/(0.25*m_wallCutOff)));
      }
    }
  }
  
  if (m_epsilon < 0.0 || m_epsilon != m_epsilon) 
  {
      CFLog(INFO, "eps: " << m_epsilon << ", eps0: " << m_epsilon0 << ", s: " << m_s << "\n");
      m_epsilon = 0.0;
  }
  
  
  
  
  // add the contribution of the faces
  const CFuint nbrFaces = m_cell->nbNeighborGeos();

  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
    m_nbrFaceFlxPnts = (*m_faceFlxPntConn)[iFace].size();
  
    if (!((*m_isFaceOnBoundaryCell)[iFace]) || m_LLAVBCZero)
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];
        m_nbrSolDep = ((*m_flxSolDep)[currFlxIdx]).size();
        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
          const CFuint dim = (*m_flxPntFlxDim)[currFlxIdx];

           m_extrapolatedFluxes[currFlxIdx] += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(m_contFlx[solIdx][dim]);
        }
      }
    }
    else
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];
        m_nbrSolDep = ((*m_flxSolDep)[currFlxIdx]).size();
        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
          const CFuint dim = (*m_flxPntFlxDim)[currFlxIdx];

           m_extrapolatedFluxes[currFlxIdx] += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(m_contFlxWoLLAV[solIdx][dim]);
           
           const CFuint stateID = (*m_cellStates)[solIdx]->getLocalID();
           CFreal eps;
           
           if (!m_useMax) 
           {
             eps = m_epsilon;
           }
           else
           {
             const CFreal maxEps = max(m_epsilon, m_cellEpsilons[m_cell->getID()]);
             eps = maxEps;
           }

           m_nodeEpsilons[stateID] += m_LLAVRelax*m_corrFctLLAV[solIdx][currFlxIdx]/3.0*eps;
        }
      } 
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstruction::setFaceData(CFuint faceID)
{   
  // get the face flux point normals
  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
  
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_faceJacobVecs[iFlxPnt][iDim] = flxPntNormals[m_face->getID()*m_nbFaceFlxPntsMax*m_dim+iFlxPnt*m_dim+iDim];
    }
  }

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
    m_unitNormalFlxPnts[iFlxPnt] = m_mappedFaceNormalDir*m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
  }

  // compute inverse characteristic lengths
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_faceInvCharLengths[iFlx] = 2.*m_faceJacobVecAbsSizeFlxPnts[iFlx]/(m_cellVolume[LEFT] + m_cellVolume[RIGHT]);
  }
  
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      const CFuint stateID = (*(m_states[iSide]))[iState]->getLocalID();
      m_cellGrads[iSide][iState] = &gradients[stateID];
      m_cellGradsAV[iSide][iState] = &gradientsAV[stateID];
    }
  }
  
//  // loop over flx pnts to extrapolate the states to the flux points
//  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
//  {   
//    // reset the states in the flx pnts
//    m_epsilonLR[LEFT][iFlxPnt] = m_cellEpsilons[m_cells[LEFT]->getID()];
//    m_epsilonLR[RIGHT][iFlxPnt] = m_cellEpsilons[m_cells[RIGHT]->getID()];
//  }
  
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {     
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];
    
    m_epsilonLR[LEFT][iFlxPnt] = 0.0;
    m_epsilonLR[RIGHT][iFlxPnt] = 0.0;

    // extrapolate the left and right states to the flx pnts
    m_nbrSolDep = ((*m_flxSolDep)[flxPntIdxL]).size();
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdxL = (*m_flxSolDep)[flxPntIdxL][iSol];
      const CFuint solIdxR = (*m_flxSolDep)[flxPntIdxR][iSol];
 
      // add the contributions of the current sol pnt
      m_epsilonLR[LEFT][iFlxPnt] += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*m_solEpsilons[LEFT][solIdxL];
      m_epsilonLR[RIGHT][iFlxPnt] += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*m_solEpsilons[RIGHT][solIdxR];
    }
    
    m_epsilonLR[LEFT][iFlxPnt] = max(0.0,m_epsilonLR[LEFT][iFlxPnt]);
    m_epsilonLR[RIGHT][iFlxPnt] = max(0.0,m_epsilonLR[RIGHT][iFlxPnt]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstruction::computeLLAVCellFluxJacobianAna(const CFreal resFactor)
{
  CFLog(VERBOSE, "computeLLAVCellFluxJacobianAna\n");
    
  // reset the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
        {
          for (CFuint iDim = 0; iDim < m_dim; ++iDim)
          { 
            m_llavFluxJacobian[m_pertSide][m_pertSol][m_pertVar][LEFT][jSol][iDim] = 0.0;
            m_llavFluxJacobian[m_pertSide][m_pertSol][m_pertVar][RIGHT][jSol][iDim] = 0.0;
          }
        }
      }
    }
  }

  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;

    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      ///@todo check this!
      CFreal w = 1.0;//0.0;
      CFreal wOther = 0.0;
        
//      for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
//      {
//        // get node local index
//        const CFuint nodeIdx = (*(m_cellNodes[m_pertSide]))[iNode]->getLocalID();
//      
//        w += m_nodePolyValsAtSolPnts[m_pertSol][iNode]/m_nbNodeNeighbors[nodeIdx];
//        
//        for (CFuint jNode = 0; jNode < m_nbrCornerNodes; ++jNode)
//        {
//          const CFuint nodeIdxOther = (*(m_cellNodes[iOtherSide]))[jNode]->getLocalID();
//            
//          if (nodeIdx == nodeIdxOther)
//          {
//            wOther += m_nodePolyValsAtSolPnts[m_pertSol][iNode]/m_nbNodeNeighbors[nodeIdx];
//          }
//        }
//      }

      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {          
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        { 
          CFreal q_n = 0.0;
              
          for (CFuint jDim = 0; jDim < m_dim; ++jDim)
          {
            q_n += (*(m_cellGradsAV[m_pertSide][m_pertSol]))[m_pertVar][jDim] * m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol][jDim];  
          }
          
          const CFreal q_n_w = q_n * w;
          const CFreal q_n_wOther = q_n * wOther;
            
          for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
          {
            for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
            {
              m_llavFluxJacobian[m_pertSide][jSol][iEq][m_pertSide][m_pertSol][iDim][m_pertVar] += q_n_w * m_epsJacobian[m_pertSide][jSol][iEq];
                
              m_llavFluxJacobian[iOtherSide][jSol][iEq][m_pertSide][m_pertSol][iDim][m_pertVar] += q_n_wOther * m_epsJacobian[iOtherSide][jSol][iEq];
            }
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstruction::computeLLAVRiemannFluxJacobianAna(const CFreal resFactor)
{
  CFLog(VERBOSE, "computeLLAVRiemannFluxJacobianAna\n");
  
  // reset the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
        {
          m_llavRiemannFluxJacobian[m_pertSide][m_pertSol][m_pertVar][iFlx] = 0.0;
        }
      }
    }
  }
  
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    const CFint sideTerm = m_pertSide == LEFT ? 1 : -1;
      
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      // local flux point index
      const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlx];
      
      // damping factor
      const CFreal dampFactor = m_dampCoeffDiff*m_faceInvCharLengths[iFlx];
      
      const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlx]+m_epsilonLR[RIGHT][iFlx]);
      
      const CFreal eps_dampF_side = dampFactor * sideTerm * epsilon;
          
      ///@todo check this!
      CFreal w = 1.0;//0.0;
      
//      for (CFuint iNode = 0; iNode < m_faceNodes->size(); ++iNode)
//      {
//        for (CFuint iNodeCell = 0; iNodeCell < m_nbrCornerNodes; ++iNodeCell)
//        {
//	  if ((*m_faceNodes)[iNode]->getLocalID() == (*(m_cellNodes[m_pertSide]))[iNodeCell]->getLocalID())
//	  {
//	    // get node local index
//            const CFuint nodeIdx = (*m_cellNodesConn)(m_cells[m_pertSide]->getID(),iNodeCell);
//	    
//            w += m_nodePolyValsAtFlxPnts[flxIdx][iNodeCell]/m_nbNodeNeighbors[nodeIdx];
//	  }
//	}
//      }
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_avgGradAV[iVar]) = (*(m_cellGradFlxPntAV[LEFT][iFlx][iVar]) + *(m_cellGradFlxPntAV[RIGHT][iFlx][iVar]))/2.0;
        
        // compute damping term for LLAV
        const RealVector dGradVarXNormalAV = ((*(m_cellStatesFlxPnt[LEFT][iFlx]))[iVar] - (*(m_cellStatesFlxPnt[RIGHT][iFlx]))[iVar])*m_unitNormalFlxPnts[iFlx];
        *m_avgGradAV[iVar] -= dampFactor*dGradVarXNormalAV;
            
        CFreal q_n = 0.0;
            
        for (CFuint jDim = 0; jDim < m_dim; ++jDim)
        {
          q_n += (*(m_avgGradAV[iVar]))[jDim] * m_unitNormalFlxPnts[iFlx][jDim];
        }
        
        const CFreal q_n_w = q_n * w;
          
        for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
        {
          for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
          {          
            m_llavRiemannFluxJacobian[m_pertSide][m_pertSol][m_pertVar][iFlx][iVar] += q_n_w * m_epsJacobian[m_pertSide][m_pertSol][m_pertVar];
          }
        }
      }
      
      for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
      {
        const CFreal l_eps_dampF_side = (*m_solPolyValsAtFlxPnts)[flxIdx][m_pertSol] * eps_dampF_side;
          
        for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
        { 
          m_llavRiemannFluxJacobian[m_pertSide][m_pertSol][m_pertVar][iFlx][m_pertVar] -= l_eps_dampF_side;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstruction::computeBothJacobsDiffFaceTerm()
{
  CFLog(VERBOSE, "computeBothJacobsDiffFaceTerm\n");
    
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;
  
  CFuint solIdx = 0;
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[m_pertSide])[iSol]->getLocalID());
    }
  }

  //// compute the needed flux jacobians
  
  initJacobianComputation();
  
  computeCellFluxJacobianNum(resFactor);
  
  computeRiemannFluxJacobianNum(resFactor);
  
  computeFluxToGradJacobianNum(resFactor);
  
  computeRiemannFluxToGradJacobianNum(resFactor);
  
  computeGradToStateJacobianAna();
  
  computeGradVarsToStateJacobianNum();
  
  computeEpsToStateJacobianAna();
  
  computeLLAVCellFluxJacobianAna(resFactor);
  
  computeLLAVRiemannFluxJacobianAna(resFactor);
  
  //// add the total jacobians to the system jacobian
  
  // loop over left and right cell to add the discontinuous (cell) part to the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // make sure this is only done once per cell
    if (!m_cellFlags[m_cells[m_pertSide]->getID()]) 
    {
      // term depending on iSide
      const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

      // loop over the states to which to derive (l)
      for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
      {
        // loop over the variables in the state (k)
        for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
        {
          const CFuint nbDepGradVar = m_nbrVarToGradVarDep[m_pertVar];
            
          // add the discontinuous part of the jacobian related to the sol pnt (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
          {
            const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];
            
            m_tempFlux = 0.;

            // (d)
            for (CFuint iDim = 0; iDim < m_dim; ++iDim)
            {
              const CFreal polyCoef = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][m_pertSol]; 
          
              m_tempFlux += m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][iDim] * polyCoef;
            }
            
            acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
          
          // add the discontinuous gradient part of the jacobian (m)
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
          {
            const CFuint kSolIdx = (*m_solSolDep)[m_pertSol][kSolPnt];
          
            // (i)
            for (CFuint jSol = 0; jSol < m_nbrSolSolDep; ++jSol)
            {
              const CFuint jSolIdx = (*m_solSolDep)[kSolIdx][jSol];
              
              m_tempFlux = 0.0;
                
              // (d)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal dl = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdx];
                
                /// llav jacob to state part //// should actually be added for all jSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][m_pertSol][m_pertVar][m_pertSide][jSolIdx][iDim] * dl;
                
                // (b)
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal dl_dqdu = dl * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                    
                  // (p)
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal dl_dqdu_dudu = m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][var] * dl_dqdu;
                      
                    m_tempFlux += m_gradientFluxJacobian[m_pertSide][kSolIdx][var][jDim][iDim] * dl_dqdu_dudu;
                  }
                  
                  CFreal llavPart = m_solEpsilons[m_pertSide][kSolIdx] * m_neighbCellFluxProjVects[m_pertSide][iDim][kSolIdx][jDim];
                  llavPart *= dl_dqdu;
                  //if(m_cells[0]->getID()==1) CFLog(INFO, "pertSol: " << m_pertSol << ", pertVar: " << m_pertVar << ", llavJC: " << llavPart << "\n");
                  // add part of analytical LLAV jacobian
                  m_tempFlux[m_pertVar] += llavPart;
                }
              }
            
              acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
          
          // add the discontinuous part of the jacobian related to the flx pnt (f)
          for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
          {
            const CFuint flxIdx = (*m_solFlxDep)[m_pertSol][iFlxPnt];
            
            // (df)
            const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
            
            m_temp = m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][dim] * (*m_solPolyValsAtFlxPnts)[flxIdx][m_pertSol];
            
            // add the second part of the discontinuous part of the jacobian (i)
            m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[flxIdx][jSolPnt];

              // get the divergence of the correction function
              const CFreal divh = m_corrFctDiv[jSolIdx][flxIdx];
                           
              m_tempFlux = -m_temp * divh;
            
              acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
          
          // add the second part of the discontinuous gradient part of the jacobian (m)
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
          {
            const CFuint kSolIdx = (*m_solSolDep)[m_pertSol][kSolPnt];
              
            for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
            {   
              const CFuint flxIdx = (*m_solFlxDep)[kSolIdx][iFlxPnt];
              
              // (df)
              const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
              
              const CFreal l = (*m_solPolyValsAtFlxPnts)[flxIdx][kSolIdx];            
              
              // (i)
              m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
              for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
              {
                const CFuint jSolIdx = (*m_flxSolDep)[flxIdx][jSolPnt];
              
                m_tempFlux = 0.0;

                const CFreal divh_l = -m_corrFctDiv[jSolIdx][flxIdx] * l;
                
                /// llav jacob to state part //// actually should loop over all kSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][m_pertSol][m_pertVar][m_pertSide][jSolIdx][dim] * divh_l;
              
                // (b)
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqdu = divh_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                  
                  // (p)
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal divh_l_dqdu_dudu = divh_l_dqdu * m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][var];             
                      
                    m_tempFlux += divh_l_dqdu_dudu * m_gradientFluxJacobian[m_pertSide][kSolIdx][var][jDim][dim];
                  }
                  
                  CFreal llavPart = divh_l_dqdu * m_solEpsilons[m_pertSide][kSolIdx];
                  
                  llavPart *= m_neighbCellFluxProjVects[m_pertSide][dim][kSolIdx][jDim];
                  //if(m_cells[0]->getID()==1) CFLog(INFO, "pertSol: " << m_pertSol << ", pertVar: " << m_pertVar << ", llavJF: " << llavPart << "\n");
                  // add part of analytical LLAV Jacobian
                  m_tempFlux[m_pertVar] += llavPart;
                }
                
                acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
              }
            }
          }
        }
      }
    }
  }
  
  // loop over left and right cell to add the riemann flux (face) part to the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // term depending on iSide
    const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

    // term depending on iOtherSide
    const CFuint otherSideTerm = iOtherSide*m_nbrSolPnts;
    
    // loop over the variables in the state (k)
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    { 
      const CFuint nbDepGradVar = m_nbrVarToGradVarDep[m_pertVar];
        
      // loop over face flx pnts (f)
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // local flux point indices in the left and right cell
        const CFuint flxPntIdxThis = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
        const CFuint flxPntIdxOther = (*m_faceFlxPntConnPerOrient)[m_orient][iOtherSide][iFlxPnt];
        
        m_temp = m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]*m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        m_tempOther = m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]*m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];
        
        const CFreal halfFaceJacob = 0.5 * m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        
        // loop over the states to perturb the states (l)
        m_nbrSolDep = ((*m_flxSolDep)[flxPntIdxThis]).size();
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
          
          m_temp2 = m_temp * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx];
          m_tempOther2 = m_tempOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx];
            
          // add the second part of the discontinuous part of the jacobian (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
          {
            const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
            const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

            // get the divergence of the correction function on this side
            CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
                          
            // add part on this side of face
            m_tempFlux = m_temp2 * divh;
              
            acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
                        
            // add cross-cell part 
            m_tempFlux = m_tempOther2 * divh;   
              
            acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
        }
        
        // loop over the states to perturb the states (l) for LLAV to state part
        for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
        { 
          m_temp2 = m_llavRiemannFluxJacobian[m_pertSide][m_pertSol][m_pertVar][iFlxPnt] * m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
          m_tempOther2 = m_llavRiemannFluxJacobian[m_pertSide][m_pertSol][m_pertVar][iFlxPnt] * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];
          
          // add the LLAV interface part of the jacobian (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
          {
            const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
            const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

            // get the divergence of the correction function on this side
            CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
                          
            // add part on this side of face
            m_tempFlux = m_temp2 * divh;
              
            acc.addValues(jSolIdxThis+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
                        
            // add cross-cell part 
            m_tempFlux = m_tempOther2 * divh;   
              
            acc.addValues(jSolIdxOther+otherSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
        }
        
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        { 
          m_needToAddSolPnt[iSol] = true;
        }
            
        // (i)
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
        {
          const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
          const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

          // get the divergence of the correction function on this side
          CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
            
          const CFreal divh_halfFaceJacob = divh * halfFaceJacob;
          
          // loop over the states to perturb the states (l)
          for (CFuint m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {   
            const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
            
            m_needToAddSolPnt[pertSolIdx] = false;
              
            // add part on this side of face
            m_tempFlux = 0.0;

            // (m)
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
              const CFreal divh_halfFaceJacob_l = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
              const CFreal divh_halfFaceJacob_lOther = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];
              
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacob_l;
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacob_lOther;
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacob_l_dqdu = divh_halfFaceJacob_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                const CFreal divh_halfFaceJacob_l_dqduOther = divh_halfFaceJacob_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];
                
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  m_temp =  m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var];
                    
                  m_tempFlux += m_temp * divh_halfFaceJacob_l_dqdu;
                    
                  m_tempFlux += m_temp * divh_halfFaceJacob_l_dqduOther;
                }
                
                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
                
                const CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                
                // add part of analytical LLAV jacobian
                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacob_l_dqdu;    
                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacob_l_dqduOther;  
              }         
            }
              
            acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
            
            const CFreal divh_halfFaceJacobOther = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];

            // add cross-cell part 
            m_tempFlux = 0.0;

            // (m)
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
              const CFreal divh_halfFaceJacobOther_lThis = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
              const CFreal divh_halfFaceJacobOther_lOther = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];
              
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacobOther_lThis;
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacobOther_lOther;
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacobOther_lThis_dqduThis = divh_halfFaceJacobOther_lThis * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                const CFreal divh_halfFaceJacobOther_lOther_dqduOther = divh_halfFaceJacobOther_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];
              
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  m_temp = m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var];
                    
                  m_tempFlux += m_temp * divh_halfFaceJacobOther_lThis_dqduThis;
                    
                  m_tempFlux += m_temp * divh_halfFaceJacobOther_lOther_dqduOther;
                }
                
                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
                
                const CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                
                // add part of analytical LLAV jacobian
                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacobOther_lThis_dqduThis;    
                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacobOther_lOther_dqduOther;
              }         
            }
            
            acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
          
          // loop over the states to perturb the states (l)
          for (CFuint pertSolIdx = 0; pertSolIdx < m_nbrSolPnts; ++pertSolIdx)
          {
            CFuint dependingKSol = 1000;
              
            if (m_needToAddSolPnt[pertSolIdx])
            {
              // add part on this side of face
              m_tempFlux = 0.0;

              // (m)
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              
                for (CFuint lSol = 0; lSol < m_nbrSolSolDep; ++lSol)
                {
                  const CFuint lSolIdx = (*m_solSolDep)[pertSolIdx][lSol]; 
                
                  if (lSolIdx == kSolIdx)
                  {
                    dependingKSol = kSolIdx;
                    break;
                  }
                }
              }
  
              const CFreal divh_halfFaceJacob_l = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][dependingKSol];
              
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacob_l;
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacob_l_dqdu = divh_halfFaceJacob_l * m_gradientStateJacobian[m_pertSide][dependingKSol][m_pertSide][pertSolIdx][iDim];
                
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  const CFreal divh_halfFaceJacob_l_dqdu_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * divh_halfFaceJacob_l_dqdu;
                    
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * divh_halfFaceJacob_l_dqdu_dudu;
                }
                
                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
                
                CFreal llavPart =  epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                
                llavPart *= divh_halfFaceJacob_l_dqdu;
                
                // add part of analytical LLAV jacobian
                m_tempFlux[m_pertVar] += llavPart;     
              } 
              
              acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
              // get the divergence of the correction function on other side
              divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
            
              const CFreal divh_halfFaceJacobOther = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];

              // add cross-cell part 
              m_tempFlux = 0.0;
              
              const CFreal divh_halfFaceJacobOther_lThis = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][dependingKSol];
                
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacobOther_lThis;
              
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacobOther_lThis_dqduThis = divh_halfFaceJacobOther_lThis * m_gradientStateJacobian[m_pertSide][dependingKSol][m_pertSide][pertSolIdx][iDim];
              
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  const CFreal divh_halfFaceJacobOther_lThis_dqduThis_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * divh_halfFaceJacobOther_lThis_dqduThis;
                    
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * divh_halfFaceJacobOther_lThis_dqduThis_dudu;
                }
                
                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
                
                CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                llavPart *= divh_halfFaceJacobOther_lThis_dqduThis;
                
                // add part of analytical LLAV jacobian
                m_tempFlux[m_pertVar] += llavPart;    
              }         
            
              acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
        }
        
        //// add the cross-element gradient part
        
        // loop over the states to perturb the states (l)
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
  
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
          {
            const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
            
            // add the first and second part of the discontinuous gradient part of the jacobian
            for (CFuint iInfluencedFlx = 0; iInfluencedFlx < m_nbrFlxDep; ++iInfluencedFlx)
            {
              const CFuint iInfluencedFlxIdx = (*m_solFlxDep)[kSolIdxOther][iInfluencedFlx];

              const CFuint dimOther = (*m_flxPntFlxDim)[iInfluencedFlxIdx];
              
              const CFreal lOther = (*m_solPolyValsAtFlxPnts)[iInfluencedFlxIdx][kSolIdxOther];
                
              for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
              {   
                const CFuint jSolIdx = (*m_solSolDep)[kSolIdxOther][jSolPnt];
                  
                m_tempFlux = 0.0;
                  
                // get the divergence of the correction function on this side
                const CFreal divh_lOther = -m_corrFctDiv[jSolIdx][iInfluencedFlxIdx] * lOther; 
                
                /// llav jacob to state part //// actually should go over all kSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iOtherSide][jSolIdx][dimOther] * divh_lOther;
              
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqduOther = divh_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                      
                    const CFreal divh_l_dqduOther_dudu = divh_l_dqduOther * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var];
                      
                    m_tempFlux += m_gradientFluxJacobian[iOtherSide][kSolIdxOther][var][jDim][dimOther] * divh_l_dqduOther_dudu;
                  }
                  
                  CFreal llavPart = divh_l_dqduOther * m_solEpsilons[iOtherSide][kSolIdxOther];
                  
                  llavPart *= m_neighbCellFluxProjVects[iOtherSide][dimOther][kSolIdxOther][jDim];
                  
                  // add part of analytical LLAV Jacobian
                  m_tempFlux[m_pertVar] += llavPart;
                }
                
                acc.addValues(jSolIdx+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
              }
            }
          }
        }
        
        // loop over the states to perturb the states (l)
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
             
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
          {
            const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_solSolDep)[kSolIdxOther][jSolPnt];
                
              m_tempFlux = 0.0;
                
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal lOther = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdxOther];
                
                /// llav jacob to state part //// actually should go over all jSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iOtherSide][jSolIdx][iDim] * lOther;
                  
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                { 
                  const CFreal l_dqduOther = lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal l_dqduOther_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * l_dqduOther;
                      
                    m_tempFlux += m_gradientFluxJacobian[iOtherSide][kSolIdxOther][var][jDim][iDim] * l_dqduOther_dudu;
                  }
                  
                  CFreal llavPart = l_dqduOther * m_solEpsilons[iOtherSide][kSolIdxOther];
                  
                  llavPart *= m_neighbCellFluxProjVects[iOtherSide][iDim][kSolIdxOther][jDim];
                  
                  // add part of analytical LLAV Jacobian
                  m_tempFlux[m_pertVar] += llavPart;
                }
              }
                
              acc.addValues(jSolIdx+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
        }
      }
    }
  }
  
//  if (m_cells[LEFT]->getID() == 1) 
//  {
//      //CFLog(INFO, "ACC: " << acc.getValue(0,4,3,3) << "\n");
//      //acc.printToScreen();
//  }
//  if (m_cells[RIGHT]->getID() == 1) 
//  {
//      //CFLog(INFO, "ACC: " << acc.getValue(4,4,3,3) << "\n");
//      //acc.printToScreen();
//  }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstruction::computeOneJacobDiffFaceTerm(const CFuint side)
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;
  
  CFuint solIdx = 0;
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[m_pertSide])[iSol]->getLocalID());
    }
  }

  //// compute the needed flux jacobians
  
  initJacobianComputation();
  
  computeCellFluxJacobianNum(resFactor);
  
  computeRiemannFluxJacobianNum(resFactor);
  
  computeFluxToGradJacobianNum(resFactor);
  
  computeRiemannFluxToGradJacobianNum(resFactor);
  
  computeGradToStateJacobianAna();
  
  computeGradVarsToStateJacobianNum();
  
  computeEpsToStateJacobianAna();
  
  computeLLAVCellFluxJacobianAna(resFactor);
  
  computeLLAVRiemannFluxJacobianAna(resFactor);
  
  //// add the total jacobians to the system jacobian
  
  // loop over left and right cell to add the discontinuous (cell) part to the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // make sure this is only done once per cell
    if (!m_cellFlags[m_cells[m_pertSide]->getID()]) 
    {
      // term depending on iSide
      const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

      // loop over the states to which to derive (l)
      for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
      {
        // loop over the variables in the state (k)
        for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
        {
          const CFuint nbDepGradVar = m_nbrVarToGradVarDep[m_pertVar];
            
          // add the discontinuous part of the jacobian related to the sol pnt (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
          {
            const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];
            
            m_tempFlux = 0.;

            // (d)
            for (CFuint iDim = 0; iDim < m_dim; ++iDim)
            {
              const CFreal polyCoef = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][m_pertSol]; 
          
              m_tempFlux += m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][iDim] * polyCoef;
            }
            
            acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
          
          // add the discontinuous gradient part of the jacobian (m)
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
          {
            const CFuint kSolIdx = (*m_solSolDep)[m_pertSol][kSolPnt];
          
            // (i)
            for (CFuint jSol = 0; jSol < m_nbrSolSolDep; ++jSol)
            {
              const CFuint jSolIdx = (*m_solSolDep)[kSolIdx][jSol];
              
              m_tempFlux = 0.0;
                
              // (d)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal dl = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdx];
                
                /// llav jacob to state part //// should actually be added for all jSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][m_pertSol][m_pertVar][m_pertSide][jSolIdx][iDim] * dl;
                
                // (b)
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal dl_dqdu = dl * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                    
                  // (p)
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal dl_dqdu_dudu = m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][var] * dl_dqdu;
                      
                    m_tempFlux += m_gradientFluxJacobian[m_pertSide][kSolIdx][var][jDim][iDim] * dl_dqdu_dudu;
                  }
                  
                  CFreal llavPart = m_solEpsilons[m_pertSide][kSolIdx] * m_neighbCellFluxProjVects[m_pertSide][iDim][kSolIdx][jDim];
                  llavPart *= dl_dqdu;
                  //if(m_cells[0]->getID()==1) CFLog(INFO, "pertSol: " << m_pertSol << ", pertVar: " << m_pertVar << ", llavJC: " << llavPart << "\n");
                  // add part of analytical LLAV jacobian
                  m_tempFlux[m_pertVar] += llavPart;
                }
              }
            
              acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
          
          // add the discontinuous part of the jacobian related to the flx pnt (f)
          for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
          {
            const CFuint flxIdx = (*m_solFlxDep)[m_pertSol][iFlxPnt];
            
            // (df)
            const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
            
            m_temp = m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][dim] * (*m_solPolyValsAtFlxPnts)[flxIdx][m_pertSol];
            
            // add the second part of the discontinuous part of the jacobian (i)
            m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[flxIdx][jSolPnt];

              // get the divergence of the correction function
              const CFreal divh = m_corrFctDiv[jSolIdx][flxIdx];
                           
              m_tempFlux = -m_temp * divh;
            
              acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
          
          // add the second part of the discontinuous gradient part of the jacobian (m)
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
          {
            const CFuint kSolIdx = (*m_solSolDep)[m_pertSol][kSolPnt];
              
            for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
            {   
              const CFuint flxIdx = (*m_solFlxDep)[kSolIdx][iFlxPnt];
              
              // (df)
              const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
              
              const CFreal l = (*m_solPolyValsAtFlxPnts)[flxIdx][kSolIdx];            
              
              // (i)
              m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
              for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
              {
                const CFuint jSolIdx = (*m_flxSolDep)[flxIdx][jSolPnt];
              
                m_tempFlux = 0.0;

                const CFreal divh_l = -m_corrFctDiv[jSolIdx][flxIdx] * l;
                
                /// llav jacob to state part //// actually should loop over all kSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][m_pertSol][m_pertVar][m_pertSide][jSolIdx][dim] * divh_l;
              
                // (b)
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqdu = divh_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                  
                  // (p)
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal divh_l_dqdu_dudu = divh_l_dqdu * m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][var];             
                      
                    m_tempFlux += divh_l_dqdu_dudu * m_gradientFluxJacobian[m_pertSide][kSolIdx][var][jDim][dim];
                  }
                  
                  CFreal llavPart = divh_l_dqdu * m_solEpsilons[m_pertSide][kSolIdx];
                  
                  llavPart *= m_neighbCellFluxProjVects[m_pertSide][dim][kSolIdx][jDim];
                  //if(m_cells[0]->getID()==1) CFLog(INFO, "pertSol: " << m_pertSol << ", pertVar: " << m_pertVar << ", llavJF: " << llavPart << "\n");
                  // add part of analytical LLAV Jacobian
                  m_tempFlux[m_pertVar] += llavPart;
                }
                
                acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
              }
            }
          }
        }
      }
    }
  }
  
  // loop over left and right cell to add the riemann flux (face) part to the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // term depending on iSide
    const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

    // term depending on iOtherSide
    const CFuint otherSideTerm = iOtherSide*m_nbrSolPnts;
    
    // loop over the variables in the state (k)
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    { 
      const CFuint nbDepGradVar = m_nbrVarToGradVarDep[m_pertVar];
        
      // loop over face flx pnts (f)
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // local flux point indices in the left and right cell
        const CFuint flxPntIdxThis = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
        const CFuint flxPntIdxOther = (*m_faceFlxPntConnPerOrient)[m_orient][iOtherSide][iFlxPnt];
        
        m_temp = m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]*m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        m_tempOther = m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]*m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];
        
        const CFreal halfFaceJacob = 0.5 * m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        
        // loop over the states to perturb the states (l)
        m_nbrSolDep = ((*m_flxSolDep)[flxPntIdxThis]).size();
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
          
          m_temp2 = m_temp * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx];
          m_tempOther2 = m_tempOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx];
            
          // add the second part of the discontinuous part of the jacobian (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
          {
            const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
            const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

            // get the divergence of the correction function on this side
            CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
                          
            // add part on this side of face
            m_tempFlux = m_temp2 * divh;
              
            acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
                        
            // add cross-cell part 
            m_tempFlux = m_tempOther2 * divh;   
              
            acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
        }
        
        // loop over the states to perturb the states (l) for LLAV to state part
        for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
        { 
          m_temp2 = m_llavRiemannFluxJacobian[m_pertSide][m_pertSol][m_pertVar][iFlxPnt] * m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
          m_tempOther2 = m_llavRiemannFluxJacobian[m_pertSide][m_pertSol][m_pertVar][iFlxPnt] * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];
          
          // add the LLAV interface part of the jacobian (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
          {
            const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
            const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

            // get the divergence of the correction function on this side
            CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
                          
            // add part on this side of face
            m_tempFlux = m_temp2 * divh;
              
            acc.addValues(jSolIdxThis+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
                        
            // add cross-cell part 
            m_tempFlux = m_tempOther2 * divh;   
              
            acc.addValues(jSolIdxOther+otherSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
        }
        
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        { 
          m_needToAddSolPnt[iSol] = true;
        }
            
        // (i)
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
        {
          const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
          const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

          // get the divergence of the correction function on this side
          CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
            
          const CFreal divh_halfFaceJacob = divh * halfFaceJacob;
          
          // loop over the states to perturb the states (l)
          for (CFuint m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {   
            const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
            
            m_needToAddSolPnt[pertSolIdx] = false;
              
            // add part on this side of face
            m_tempFlux = 0.0;

            // (m)
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
              const CFreal divh_halfFaceJacob_l = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
              const CFreal divh_halfFaceJacob_lOther = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];
              
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacob_l;
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacob_lOther;
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacob_l_dqdu = divh_halfFaceJacob_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                const CFreal divh_halfFaceJacob_l_dqduOther = divh_halfFaceJacob_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];
                
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  m_temp =  m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var];
                    
                  m_tempFlux += m_temp * divh_halfFaceJacob_l_dqdu;
                    
                  m_tempFlux += m_temp * divh_halfFaceJacob_l_dqduOther;
                }
                
                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
                
                const CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                
                // add part of analytical LLAV jacobian
                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacob_l_dqdu;    
                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacob_l_dqduOther;  
              }         
            }
              
            acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
            
            const CFreal divh_halfFaceJacobOther = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];

            // add cross-cell part 
            m_tempFlux = 0.0;

            // (m)
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
              const CFreal divh_halfFaceJacobOther_lThis = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
              const CFreal divh_halfFaceJacobOther_lOther = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];
              
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacobOther_lThis;
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacobOther_lOther;
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacobOther_lThis_dqduThis = divh_halfFaceJacobOther_lThis * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                const CFreal divh_halfFaceJacobOther_lOther_dqduOther = divh_halfFaceJacobOther_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];
              
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  m_temp = m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var];
                    
                  m_tempFlux += m_temp * divh_halfFaceJacobOther_lThis_dqduThis;
                    
                  m_tempFlux += m_temp * divh_halfFaceJacobOther_lOther_dqduOther;
                }
                
                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
                
                const CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                
                // add part of analytical LLAV jacobian
                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacobOther_lThis_dqduThis;    
                m_tempFlux[m_pertVar] += llavPart * divh_halfFaceJacobOther_lOther_dqduOther;
              }         
            }
            
            acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
          
          // loop over the states to perturb the states (l)
          for (CFuint pertSolIdx = 0; pertSolIdx < m_nbrSolPnts; ++pertSolIdx)
          {
            CFuint dependingKSol = 1000;
              
            if (m_needToAddSolPnt[pertSolIdx])
            {
              // add part on this side of face
              m_tempFlux = 0.0;

              // (m)
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              
                for (CFuint lSol = 0; lSol < m_nbrSolSolDep; ++lSol)
                {
                  const CFuint lSolIdx = (*m_solSolDep)[pertSolIdx][lSol]; 
                
                  if (lSolIdx == kSolIdx)
                  {
                    dependingKSol = kSolIdx;
                    break;
                  }
                }
              }
  
              const CFreal divh_halfFaceJacob_l = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][dependingKSol];
              
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacob_l;
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacob_l_dqdu = divh_halfFaceJacob_l * m_gradientStateJacobian[m_pertSide][dependingKSol][m_pertSide][pertSolIdx][iDim];
                
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  const CFreal divh_halfFaceJacob_l_dqdu_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * divh_halfFaceJacob_l_dqdu;
                    
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * divh_halfFaceJacob_l_dqdu_dudu;
                }
                
                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
                
                CFreal llavPart =  epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                
                llavPart *= divh_halfFaceJacob_l_dqdu;
                
                // add part of analytical LLAV jacobian
                m_tempFlux[m_pertVar] += llavPart;     
              } 
              
              acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
              // get the divergence of the correction function on other side
              divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
            
              const CFreal divh_halfFaceJacobOther = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];

              // add cross-cell part 
              m_tempFlux = 0.0;
              
              const CFreal divh_halfFaceJacobOther_lThis = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][dependingKSol];
                
              /// llav jacob to state part
              //m_tempFlux += m_llavRiemannFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iFlxPnt] * divh_halfFaceJacobOther_lThis;
              
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacobOther_lThis_dqduThis = divh_halfFaceJacobOther_lThis * m_gradientStateJacobian[m_pertSide][dependingKSol][m_pertSide][pertSolIdx][iDim];
              
                // (p)
                for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                {
                  const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                  
                  const CFreal divh_halfFaceJacobOther_lThis_dqduThis_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * divh_halfFaceJacobOther_lThis_dqduThis;
                    
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][var][iDim] * divh_halfFaceJacobOther_lThis_dqduThis_dudu;
                }
                
                const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
                
                CFreal llavPart = epsilon * m_unitNormalFlxPnts[iFlxPnt][iDim];
                llavPart *= divh_halfFaceJacobOther_lThis_dqduThis;
                
                // add part of analytical LLAV jacobian
                m_tempFlux[m_pertVar] += llavPart;    
              }         
            
              acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
        }
        
        //// add the cross-element gradient part
        
        // loop over the states to perturb the states (l)
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
  
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
          {
            const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
            
            // add the first and second part of the discontinuous gradient part of the jacobian
            for (CFuint iInfluencedFlx = 0; iInfluencedFlx < m_nbrFlxDep; ++iInfluencedFlx)
            {
              const CFuint iInfluencedFlxIdx = (*m_solFlxDep)[kSolIdxOther][iInfluencedFlx];

              const CFuint dimOther = (*m_flxPntFlxDim)[iInfluencedFlxIdx];
              
              const CFreal lOther = (*m_solPolyValsAtFlxPnts)[iInfluencedFlxIdx][kSolIdxOther];
                
              for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
              {   
                const CFuint jSolIdx = (*m_solSolDep)[kSolIdxOther][jSolPnt];
                  
                m_tempFlux = 0.0;
                  
                // get the divergence of the correction function on this side
                const CFreal divh_lOther = -m_corrFctDiv[jSolIdx][iInfluencedFlxIdx] * lOther; 
                
                /// llav jacob to state part //// actually should go over all kSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iOtherSide][jSolIdx][dimOther] * divh_lOther;
              
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqduOther = divh_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                      
                    const CFreal divh_l_dqduOther_dudu = divh_l_dqduOther * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var];
                      
                    m_tempFlux += m_gradientFluxJacobian[iOtherSide][kSolIdxOther][var][jDim][dimOther] * divh_l_dqduOther_dudu;
                  }
                  
                  CFreal llavPart = divh_l_dqduOther * m_solEpsilons[iOtherSide][kSolIdxOther];
                  
                  llavPart *= m_neighbCellFluxProjVects[iOtherSide][dimOther][kSolIdxOther][jDim];
                  
                  // add part of analytical LLAV Jacobian
                  m_tempFlux[m_pertVar] += llavPart;
                }
                
                acc.addValues(jSolIdx+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
              }
            }
          }
        }
        
        // loop over the states to perturb the states (l)
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
             
          for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
          {
            const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_solSolDep)[kSolIdxOther][jSolPnt];
                
              m_tempFlux = 0.0;
                
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal lOther = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdxOther];
                
                /// llav jacob to state part //// actually should go over all jSol
                m_tempFlux += m_llavFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iOtherSide][jSolIdx][iDim] * lOther;
                  
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                { 
                  const CFreal l_dqduOther = lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal l_dqduOther_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * l_dqduOther;
                      
                    m_tempFlux += m_gradientFluxJacobian[iOtherSide][kSolIdxOther][var][jDim][iDim] * l_dqduOther_dudu;
                  }
                  
                  CFreal llavPart = l_dqduOther * m_solEpsilons[iOtherSide][kSolIdxOther];
                  
                  llavPart *= m_neighbCellFluxProjVects[iOtherSide][iDim][kSolIdxOther][jDim];
                  
                  // add part of analytical LLAV Jacobian
                  m_tempFlux[m_pertVar] += llavPart;
                }
              }
                
              acc.addValues(jSolIdx+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
        }
      }
    }
  }
  
//  if (m_cells[LEFT]->getID() == 1) 
//  {
//      //CFLog(INFO, "ACC: " << acc.getValue(0,4,3,3) << "\n");
//      //acc.printToScreen();
//  }
//  if (m_cells[RIGHT]->getID() == 1) 
//  {
//      //CFLog(INFO, "ACC: " << acc.getValue(4,4,3,3) << "\n");
//      //acc.printToScreen();
//  }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;

  // setup parent class
  ConvDiffLLAVJacobFluxReconstruction::setup();
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type orderLLAVCorr = static_cast<CFPolyOrder::Type> (m_LLAVCorrFctOrder);
  
  const CFuint totNbFlxPnts = m_flxPntsLocalCoords->size();
  
  m_corrFctLLAV.resize(m_nbrSolPnts);
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_corrFctLLAV[iSol].resize(totNbFlxPnts);
  }
  
  // compute the divergence of the correction function
  m_corrFctComputer->computeCorrectionFunction(orderLLAVCorr,m_LLAVCorrFctFactor,frLocalData[0],m_corrFctLLAV);
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  // get the number of cells in the mesh
  const CFuint nbrCells = (*elemType)[0].getEndIdx();
  
  const CFuint nbStates = nbrCells*m_nbrSolPnts;
  
  m_nodeEpsilons.resize(nbStates);
  m_nbNodeNeighbors.resize(0);
  
  m_alphaValues.resize(2);
  m_alphaValues[LEFT].resize(m_nbrSolPnts);
  m_alphaValues[RIGHT].resize(m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  ConvDiffLLAVJacobFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
