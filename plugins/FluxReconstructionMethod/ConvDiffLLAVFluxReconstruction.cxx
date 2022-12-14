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

#include "FluxReconstructionMethod/ConvDiffLLAVFluxReconstruction.hh"
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

MethodCommandProvider< ConvDiffLLAVFluxReconstruction,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
convDiffLLAVRHSFluxReconstructionProvider("ConvDiffLLAVRHS");

//////////////////////////////////////////////////////////////////////////////
  
ConvDiffLLAVFluxReconstruction::ConvDiffLLAVFluxReconstruction(const std::string& name) :
  DiffRHSFluxReconstruction(name),
  m_updateVarSet(CFNULL),
  m_order(),
  m_transformationMatrix(),
  m_statesPMinOne(),
  m_epsilon(),
  m_solEpsilons(),
  m_epsilonLR(),
  m_epsilon0(),
  m_s0(),
  m_s(),
  m_kappa(),
  m_showrate(),
  m_peclet(),
  m_nodeEpsilons(),
  m_nbNodeNeighbors(),
  m_cellEpsilons(),
  m_cellNodes(),
  m_nbrCornerNodes(),
  m_faceNodes(),
  m_flagComputeNbNghb(),
  m_nodePolyValsAtFlxPnts(),
  m_nodePolyValsAtSolPnts(),
  m_cellNodesConn(CFNULL),
  m_elemIdx(),
  m_facesCell(),
  m_jacob(),
  m_freezeLimiterRes(),
  m_freezeLimiterIter(),
  m_useMax(),
  m_totalEps(),
  m_totalEpsGlobal(),
  m_Smax(),
  m_SmaxGlobal(),
  m_subcellRes(),
  m_unitNormalFlxPnts2(),
  socket_artVisc("artVisc"),
  socket_monPhysVar("monPhysVar"),
  socket_smoothness("smoothness"),
  socket_solPntNormals("solPntNormals"),
  socket_flxPntNormals("flxPntNormals"),
  socket_cellVolumes("cellVolumes"),
  socket_volumes("volumes"),
  socket_wallDistance("wallDistance"),
  m_epsBackUp(),
  m_tempSolPntVec(),
  m_tempSolPntVec2(),
  m_cellGradsAV(),
  m_cellGradFlxPntAV(),
  m_avgGradAV(),
  m_cellFlags(),
  m_cellBuilders(),
  m_cellBuilder(CFNULL),
  m_isFaceOnBoundary(),
  m_nghbrCellSide(),
  m_currCellSide(),
  m_faceOrients(),
  m_faceBCIdx(),
  m_isFaceOnBoundaryCell(CFNULL),
  m_nghbrCellSideCell(CFNULL),
  m_currCellSideCell(CFNULL),
  m_faceOrientsCell(CFNULL),
  m_faceBCIdxCell(CFNULL),
  m_unpertCellDiffRes(),
  m_gradUpdates(),
  m_projectedCorrL(),
  m_projectedCorrR(),
  m_contFlxWoLLAV(),
  m_updateToSolutionVecTrans(CFNULL)
  {
    addConfigOptionsTo(this);
    
    m_kappa = 1.0;
    setParameter( "Kappa", &m_kappa);
    
    m_showrate = 1;
    setParameter( "ShowRate", &m_showrate);
    
    m_peclet = 2.0;
    setParameter( "Peclet", &m_peclet);
    
    m_s0 = 0.0;
    setParameter( "S0", &m_s0);
    
    m_dampingCoeff = 1.0;
    setParameter( "DampingCoeff", &m_dampingCoeff);
    
    m_freezeLimiterRes = -20.;
    setParameter( "FreezeLimiterRes", &m_freezeLimiterRes);
  
    m_freezeLimiterIter = MathTools::MathConsts::CFuintMax();
    setParameter( "FreezeLimiterIter", &m_freezeLimiterIter);
    
    m_monitoredVar = 0;
    setParameter( "MonitoredVar", &m_monitoredVar);
    
    m_addUpdCoeff = true;
    setParameter( "AddUpdateCoeff", &m_addUpdCoeff);
    
    m_monitoredPhysVar = MathTools::MathConsts::CFuintMax();
    setParameter( "MonitoredPhysVar", &m_monitoredPhysVar);
    
    m_printLLAV = true;
    setParameter( "printLLAV", &m_printLLAV);
    
    m_LLAVBCZero = true;
    setParameter( "LLAVBCZero", &m_LLAVBCZero);
    
    m_wallCutOff = 0.0;
    setParameter( "WallCutOffDistance", &m_wallCutOff);
    
    m_useWallCutOff = false;
    setParameter( "UseWallCutOff", &m_useWallCutOff);
  }

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  DiffRHSFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal,Config::DynamicOption<> >("Kappa","Kappa factor of artificial viscosity.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("ShowRate","Showrate of LLAV information.");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("Peclet","Peclet number to be used for artificial viscosity.");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("S0","Reference smoothness factor, will be multiplied by -log(P).");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("DampingCoeff","Damping coefficient for reculculation of eps (0<coeff<1).");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("FreezeLimiterRes","Residual after which to freeze the residual.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("FreezeLimiterIter","Iteration after which to freeze the residual.");
    
  options.addConfigOption< CFuint,Config::DynamicOption<> >("MonitoredVar","Index of the monitored var for positivity preservation.");
  
  options.addConfigOption< bool >("AddUpdateCoeff","Boolean telling whether the update coefficient based on the artificial flux is added.");
  
  options.addConfigOption< bool >("printLLAV","Boolean telling whether to print the output of LLAV (default true).");
  
  options.addConfigOption< bool >("LLAVBCZero","Boolean telling whether to set LLAV to zero on the bnd.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("MonitoredPhysVar","Index of the monitored physical var for positivity preservation, if not specified MonitoredVar is used instead.");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("WallCutOffDistance","Distance from wall at which to cut off LLAV.");
  
  options.addConfigOption< bool >("UseWallCutOff","Boolean telling whether to use wall distance cut off of LLAV.");
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  ConvDiffLLAVFluxReconstruction::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = DiffRHSFluxReconstruction::providesSockets();
  result.push_back(&socket_artVisc);
  result.push_back(&socket_monPhysVar);
  result.push_back(&socket_smoothness);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
ConvDiffLLAVFluxReconstruction::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = DiffRHSFluxReconstruction::needsSockets();
  result.push_back(&socket_solPntNormals);
  result.push_back(&socket_flxPntNormals);
  result.push_back(&socket_cellVolumes);
  result.push_back(&socket_volumes);
  if (m_useWallCutOff) result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ConvDiffLLAVFluxReconstruction::execute()\n");
  
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

  //// Loop over faces to add the face part to the gradients
  
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

      // release the cells
      m_cellBuilders[LEFT ]->releaseGE();
      m_cellBuilders[RIGHT]->releaseGE();
      
      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }

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
      
      // get the nodes in this cell
      m_cellNodes  = m_cell->getNodes();

//       // if the states in the cell are parallel updatable, compute the resUpdates (-divFC)
//       if ((*m_cellStates)[0]->isParUpdatable())
//       {
	// compute the states projected on order P-1
	computeProjStates(m_statesPMinOne);

	// compute the artificial viscosity
	computeEpsilon();

	// store epsilon
	storeEpsilon();
        
        //if (m_cell->getID() == 0) CFLog(INFO, "eps0: " << m_epsilon0 << "\n");
        //if (m_cell->getID() == 0) CFLog(INFO, "s: " << m_s << "\n");
        //if (m_cell->getID() == 0) CFLog(INFO, "eps: " << m_epsilon << ", s: " << m_s << ", s0: " << m_s0 << ", k: " << m_kappa << "\n");
//       } 

      // add the cell part to the gradients
      computeGradients();

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
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
      m_cellVolume[LEFT] = cellVolumes[m_cells[LEFT]->getID()];
      m_cellVolume[RIGHT] = cellVolumes[m_cells[RIGHT]->getID()];
      
      cf_assert(m_cellVolume[LEFT] > 0.0);
      cf_assert(m_cellVolume[RIGHT] > 0.0);
      
	// build the neighbouring cells
        const CFuint cellIDL = m_face->getNeighborGeo(LEFT)->getID();
        geoDataCBL.idx = cellIDL;
        m_cells[LEFT] = m_cellBuilders[LEFT ]->buildGE();
        const CFuint cellIDR = m_face->getNeighborGeo(RIGHT)->getID();
        geoDataCBR.idx = cellIDR;
        m_cells[RIGHT] = m_cellBuilders[RIGHT]->buildGE();

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
	computeCorrection(LEFT, m_divContFlx);

	// update RHS
	updateRHS();
	
	// compute the correction for the right neighbour
	computeCorrection(RIGHT, m_divContFlx);
	
	// update RHS
	updateRHS();
        
        // compute needed cell contributions: what used to be cell loop is incorporated here!!
        if (!m_cellFlags[cellIDL])
        {
          computeUnpertCellDiffResiduals(LEFT);

	  // update RHS
	  // get the datahandle of the rhs
          DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

          // get residual factor
          const CFreal resFactor = getMethodData().getResFactor();

          // update rhs
          for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
          {
            CFuint resID = m_nbrEqs*( (*m_states[LEFT])[iState]->getLocalID() );
            for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
            {
              rhs[resID+iVar] += resFactor*m_unpertCellDiffRes[LEFT][m_nbrEqs*iState+iVar];
            }
          }
        }
        if (!m_cellFlags[cellIDR])
        {
          computeUnpertCellDiffResiduals(RIGHT);

	  // update RHS
	  // get the datahandle of the rhs
          DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

          // get residual factor
          const CFreal resFactor = getMethodData().getResFactor();

          // update rhs
          for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
          {
            CFuint resID = m_nbrEqs*( (*m_states[RIGHT])[iState]->getLocalID() );
            for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
            {
              rhs[resID+iVar] += resFactor*m_unpertCellDiffRes[RIGHT][m_nbrEqs*iState+iVar];
            }
          }
        }

        // release the cells
        m_cellBuilders[LEFT ]->releaseGE();
        m_cellBuilders[RIGHT]->releaseGE();
        
        m_cellFlags[cellIDL] = true;
        m_cellFlags[cellIDR] = true;
      
      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::setFaceDataForGradients(CFuint faceID)
{   
  // get the face flux point normals
  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
  
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_faceJacobVecs[iFlxPnt][iDim] = flxPntNormals[m_face->getID()*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
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
    m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
  }

  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {     
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];
    
    // reset states in flx pnt
    *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdxL = (*m_flxSolDep)[flxPntIdxL][iSol];
      const CFuint solIdxR = (*m_flxSolDep)[flxPntIdxR][iSol];
 
      // add the contributions of the current sol pnt
      *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*(*((*(m_states[LEFT]))[solIdxL]));
      *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*(*((*(m_states[RIGHT]))[solIdxR]));
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::setFaceData(CFuint faceID)
{   
  // get the face flux point normals
  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
  
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_faceJacobVecs[iFlxPnt][iDim] = flxPntNormals[m_face->getID()*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
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
    m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
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
  
  // get the AV properties
  m_faceNodes = m_face->getNodes();
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {   
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt]; 
    
    // reset the states in the flx pnts
    m_epsilonLR[LEFT][iFlxPnt] = 0.0;
    m_epsilonLR[RIGHT][iFlxPnt] = 0.0;
    
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_cellNodes = m_cells[iSide]->getNodes();
      CFuint flxIdx;
      iSide == LEFT ? flxIdx = flxPntIdxL : flxIdx = flxPntIdxR;

      for (CFuint iNode = 0; iNode < m_faceNodes->size(); ++iNode)
      {
	for (CFuint iNodeCell = 0; iNodeCell < m_nbrCornerNodes; ++iNodeCell)
        {
            
            //if(m_cells[iSide]->getID() == 0) CFLog(INFO,"faceNodeID: " << (*m_faceNodes)[iNode]->getLocalID() << ", cellNodeID: " << (*m_cellNodes)[iNodeCell]->getLocalID() << "\n");
                    
	  if ((*m_faceNodes)[iNode]->getLocalID() == (*m_cellNodes)[iNodeCell]->getLocalID())
	  {
	    // get node local index
            const CFuint nodeIdx = (*m_cellNodesConn)(m_cells[iSide]->getID(),iNodeCell);
	    
            m_epsilonLR[iSide][iFlxPnt] += m_nodePolyValsAtFlxPnts[flxIdx][iNodeCell]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
            
            //if(m_cells[iSide]->getID() == 0) CFLog(INFO,"node eps: " << m_nodeEpsilons[nodeIdx] <<", polyVal: " << m_nodePolyValsAtFlxPnts[flxIdx][iNodeCell] << ", nbNeighb: " << m_nbNodeNeighbors[nodeIdx] << "\n");
	  }
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeFlxPntStatesAndGrads()
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
      *(m_cellGradFlxPntAV[LEFT][iFlxPnt][iVar]) = 0.0;
      *(m_cellGradFlxPntAV[RIGHT][iFlxPnt][iVar]) = 0.0;
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
        *(m_cellGradFlxPntAV[LEFT][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*((*(m_cellGradsAV[LEFT][solIdxL]))[iVar]);
	*(m_cellGradFlxPntAV[RIGHT][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*((*(m_cellGradsAV[RIGHT][solIdxR]))[iVar]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeInterfaceFlxCorrection()
{
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    // compute the average sol and grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
      
      *(m_avgGradAV[iVar]) = (*(m_cellGradFlxPntAV[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPntAV[RIGHT][iFlxPnt][iVar]))/2.0;
             
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
    }
    
    prepareFluxComputation();
     
    // compute diffusive flux
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFlux[iFlxPnt]);
    
    // compute the convective riemann flux
    m_flxPntRiemannFlux[iFlxPnt] -= m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),
									*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]),
									m_unitNormalFlxPnts[iFlxPnt]);
    
    // compute artificial part
    // get epsilon
    const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
    
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_flxPntRiemannFlux[iFlxPnt][iVar] += epsilon*((*(m_avgGradAV[iVar]))[iDim])*m_unitNormalFlxPnts[iFlxPnt][iDim];
        
        //if (m_cells[LEFT]->getID() == 0 || m_cells[RIGHT]->getID() == 0) printf("second flx: %d, var: %d, dim: %d, eps: %e, grad: %e, n: %e\n",iFlxPnt,iVar,iDim,epsilon,((*(m_avgGradAV[iVar]))[iDim]),m_unitNormalFlxPnts[iFlxPnt][iDim]); 
      }
    }
     
    // compute FI in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  
  // add convective part
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    waveSpeedUpd[iSide] = 0.0;

    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFreal jacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceIntegrationCoefs)[iFlx];

      // transform update states to physical data to calculate eigenvalues
      m_updateVarSet->computePhysicalData(*(m_cellStatesFlxPnt[iSide][iFlx]), m_pData);

      waveSpeedUpd[iSide] += jacobXIntCoef * m_updateVarSet->getMaxAbsEigenValue(m_pData,m_unitNormalFlxPnts[iFlx]);
    }
  }
  
  // add artificial part if needed
  if(m_addUpdCoeff)
  {
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
      {
        const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                         m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                         (*m_faceIntegrationCoefs)[iFlx]*
                                         m_cflConvDiffRatio;
        
        const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlx]+m_epsilonLR[RIGHT][iFlx]);

        // transform update states to physical data to calculate eigenvalues
        waveSpeedUpd[iSide] += epsilon*jacobXJacobXIntCoef/m_cellVolume[iSide];
      }
    } 
  }
  
  //////Diffusive part is added in NS!!!!
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeUnpertCellDiffResiduals(const CFuint side)
{
  // get the sol pnt normals
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
  
  // store the sol pnt normals
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    const CFuint solID = (*(m_states[side]))[iState]->getLocalID();
      
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      for (CFuint jDim = 0; jDim < m_dim; ++jDim)
      {
        m_cellFluxProjVects[iDim][iState][jDim] = solPntNormals[solID*(m_dim+m_ndimplus)*m_dim+iDim*m_dim+jDim];
      }
    }
  }
  
  m_cellNodes = m_cells[side]->getNodes();
  
  // get datahandle
  DataHandle< CFreal > artVisc = socket_artVisc.getDataHandle();
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {   
    // reset the states in the flx pnts
    m_solEpsilons[side][iSol] = 0.0;

    // loop over the sol pnts to compute the states and grads in the flx pnts
    for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
    {
      // get node local index
      const CFuint nodeIdx = (*m_cellNodes)[iNode]->getLocalID();
      
      m_solEpsilons[side][iSol] += m_nodePolyValsAtSolPnts[iSol][iNode]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
    }
        
    artVisc[(((*(m_states[side]))[iSol]))->getLocalID()] = m_solEpsilons[side][iSol];
    
    //if (m_cells[side]->getID() == 0) CFLog(INFO, "eps: " << m_solEpsilons[side][iSol] << "\n");
  }
  
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
      *(m_tempGrad[iVar]) = (*(m_cellGrads[side][iSolPnt]))[iVar];
      *(m_avgGradAV[iVar]) = (*(m_cellGradsAV[side][iSolPnt]))[iVar];
    }
    
    m_updateVarSet->computePhysicalData(*((*(m_states[side]))[iSolPnt]), m_pData); 

    m_avgSol = *((*(m_states[side]))[iSolPnt]->getData());

    prepareFluxComputation();

    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      // add diffusive part 
      computeFlux(m_avgSol,m_tempGrad,m_cellFluxProjVects[iDim][iSolPnt],0,m_contFlxWoLLAV[iSolPnt][iDim]);
      
      // add convective part
      m_contFlxWoLLAV[iSolPnt][iDim] -= m_updateVarSet->getFlux()(m_pData,m_cellFluxProjVects[iDim][iSolPnt]);
      
      m_contFlx[iSolPnt][iDim] = m_contFlxWoLLAV[iSolPnt][iDim];
      
      // add artificial part
      for (CFuint iDim2 = 0; iDim2 < m_dim; ++iDim2)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_contFlx[iSolPnt][iDim][iVar] += m_solEpsilons[side][iSolPnt]*((*m_avgGradAV[iVar])[iDim2])*m_cellFluxProjVects[iDim][iSolPnt][iDim2];
          
//       if (m_cells[side]->getID() == 0) CFLog(INFO,"first iSol: " << iSolPnt << ", iDir: " << iDim << ", iEq: " << iVar << 
//               ", gradAV: " << ((*m_avgGradAV[iVar])[iDim2]) << ", eps: " << m_solEpsilons[side][iSolPnt] << ", n: " << m_cellFluxProjVects[iDim][iSolPnt][iDim2] << "\n");  
        }
      }
    }

//    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
//    {
//      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];
//      const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
//
//      m_extrapolatedFluxes[flxIdx] += (*m_solPolyValsAtFlxPnts)[flxIdx][iSolPnt]*(m_contFlx[iSolPnt][dim]);
//    }
  }
  
  // add the contribution of the faces
  const CFuint nbrFaces = m_cells[side]->nbNeighborGeos();

  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if (!((*m_isFaceOnBoundaryCell)[iFace]) || m_LLAVBCZero)
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];

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

        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
          const CFuint dim = (*m_flxPntFlxDim)[currFlxIdx];

           m_extrapolatedFluxes[currFlxIdx] += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(m_contFlxWoLLAV[solIdx][dim]);
        }
      } 
    }
  }

  // Loop over solution pnts to calculate the divergence of the discontinuous flux
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // reset the divergence of FC
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_unpertCellDiffRes[side][m_nbrEqs*iSolPnt+iEq] = 0.0;
    }

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
          m_unpertCellDiffRes[side][m_nbrEqs*iSolPnt+iEq] += polyCoef*(m_contFlx[jSolIdx][iDir+m_ndimplus][iEq]);
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
        m_unpertCellDiffRes[side][m_nbrEqs*iSolPnt+iVar] += -m_extrapolatedFluxes[flxIdx][iVar] * divh; 
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeGradientFaceCorrections()
{ 
  // get the face flux point normals
  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
  
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_faceJacobVecs[iFlxPnt][iDim] = flxPntNormals[m_face->getID()*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
    }
  }

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
    
    // reset states in flx pnt
    *(m_cellStatesFlxPnt[LEFT][iFlx]) = 0.0;
    *(m_cellStatesFlxPnt[RIGHT][iFlx]) = 0.0;

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdxL = (*m_flxSolDep)[flxIdxL][iSol];
      const CFuint solIdxR = (*m_flxSolDep)[flxIdxR][iSol];
 
      // add the contributions of the current sol pnt
      *(m_cellStatesFlxPnt[LEFT][iFlx]) += (*m_solPolyValsAtFlxPnts)[flxIdxL][solIdxL]*(*((*(m_states[LEFT]))[solIdxL]));
      *(m_cellStatesFlxPnt[RIGHT][iFlx]) += (*m_solPolyValsAtFlxPnts)[flxIdxR][solIdxR]*(*((*(m_states[RIGHT]))[solIdxR]));
    }

    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFreal avgSol = ((*m_cellStatesFlxPnt[LEFT][iFlx])[iEq]+(*m_cellStatesFlxPnt[RIGHT][iFlx])[iEq])/2.0;
      m_projectedCorrL = (avgSol-(*m_cellStatesFlxPnt[LEFT][iFlx])[iEq])*(*m_faceMappedCoordDir)[m_orient][LEFT]*m_faceJacobVecs[iFlx];
      m_projectedCorrR = (avgSol-(*m_cellStatesFlxPnt[RIGHT][iFlx])[iEq])*(*m_faceMappedCoordDir)[m_orient][RIGHT]*m_faceJacobVecs[iFlx];

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
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeGradients()
{
  // get the sol pnt normals
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();

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
    const CFuint solID = (*m_cellStates)[iSolPnt]->getLocalID();
    
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // Loop over gradient directions
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        for (CFuint jDir = 0; jDir < m_dim; ++jDir)
        {
	  // project the state on a normal and reuse a RealVector variable of the class to store
	  m_projectedCorrL[jDir] = ((*(*m_cellStates)[iSolPnt])[iEq]) * solPntNormals[solID*(m_dim+m_ndimplus)*m_dim+(iDir+m_ndimplus)*m_dim+jDir];
        }
	
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

//  // get jacobian determinants at solution points
//  m_jacobDet = m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

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
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeProjStates(std::vector< RealVector >& projStates)
{
  cf_assert(m_nbrSolPnts == projStates.size());
  
  if (m_order != 1)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        m_tempSolPntVec[iSol] = (*((*m_cellStates)[iSol]))[iEq];
      }

      m_tempSolPntVec2 = m_transformationMatrix*m_tempSolPntVec;

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        projStates[iSol][iEq] = m_tempSolPntVec2[iSol];
      }
    }
  }
  else
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      CFreal stateSum = 0.0;
      
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        stateSum += (*((*m_cellStates)[iSol]))[iEq];
      }

      stateSum /= m_nbrSolPnts;

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        projStates[iSol][iEq] = stateSum;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeProjStates(std::vector< RealVector >& projStates, const CFuint side)
{
  cf_assert(m_nbrSolPnts == projStates.size());
    
  if (m_order != 1)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        m_tempSolPntVec[iSol] = (*((*m_states[side])[iSol]))[iEq];
      }

      m_tempSolPntVec = m_transformationMatrix*m_tempSolPntVec;

      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        projStates[iSol][iEq] = m_tempSolPntVec[iSol];
      }
    }
  }
  else
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      CFreal stateSum = 0.0;
      
      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        stateSum += (*((*m_states[side])[iSol]))[iEq];
      }

      stateSum /= projStates.size();

      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        projStates[iSol][iEq] = stateSum;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeEpsilon()
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
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvDiffLLAVFluxReconstruction::computePeclet()
{
  return m_peclet;
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeSmoothness()
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    CFreal stateP = (*((*m_cellStates)[iSol]))[m_monitoredVar];
    CFreal diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];
    sNum += diffStatesPPMinOne*diffStatesPPMinOne;
    sDenom += stateP*stateP;
  }
  if (sNum <= MathTools::MathConsts::CFrealEps() || sDenom <= MathTools::MathConsts::CFrealEps())
  {
    m_s = -100.0;
  }
  else
  {
    m_s = log10(sNum/sDenom);
  }
  
  // get datahandle
  DataHandle< CFreal > smoothness = socket_smoothness.getDataHandle();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    smoothness[(((*m_cellStates)[iSol]))->getLocalID()] = m_s;
  }
  
  if (m_s > m_Smax)
  {
      m_Smax = m_s;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeSmoothness(const CFuint side)
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    CFreal stateP = (*((*m_states[side])[iSol]))[m_monitoredVar];
    CFreal diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];
    sNum += diffStatesPPMinOne*diffStatesPPMinOne;
    sDenom += stateP*stateP;
  }
  if (sNum <= MathTools::MathConsts::CFrealEps() || sDenom <= MathTools::MathConsts::CFrealEps())
  {
    m_s = -100.0;
  }
  else
  {
    m_s = log10(sNum/sDenom);
  }
  
  if (m_s > m_Smax)
  {
      m_Smax = m_s;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::storeEpsilon()
{
  for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
  {
    // get node ID
    const CFuint nodeID = (*m_cellNodes)[iNode]->getLocalID();

    if (!m_useMax) 
    {
      const CFreal newEps = (1.0-m_dampingCoeff)*m_cellEpsilons[m_cell->getID()] + m_dampingCoeff*m_epsilon;
      m_nodeEpsilons[nodeID] += newEps;
      m_cellEpsilons[m_cell->getID()] = newEps;
      m_totalEps += newEps;
    }
    else
    {
      const CFreal maxEps = max(m_epsilon, m_cellEpsilons[m_cell->getID()]);
      m_nodeEpsilons[nodeID] += maxEps;
      m_cellEpsilons[m_cell->getID()] = maxEps;
      m_totalEps += maxEps;
    }
    
    if (m_flagComputeNbNghb)
    {
      m_nbNodeNeighbors[nodeID] += 1.0;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeFlux(const RealVector& values, const std::vector< RealVector* >& gradients, const RealVector& normal, const CFreal& radius, RealVector& flux)
{
  flux = m_diffusiveVarSet->getFlux(values,gradients,normal,radius);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeEpsilon0()
{ 
  const CFreal wavespeed = 1;
  
  //const CFreal deltaKsi = 2.0/(m_order+2.0);
  
  const CFreal peclet = computePeclet();
  
  m_epsilon0 = max(wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::computeEpsilon0(const CFuint side)
{ 
  const CFreal wavespeed = 1;
  
  //const CFreal deltaKsi = 2.0/(m_order+2.0);
  
  const CFreal peclet = computePeclet();
  
  m_epsilon0 = max(wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::setup()
{
  CFAUTOTRACE;

  // setup parent class
  DiffRHSFluxReconstruction::setup();
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  
  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();
  
  // get cell-node connectivity
  m_cellNodesConn = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  // resize the physical data temporary vector
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  convTerm->resizePhysicalData(m_pData);
  
  m_subcellRes = frLocalData[0]->getSubcellResolution();
  
  m_order = static_cast<CFuint>(order);
  
  // get the coefs for extrapolation of the node artificial viscosities to the flx pnts
  m_nodePolyValsAtFlxPnts = frLocalData[0]->getNodePolyValsAtPnt(*(frLocalData[0]->getFlxPntsLocalCoords()));
  
  // get the coefs for extrapolation of the node artificial viscosities to the sol pnts
  m_nodePolyValsAtSolPnts = frLocalData[0]->getNodePolyValsAtPnt(*(frLocalData[0]->getSolPntsLocalCoords()));
  
  // number of cell corner nodes
  /// @note in the future, hanging nodes should be taken into account here
  m_nbrCornerNodes = frLocalData[0]->getNbrCornerNodes();
  
  // get the number of nodes in the mesh
  const CFuint nbrNodes = MeshDataStack::getActive()->getNbNodes();
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  // get the number of cells in the mesh
  const CFuint nbrCells = (*elemType)[0].getEndIdx();
  
  // get datahandle
  DataHandle< CFreal > artVisc = socket_artVisc.getDataHandle();
  DataHandle< CFreal > monPhysVar = socket_monPhysVar.getDataHandle();
  DataHandle< CFreal > smoothness = socket_smoothness.getDataHandle();
  
  const CFuint nbStates = nbrCells*m_nbrSolPnts;

  // resize socket
  artVisc.resize(nbStates);
  monPhysVar.resize(nbStates);
  smoothness.resize(nbStates);
  
  m_nodeEpsilons.resize(nbrNodes);
  m_nbNodeNeighbors.resize(nbrNodes);
  m_cellEpsilons.resize(nbrCells);
  m_solEpsilons.resize(2);
  m_solEpsilons[LEFT].resize(m_nbrSolPnts);
  m_solEpsilons[RIGHT].resize(m_nbrSolPnts);
  m_epsilonLR.resize(2);
  m_epsilonLR[LEFT].resize(m_nbrFaceFlxPnts);
  m_epsilonLR[RIGHT].resize(m_nbrFaceFlxPnts);
  m_unitNormalFlxPnts2.resize(m_nbrFaceFlxPnts);
  m_tempSolPntVec.resize(m_nbrSolPnts);
  m_tempSolPntVec2.resize(m_nbrSolPnts);
  m_cellGradsAV.resize(2);
  m_cellGradFlxPntAV.resize(2);
  m_cellGradsAV[LEFT].resize(m_nbrSolPnts);
  m_cellGradsAV[RIGHT].resize(m_nbrSolPnts);
  m_cellGradFlxPntAV[LEFT].resize(m_nbrFaceFlxPnts);
  m_cellGradFlxPntAV[RIGHT].resize(m_nbrFaceFlxPnts);
  m_avgGradAV.resize(m_nbrEqs);
  m_cellFlags.resize(nbrCells);
  m_unpertCellDiffRes.resize(2);
  m_unpertCellDiffRes[LEFT ].resize(m_nbrSolPnts*m_nbrEqs);
  m_unpertCellDiffRes[RIGHT].resize(m_nbrSolPnts*m_nbrEqs);
  m_projectedCorrL.resize(m_dim);
  m_projectedCorrR.resize(m_dim);
  m_contFlxWoLLAV.resize(m_nbrSolPnts);
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    RealVector temp(m_nbrEqs);
    temp = 0.0;
    m_statesPMinOne.push_back(temp);
    
    m_contFlxWoLLAV[iSol].resize(m_dim+m_ndimplus);
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_contFlxWoLLAV[iSol][iDim].resize(m_nbrEqs);
    }
  }
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_cellGradFlxPntAV[LEFT][iFlx].push_back(new RealVector(m_dim));
      m_cellGradFlxPntAV[RIGHT][iFlx].push_back(new RealVector(m_dim));
    }
  }
  
  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
  {
    m_avgGradAV[iVar] = new RealVector(m_dim);
  }
  
  SafePtr<RealMatrix> vdm = frLocalData[0]->getVandermondeMatrix();
  
  SafePtr<RealMatrix> vdmInv = frLocalData[0]->getVandermondeMatrixInv();
  
  RealMatrix temp(m_nbrSolPnts,m_nbrSolPnts);
  temp = 0.0;
  if (m_dim == 2)
  {
    if (m_ndimplus==3){  //if Triag
      for (CFuint idx = 0; idx < (m_order)*(m_order+1)/2; ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }
    else{
      for (CFuint idx = 0; idx < (m_order)*(m_order); ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }  
  }
  else if (m_dim == 3)
  {
    if (m_ndimplus==4){  //if Tetra
      for (CFuint idx = 0; idx < (m_order)*(m_order+1)*(m_order+2)/6; ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }
    else{
       for (CFuint idx = 0; idx < (m_order)*(m_order)*(m_order); ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }  
  }
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_unitNormalFlxPnts2[iFlx].resize(m_dim);
  }
  
  m_transformationMatrix.resize(m_nbrSolPnts,m_nbrSolPnts);
  
  m_transformationMatrix = (*vdm)*temp*(*vdmInv);
  
  //m_s0 = -m_s0*log10(static_cast<CFreal>(m_order));
  
  m_Smax = m_s0 + m_kappa;
  
  m_SmaxGlobal = m_Smax;
  
  m_nbNodeNeighbors = 0.0;
  
  m_flagComputeNbNghb = true;
  
  cf_assert(m_monitoredVar < m_nbrEqs);
  
  // get CellToFaceGeBuilders
  m_cellBuilders.resize(2);
  m_cellBuilders[LEFT ] = getMethodData().getCellBuilder();
  m_cellBuilders[RIGHT] = getMethodData().getSecondCellBuilder();

  // get some additional data for cell building
  m_isFaceOnBoundary.resize(2);
  m_nghbrCellSide   .resize(2);
  m_currCellSide    .resize(2);
  m_faceOrients     .resize(2);
  m_faceBCIdx       .resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_isFaceOnBoundary[iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getIsFaceOnBoundary();
    m_nghbrCellSide   [iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getNeighbrCellSide ();
    m_currCellSide    [iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getCurrentCellSide ();
    m_faceOrients     [iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getFaceOrient      ();
    m_faceBCIdx       [iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getFaceBCIdx       ();
  }

  
  // get CellToFaceGeBuilder
  m_cellBuilder      = getMethodData().getCellBuilder();
  m_isFaceOnBoundaryCell = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
  m_nghbrCellSideCell    = m_cellBuilder->getGeoBuilder()->getNeighbrCellSide ();
  m_currCellSideCell     = m_cellBuilder->getGeoBuilder()->getCurrentCellSide ();
  m_faceOrientsCell      = m_cellBuilder->getGeoBuilder()->getFaceOrient      ();
  m_faceBCIdxCell        = m_cellBuilder->getGeoBuilder()->getFaceBCIdx       ();
  
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
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
  {
    deletePtr(m_avgGradAV[iVar]); 
  }
  m_avgGradAV.clear();
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      deletePtr(m_cellGradFlxPntAV[LEFT][iFlx][iGrad]);  
      deletePtr(m_cellGradFlxPntAV[RIGHT][iFlx][iGrad]);
    }
    m_cellGradFlxPntAV[LEFT][iFlx].clear();
    m_cellGradFlxPntAV[RIGHT][iFlx].clear();
  }
  m_cellGradFlxPntAV[LEFT].clear();
  m_cellGradFlxPntAV[RIGHT].clear();
  m_cellGradFlxPntAV.clear();
  
  // unsetup parent class
  DiffRHSFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
