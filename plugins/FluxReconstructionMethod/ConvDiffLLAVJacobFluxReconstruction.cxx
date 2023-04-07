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

#include "FluxReconstructionMethod/ConvDiffLLAVJacobFluxReconstruction.hh"
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
  
ConvDiffLLAVJacobFluxReconstruction::ConvDiffLLAVJacobFluxReconstruction(const std::string& name) :
  DiffRHSJacobFluxReconstruction(name),
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
  m_fluxJacobian(),
  m_riemannFluxJacobian(),
  m_flxPntRiemannFluxDiff(),
  m_flxPntRiemannFluxPert(),
  m_tempFlux(),
  m_gradientFluxJacobian(),
  m_gradVarsToStateJacobian(),
  m_gradientStateJacobian(),
  m_riemannFluxGradJacobian(),
  m_contFlxBackupDiff(),
  m_flxPntRiemannFluxDiffConv(),
  m_varToGradVarDep(),
  m_nbrVarToGradVarDep(),
  m_faceID(),
  m_needToAddSolPnt(),
  m_temp(),
  m_tempOther(),
  m_temp2(),
  m_tempOther2(),
  m_epsJacobian(),
  m_llavFluxJacobian(),
  m_llavRiemannFluxJacobian(),
  m_dampCoeffDiff(),
  m_contFlxWoLLAV(),
  m_updateToSolutionVecTrans(CFNULL),
  m_tempSolVarState(),
  m_tempSolVarState2()
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
    
    m_LLAVSubCellRedistribution = false;
    setParameter( "LLAVSubCellRedistribution", &m_LLAVSubCellRedistribution);
    
    m_LLAVRelax = 1.0;
    setParameter( "LLAVRelaxationFactor", &m_LLAVRelax);
  }

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  DiffRHSJacobFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::defineConfigOptions(Config::OptionList& options)
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
  
  options.addConfigOption< bool,Config::DynamicOption<> >("LLAVSubCellRedistribution","Option to activate subcell AV redistribution.");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("LLAVRelaxationFactor","Relaxation factor for CLLAV correction function.");
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  ConvDiffLLAVJacobFluxReconstruction::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = DiffRHSJacobFluxReconstruction::providesSockets();
  result.push_back(&socket_artVisc);
  result.push_back(&socket_monPhysVar);
  result.push_back(&socket_smoothness);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
ConvDiffLLAVJacobFluxReconstruction::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = DiffRHSJacobFluxReconstruction::needsSockets();
  result.push_back(&socket_solPntNormals);
  result.push_back(&socket_flxPntNormals);
  result.push_back(&socket_cellVolumes);
  result.push_back(&socket_volumes);
  if (m_useWallCutOff) result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ConvDiffLLAVJacobFluxReconstruction::execute()\n");
  
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
      m_cellNodes[0]  = m_cell->getNodes();

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

        // get the sol pnt jacobians
        DataHandle< CFreal > volumes = socket_volumes.getDataHandle();
        
        // get the sol pnt normals
        DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();

        // compute solution points Jacobian determinants and epsilons
	for (CFuint iSide = 0; iSide < 2; ++iSide)
        {   
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            m_solJacobDet[iSide][iSol] = volumes[(*(m_states[iSide]))[iSol]->getLocalID()];

            // reset the states in the flx pnts
            m_solEpsilons[iSide][iSol] = 0.0;

            // loop over the sol pnts to compute the states and grads in the flx pnts
            for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
            {
              // get node local index
              const CFuint nodeIdx = (*(m_cellNodes[iSide]))[iNode]->getLocalID();

              m_solEpsilons[iSide][iSol] += m_nodePolyValsAtSolPnts[iSol][iNode]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
            }
            
            const CFuint solID = (*(m_states[iSide]))[iSol]->getLocalID();
      
            for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
            {
              for (CFuint jDim = 0; jDim < m_dim; ++jDim)
              {
                m_neighbCellFluxProjVects[iSide][iDim][iSol][jDim] = solPntNormals[solID*(m_dim+m_ndimplus)*m_dim+iDim*m_dim+jDim];
              }
            }
          }
	}

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
            m_alphaValues[LEFT][iSol] *= m_LLAVRelax/leftAvAlpha;
            m_alphaValues[RIGHT][iSol] *= m_LLAVRelax/rightAvAlpha;
            
            m_alphaValues[LEFT][iSol] += 1.0-m_LLAVRelax;
            m_alphaValues[RIGHT][iSol] += 1.0-m_LLAVRelax;
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

void ConvDiffLLAVJacobFluxReconstruction::setFaceDataForGradients(CFuint faceID)
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

void ConvDiffLLAVJacobFluxReconstruction::setFaceData(CFuint faceID)
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
      m_cellNodes[iSide] = m_cells[iSide]->getNodes();
      CFuint flxIdx;
      iSide == LEFT ? flxIdx = flxPntIdxL : flxIdx = flxPntIdxR;

      for (CFuint iNode = 0; iNode < m_faceNodes->size(); ++iNode)
      {
	for (CFuint iNodeCell = 0; iNodeCell < m_nbrCornerNodes; ++iNodeCell)
        {
	  if ((*m_faceNodes)[iNode]->getLocalID() == (*(m_cellNodes[iSide]))[iNodeCell]->getLocalID())
	  {
	    // get node local index
            const CFuint nodeIdx = (*m_cellNodesConn)(m_cells[iSide]->getID(),iNodeCell);
	    
            m_epsilonLR[iSide][iFlxPnt] += m_nodePolyValsAtFlxPnts[flxIdx][iNodeCell]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
	  }
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::computeFlxPntStatesAndGrads()
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

void ConvDiffLLAVJacobFluxReconstruction::computeInterfaceFlxCorrection()
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
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFluxDiff[iFlxPnt]);
    
    m_flxPntRiemannFluxDiffConv[iFlxPnt] = m_flxPntRiemannFluxDiff[iFlxPnt];
    
    // compute the convective riemann flux
    m_flxPntRiemannFluxDiffConv[iFlxPnt] -= m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),
									    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]),
									    m_unitNormalFlxPnts[iFlxPnt]);
    
    m_flxPntRiemannFlux[iFlxPnt] = m_flxPntRiemannFluxDiffConv[iFlxPnt];    
    
    // compute artificial part
    // get epsilon
    const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
    
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_flxPntRiemannFlux[iFlxPnt][iVar] += epsilon*((*(m_avgGradAV[iVar]))[iDim])*m_unitNormalFlxPnts[iFlxPnt][iDim];
      }
    }
     
    // compute FI in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
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

void ConvDiffLLAVJacobFluxReconstruction::initJacobianComputation()
{
  CFLog(VERBOSE, "initJacobianComputation\n");
  
  // set block row and column indices, proj vectors and make a backup of discontinuous fluxes
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {    
    // Loop over solution points to calculate the discontinuous flux.
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    { 
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }

      m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());
      
      m_updateVarSet->computePhysicalData(*((*(m_states[m_pertSide]))[m_pertSol]), m_pData); 

      prepareFluxComputation();

      // calculate the discontinuous flux projected on x, y, z-directions
      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      {
        // diffusive part 
        computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxBackupDiff[m_pertSide][m_pertSol][iDim]);
        
        m_contFlxBackup[m_pertSide][m_pertSol][iDim] = m_contFlxBackupDiff[m_pertSide][m_pertSol][iDim];
        
        // convective part
        m_contFlxBackup[m_pertSide][m_pertSol][iDim] -= m_updateVarSet->getFlux()(m_pData,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::computeCellFluxJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "computeCellFluxJacobianNum\n");
    
  // loop over left and right cell to compute flux jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // cell ID of the cell at the non-perturbed side
    const CFuint otherCellID = m_cells[iOtherSide]->getID();

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
        
        // compute perturbed physical data
        m_updateVarSet->computePhysicalData(*((*(m_states[m_pertSide]))[m_pertSol]), m_pData); 

        m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());
        
        // compute perturbed fluxes
        prepareFluxComputation();

        // calculate the discontinuous flux projected on x, y, z-directions
        for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
        {       
          // diffusive part
          computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxNeighb[m_pertSide][m_pertSol][iDim]);
                    
          // convective part
          m_contFlxNeighb[m_pertSide][m_pertSol][iDim] -= m_updateVarSet->getFlux()(m_pData,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol]);

          // compute the flux current jacobian term
          // compute the finite difference derivative of the face term (note an implicit minus sign is added here, by the ordering of arguments)
          m_numJacob->computeDerivative(m_contFlxNeighb[m_pertSide][m_pertSol][iDim],m_contFlxBackup[m_pertSide][m_pertSol][iDim],m_tempFlux);

          // multiply residual update derivatives with residual factor so it is taken into the final jacobian
          m_tempFlux *= resFactor;
        
          // store the flux jacobian
          m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][iDim] = m_tempFlux;
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::computeRiemannFluxJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "computeRiemannFluxJacobianNum\n");
    
  // loop over the face flux points to compute the Riemann flux jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // dereference state
      State& pertState = *(m_cellStatesFlxPnt[m_pertSide][iFlxPnt]);
      
      // compute the average grad to use the BR2 scheme
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;            
      }
                
      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
        
        // compute the average sol
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {        
          m_avgSol[iVar] = (pertState[iVar] + (*(m_cellStatesFlxPnt[iOtherSide][iFlxPnt]))[iVar])/2.0; 
        }
    
        prepareFluxComputation();
     
        // compute diffusive flux
        computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFluxPert[iFlxPnt]);
        
        if (m_pertSide == LEFT)
        {
          // compute the convective riemann flux
          m_flxPntRiemannFluxPert[iFlxPnt] -= m_riemannFluxComputer->computeFlux(pertState,
									              *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]),
									              m_unitNormalFlxPnts[iFlxPnt]);
        }
        else
        {
          // compute the convective riemann flux
          m_flxPntRiemannFluxPert[iFlxPnt] -= m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),
                                                                                      pertState,
									              m_unitNormalFlxPnts[iFlxPnt]);
        }
       
        // compute the flux current jacobian term
        // compute the finite difference derivative of the face term (note an implicit minus sign is added here, by the ordering of arguments)
        m_numJacob->computeDerivative(m_flxPntRiemannFluxPert[iFlxPnt],m_flxPntRiemannFluxDiffConv[iFlxPnt],m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]);

        // multiply residual update derivatives with residual factor so it is taken into the final jacobian
        m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar] *= resFactor;

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::computeFluxToGradJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "computeFluxToGradJacobianNum\n");
    
  // loop over left and right cell to compute flux jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        for (CFuint pertDim = 0; pertDim < m_dim; ++pertDim)
        {
          // perturb physical variable in state
          m_numJacob->perturb(m_pertVar,(*(m_tempGrad[m_pertVar]))[pertDim]);
        
          // compute perturbed fluxes
          prepareFluxComputation();

          // calculate the discontinuous flux projected on x, y, z-directions
          for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
          {       
            // diffusive part
            computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxNeighb[m_pertSide][m_pertSol][iDim]);

            // compute the flux current jacobian term
            // compute the finite difference derivative of the face term (note an implicit minus sign is added here, by the ordering of arguments)
            m_numJacob->computeDerivative(m_contFlxNeighb[m_pertSide][m_pertSol][iDim],m_contFlxBackupDiff[m_pertSide][m_pertSol][iDim],m_tempFlux);

            // multiply residual update derivatives with residual factor so it is taken into the final jacobian
            m_tempFlux *= resFactor;
        
            // store the flux jacobian
            m_gradientFluxJacobian[m_pertSide][m_pertSol][m_pertVar][pertDim][iDim] = m_tempFlux;
          }

          // restore physical variable in state
          m_numJacob->restore((*(m_tempGrad[m_pertVar]))[pertDim]);
        }
      }
    }
  }  
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::computeRiemannFluxToGradJacobianNum(const CFreal resFactor)
{  
  CFLog(VERBOSE, "computeRiemannFluxToGradJacobianNum\n");
    
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {      
    // compute the average sol
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {        
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
    }
    
    // compute the average grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;            
    }
           
    // loop over the variables in the state
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {
      for (CFuint pertDir = 0; pertDir < m_dim; ++pertDir)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,(*(m_avgGrad[m_pertVar]))[pertDir]);

        prepareFluxComputation();
     
        // compute diffusive flux
        computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFluxPert[iFlxPnt]);
     
        // compute the flux current jacobian term
        // compute the finite difference derivative of the face term (note an implicit minus sign is added here, by the ordering of arguments)
        m_numJacob->computeDerivative(m_flxPntRiemannFluxPert[iFlxPnt],m_flxPntRiemannFluxDiff[iFlxPnt],m_riemannFluxGradJacobian[iFlxPnt][m_pertVar][pertDir]);

        // multiply residual update derivatives with residual factor so it is taken into the final jacobian
        m_riemannFluxGradJacobian[iFlxPnt][m_pertVar][pertDir] *= resFactor;

        // restore physical variable in state
        m_numJacob->restore((*(m_avgGrad[m_pertVar]))[pertDir]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::computeGradToStateJacobianAna()
{
  CFLog(VERBOSE, "computeGradToStateJacobianAna\n");
  
  // reset the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
      {
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        { 
          m_gradientStateJacobian[m_pertSide][m_pertSol][LEFT][jSol][iDim] = 0.0;
          m_gradientStateJacobian[m_pertSide][m_pertSol][RIGHT][jSol][iDim] = 0.0;
        }
      }
    }
  }
  
  // get the face flux point normals
  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
    
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // loop over the sol pnts of which the gradients will be derived
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // loop over the sol pnts to which to derive
      for (CFuint jSol = 0; jSol < m_nbrSolSolDep; ++jSol)
      {
        const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSol];  
        
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {                    
          for (CFuint jDim = 0; jDim < m_dim; ++jDim)
          {
            m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][m_pertSol][iDim] += m_neighbCellFluxProjVects[m_pertSide][jDim+m_ndimplus][m_pertSol][iDim] * (*m_solPolyDerivAtSolPnts)[jSolIdx][jDim][m_pertSol];
          }
        }
      }
    }
    
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      // local flux point indices in the left and right cell
      const CFuint flxPntIdxThis = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlx];
      const CFuint flxPntIdxOther = (*m_faceFlxPntConnPerOrient)[m_orient][iOtherSide][iFlx];
      
      for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
      {
        const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol]; 
    
        for (CFuint jSol = 0; jSol < m_nbrSolDep; ++jSol)
        {
          const CFuint jSolIdx = (*m_flxSolDep)[flxPntIdxThis][jSol];  
          const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSol];  
        
          for (CFuint iDim = 0; iDim < m_dim; ++iDim)
          {
            m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][pertSolIdx][iDim] -= 0.5 * m_corrFctDiv[jSolIdx][flxPntIdxThis] * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx] *
                    (*m_faceMappedCoordDir)[m_orient][m_pertSide]*m_faceJacobVecs[iFlx][iDim];
         
            m_gradientStateJacobian[iOtherSide][jSolIdxOther][m_pertSide][pertSolIdx][iDim] += 0.5 * m_corrFctDiv[jSolIdxOther][flxPntIdxOther] * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx] *
                    (*m_faceMappedCoordDir)[m_orient][iOtherSide]*m_faceJacobVecs[iFlx][iDim];
          }
        }   
      }
    }
  
    // Add the contribution of the correction to the gradients for each face
    // compute other face contributions to the gradients
    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[m_pertSide].size();
  
    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
    {
      // get local face index
      const CFuint faceIdx = m_otherFaceLocalIdxs[m_pertSide][iFace];
      
      const CFuint faceID = (*m_faces[m_pertSide])[faceIdx]->getID();
      
      if ((*m_isFaceOnBoundary[m_pertSide])[faceIdx])
      {  
        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {
          const CFuint currFlxIdx = (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
          
          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {
            const CFuint pertSolIdx = (*m_flxSolDep)[currFlxIdx][m_pertSol];        
        
            for (CFuint jSol = 0; jSol < m_nbrSolDep; ++jSol)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[currFlxIdx][jSol]; 
    
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][pertSolIdx][iDim] -= 0.5 * m_corrFctDiv[jSolIdx][currFlxIdx] * (*m_solPolyValsAtFlxPnts)[currFlxIdx][pertSolIdx] *
                          (*m_faceLocalDir)[faceIdx]*flxPntNormals[faceID*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
              }
            }
          }
        }
      }
      else
      {
        // Get orientation of face
        const CFuint orient = (*m_faceOrients[m_pertSide])[faceIdx];
        
        // cell side with respect to this face
        const CFuint cellSide = (*m_currCellSide[m_pertSide])[faceIdx];
      
        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {
          const CFuint currFlxIdx = (*m_faceFlxPntConnPerOrient)[orient][cellSide][iFlxPnt];
          
          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {
            const CFuint pertSolIdx = (*m_flxSolDep)[currFlxIdx][m_pertSol];        
        
            for (CFuint jSol = 0; jSol < m_nbrSolDep; ++jSol)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[currFlxIdx][jSol]; 
    
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][pertSolIdx][iDim] -= 0.5 * m_corrFctDiv[jSolIdx][currFlxIdx] * (*m_solPolyValsAtFlxPnts)[currFlxIdx][pertSolIdx] *
                          ((*m_faceMappedCoordDir)[orient][cellSide])*flxPntNormals[faceID*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
              }
            }
          }
        } 
      }
    }
  }
  
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      
      
      for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
      {
        const CFreal invJacobL = 1./m_solJacobDet[LEFT][jSol];
        const CFreal invJacobR = 1./m_solJacobDet[RIGHT][jSol];
          
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        { 
          m_gradientStateJacobian[LEFT][jSol][m_pertSide][m_pertSol][iDim] *= invJacobL;
          m_gradientStateJacobian[RIGHT][jSol][m_pertSide][m_pertSol][iDim] *= invJacobR;
        }
//        if (m_cells[LEFT]->getID()==1) CFLog(INFO,"side: " << m_pertSide << ", sol: " << m_pertSol << ", to side: 0, sol: " << jSol << ": " << m_gradientStateJacobian[LEFT][jSol][m_pertSide][m_pertSol] <<"\n");
//      if (m_cells[LEFT]->getID()==1) CFLog(INFO,"side: " << m_pertSide << ", sol: " << m_pertSol << ", to side: 1, sol: " << jSol << ": " << m_gradientStateJacobian[RIGHT][jSol][m_pertSide][m_pertSol] <<"\n");
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::computeGradVarsToStateJacobianNum()
{
  CFLog(VERBOSE, "computeGradVarsToStateJacobianNum\n");

  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][iEq] = 0.0;
        }
        m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][m_pertVar] = 1.0;
      }
    }
  }
  
  ////@TODO find a better way to to this
  
  // get current iteration
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();
  
  if (iter == 1 && m_faceID == innerFacesStartIdxs[0])
  {
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {
      m_varToGradVarDep[m_pertVar].resize(0);
      
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        if (m_gradVarsToStateJacobian[0][0][m_pertVar][iEq] != 0.0) m_varToGradVarDep[m_pertVar].push_back(iEq);
      }
      m_nbrVarToGradVarDep[m_pertVar] = m_varToGradVarDep[m_pertVar].size();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::computeLLAVCellFluxJacobianAna(const CFreal resFactor)
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
          for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
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
      CFreal w = 0.0;
      CFreal wOther = 0.0;
        
      for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
      {
        // get node local index
        const CFuint nodeIdx = (*(m_cellNodes[m_pertSide]))[iNode]->getLocalID();
      
        w += m_nodePolyValsAtSolPnts[m_pertSol][iNode]/m_nbNodeNeighbors[nodeIdx];
        
        for (CFuint jNode = 0; jNode < m_nbrCornerNodes; ++jNode)
        {
          const CFuint nodeIdxOther = (*(m_cellNodes[iOtherSide]))[jNode]->getLocalID();
            
          if (nodeIdx == nodeIdxOther)
          {
            wOther += m_nodePolyValsAtSolPnts[m_pertSol][iNode]/m_nbNodeNeighbors[nodeIdx];
          }
        }
      }

      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {          
        for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
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

void ConvDiffLLAVJacobFluxReconstruction::computeLLAVRiemannFluxJacobianAna(const CFreal resFactor)
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
          
      CFreal w = 0.0;
      
      for (CFuint iNode = 0; iNode < m_faceNodes->size(); ++iNode)
      {
        for (CFuint iNodeCell = 0; iNodeCell < m_nbrCornerNodes; ++iNodeCell)
        {
	  if ((*m_faceNodes)[iNode]->getLocalID() == (*(m_cellNodes[m_pertSide]))[iNodeCell]->getLocalID())
	  {
	    // get node local index
            const CFuint nodeIdx = (*m_cellNodesConn)(m_cells[m_pertSide]->getID(),iNodeCell);
	    
            w += m_nodePolyValsAtFlxPnts[flxIdx][iNodeCell]/m_nbNodeNeighbors[nodeIdx];
	  }
	}
      }
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {     
        m_tempSolVarState = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_cellStatesFlxPnt[LEFT][iFlx]));
        m_tempSolVarState2 = static_cast<RealVector&>(*m_updateToSolutionVecTrans->transform(m_cellStatesFlxPnt[RIGHT][iFlx]));
    
        *(m_avgGradAV[iVar]) = (*(m_cellGradFlxPntAV[LEFT][iFlx][iVar]) + *(m_cellGradFlxPntAV[RIGHT][iFlx][iVar]))/2.0;
        
        // compute damping term for LLAV
        const RealVector dGradVarXNormalAV = (m_tempSolVarState[iVar] - m_tempSolVarState2[iVar])*m_unitNormalFlxPnts[iFlx];
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

void ConvDiffLLAVJacobFluxReconstruction::computeBothJacobsDiffFaceTerm()
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
  
  if (m_addRiemannToGradJacob || m_addRiemannToGradCrossCellJacob) computeRiemannFluxToGradJacobianNum(resFactor);
  
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
          
              m_tempFlux += m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][iDim+m_ndimplus] * polyCoef;
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
                m_tempFlux += m_llavFluxJacobian[m_pertSide][m_pertSol][m_pertVar][m_pertSide][jSolIdx][iDim+m_ndimplus] * dl;
                
                // (b)
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal dl_dqdu = dl * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                    
                  // (p)
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal dl_dqdu_dudu = m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][var] * dl_dqdu;
                      
                    m_tempFlux += m_gradientFluxJacobian[m_pertSide][kSolIdx][var][jDim][iDim+m_ndimplus] * dl_dqdu_dudu;
                  }

                  CFreal llavPart = m_solEpsilons[m_pertSide][kSolIdx] * m_neighbCellFluxProjVects[m_pertSide][iDim+m_ndimplus][kSolIdx][jDim];
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
            
            if (m_addRiemannToGradJacob)
            {
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
          }
          
          if (m_addRiemannToGradCrossCellJacob)
          {
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
        }
        
        //// add the cross-element gradient part
        if (m_addFluxToGradCrossCellJacob)
        {
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
                m_tempFlux += m_llavFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iOtherSide][jSolIdx][iDim+m_ndimplus] * lOther;
                  
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                { 
                  const CFreal l_dqduOther = lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal l_dqduOther_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * l_dqduOther;
                      
                    m_tempFlux += m_gradientFluxJacobian[iOtherSide][kSolIdxOther][var][jDim][iDim+m_ndimplus] * l_dqduOther_dudu;
                  }
                  
                  CFreal llavPart = l_dqduOther * m_solEpsilons[iOtherSide][kSolIdxOther];
                  
                  llavPart *= m_neighbCellFluxProjVects[iOtherSide][iDim+m_ndimplus][kSolIdxOther][jDim];
                  
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

void ConvDiffLLAVJacobFluxReconstruction::computeOneJacobDiffFaceTerm(const CFuint side)
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
  
  if (m_addRiemannToGradJacob || m_addRiemannToGradCrossCellJacob) computeRiemannFluxToGradJacobianNum(resFactor);
  
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
          
              m_tempFlux += m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][iDim+m_ndimplus] * polyCoef;
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
                m_tempFlux += m_llavFluxJacobian[m_pertSide][m_pertSol][m_pertVar][m_pertSide][jSolIdx][iDim+m_ndimplus] * dl;
                
                // (b)
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal dl_dqdu = dl * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                    
                  // (p)
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal dl_dqdu_dudu = m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][var] * dl_dqdu;
                      
                    m_tempFlux += m_gradientFluxJacobian[m_pertSide][kSolIdx][var][jDim][iDim+m_ndimplus] * dl_dqdu_dudu;
                  }
                  
                  CFreal llavPart = m_solEpsilons[m_pertSide][kSolIdx] * m_neighbCellFluxProjVects[m_pertSide][iDim+m_ndimplus][kSolIdx][jDim];
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
            
            if (m_addRiemannToGradJacob)
            {
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
          }
          
          if (m_addRiemannToGradCrossCellJacob)
          {
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
        }
        
        //// add the cross-element gradient part
        if (m_addFluxToGradCrossCellJacob)
        {
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
                m_tempFlux += m_llavFluxJacobian[m_pertSide][pertSolIdx][m_pertVar][iOtherSide][jSolIdx][iDim+m_ndimplus] * lOther;
                  
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                { 
                  const CFreal l_dqduOther = lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < nbDepGradVar; ++iEq)
                  {
                    const CFuint var = m_varToGradVarDep[m_pertVar][iEq];
                    
                    const CFreal l_dqduOther_dudu = m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][var] * l_dqduOther;
                      
                    m_tempFlux += m_gradientFluxJacobian[iOtherSide][kSolIdxOther][var][jDim][iDim+m_ndimplus] * l_dqduOther_dudu;
                  }
                  
                  CFreal llavPart = l_dqduOther * m_solEpsilons[iOtherSide][kSolIdxOther];
                  
                  llavPart *= m_neighbCellFluxProjVects[iOtherSide][iDim+m_ndimplus][kSolIdxOther][jDim];
                  
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

void ConvDiffLLAVJacobFluxReconstruction::computeUnpertCellDiffResiduals(const CFuint side)
{ 
  // get datahandle
  DataHandle< CFreal > artVisc = socket_artVisc.getDataHandle();
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {     
    artVisc[(((*(m_states[side]))[iSol]))->getLocalID()] = m_solEpsilons[side][iSol];
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
      computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[side][iDim][iSolPnt],0,m_contFlxWoLLAV[iSolPnt][iDim]);
      
      // add convective part
      m_contFlxWoLLAV[iSolPnt][iDim] -= m_updateVarSet->getFlux()(m_pData,m_neighbCellFluxProjVects[side][iDim][iSolPnt]);
      
      m_contFlx[iSolPnt][iDim] = m_contFlxWoLLAV[iSolPnt][iDim];
      
      // add artificial part
      for (CFuint iDim2 = 0; iDim2 < m_dim; ++iDim2)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_contFlx[iSolPnt][iDim][iVar] += m_solEpsilons[side][iSolPnt]*((*m_avgGradAV[iVar])[iDim2])*m_neighbCellFluxProjVects[side][iDim][iSolPnt][iDim2];
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

void ConvDiffLLAVJacobFluxReconstruction::computeGradientFaceCorrections()
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

void ConvDiffLLAVJacobFluxReconstruction::computeGradients()
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

void ConvDiffLLAVJacobFluxReconstruction::computeProjStates(std::vector< RealVector >& projStates)
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

void ConvDiffLLAVJacobFluxReconstruction::computeProjStates(std::vector< RealVector >& projStates, const CFuint side)
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

void ConvDiffLLAVJacobFluxReconstruction::computeEpsilon()
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

CFreal ConvDiffLLAVJacobFluxReconstruction::computePeclet()
{
  return m_peclet;
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::computeSmoothness()
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

void ConvDiffLLAVJacobFluxReconstruction::computeSmoothness(const CFuint side)
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

void ConvDiffLLAVJacobFluxReconstruction::storeEpsilon()
{
  for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
  {
    // get node ID
    const CFuint nodeID = (*(m_cellNodes[0]))[iNode]->getLocalID();

    if (!m_useMax) 
    {
      const CFreal newEps = (1.0-m_dampingCoeff)*m_cellEpsilons[m_cell->getID()] + m_dampingCoeff*m_epsilon;
      m_nodeEpsilons[nodeID] += newEps;
      m_cellEpsilons[m_cell->getID()] = newEps;
      //if(m_cell->getID()==1) CFLog(INFO, "eps: " << newEps << ", eps0: " << m_epsilon0 << ", S: " << m_s << "\n");
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

void ConvDiffLLAVJacobFluxReconstruction::computeFlux(const RealVector& values, const std::vector< RealVector* >& gradients, const RealVector& normal, const CFreal& radius, RealVector& flux)
{
  flux = m_diffusiveVarSet->getFlux(values,gradients,normal,radius);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::computeEpsilon0()
{ 
  const CFreal wavespeed = 1;
  
  //const CFreal deltaKsi = 2.0/(m_order+2.0);
  
  const CFreal peclet = computePeclet();
  
  m_epsilon0 = max(wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::computeEpsilon0(const CFuint side)
{ 
  const CFreal wavespeed = 1;
  
  //const CFreal deltaKsi = 2.0/(m_order+2.0);
  
  const CFreal peclet = computePeclet();
  
  m_epsilon0 = max(wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;

  // setup parent class
  DiffRHSJacobFluxReconstruction::setup();

  m_unpertAllCellDiffRes.resize(0);
  
  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();
  
  // resize the physical data temporary vector
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  convTerm->resizePhysicalData(m_pData);
  
  // get cell-node connectivity
  m_cellNodesConn = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
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

  // get damping coeff
  m_dampCoeffDiff = getMethodData().getDiffDampCoefficient();
  
  m_tempSolVarState.resize(m_nbrEqs);
  m_tempSolVarState2.resize(m_nbrEqs);
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
  m_fluxJacobian.resize(2);
  m_fluxJacobian[LEFT].resize(m_nbrSolPnts);
  m_fluxJacobian[RIGHT].resize(m_nbrSolPnts);
  m_riemannFluxJacobian.resize(2);
  m_riemannFluxJacobian[LEFT].resize(m_nbrFaceFlxPnts);
  m_riemannFluxJacobian[RIGHT].resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFluxDiff.resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFluxDiffConv.resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFluxPert.resize(m_nbrFaceFlxPnts);
  m_tempFlux.resize(m_nbrEqs);
  m_gradientFluxJacobian.resize(2);
  m_gradientFluxJacobian[LEFT].resize(m_nbrSolPnts);
  m_gradientFluxJacobian[RIGHT].resize(m_nbrSolPnts);
  m_gradientStateJacobian.resize(2);
  m_gradientStateJacobian[LEFT].resize(m_nbrSolPnts);
  m_gradientStateJacobian[RIGHT].resize(m_nbrSolPnts);
  m_riemannFluxGradJacobian.resize(m_nbrFaceFlxPnts);
  m_gradVarsToStateJacobian.resize(2);
  m_gradVarsToStateJacobian[LEFT].resize(m_nbrSolPnts);
  m_gradVarsToStateJacobian[RIGHT].resize(m_nbrSolPnts);
  m_contFlxBackupDiff.resize(2);
  m_contFlxBackupDiff[LEFT].resize(m_nbrSolPnts);
  m_contFlxBackupDiff[RIGHT].resize(m_nbrSolPnts);
  m_nbrVarToGradVarDep.resize(m_nbrEqs);
  m_varToGradVarDep.resize(m_nbrEqs);
  m_needToAddSolPnt.resize(m_nbrSolPnts);
  m_temp.resize(m_nbrEqs);
  m_tempOther.resize(m_nbrEqs);
  m_temp2.resize(m_nbrEqs);
  m_tempOther2.resize(m_nbrEqs);
  m_epsJacobian.resize(2);
  m_epsJacobian[LEFT].resize(m_nbrSolPnts);
  m_epsJacobian[RIGHT].resize(m_nbrSolPnts);
  m_llavFluxJacobian.resize(2);
  m_llavFluxJacobian[LEFT].resize(m_nbrSolPnts);
  m_llavFluxJacobian[RIGHT].resize(m_nbrSolPnts);
  m_cellNodes.resize(2);
  m_llavRiemannFluxJacobian.resize(2);
  m_llavRiemannFluxJacobian[LEFT].resize(m_nbrSolPnts);
  m_llavRiemannFluxJacobian[RIGHT].resize(m_nbrSolPnts);
  m_contFlxWoLLAV.resize(m_nbrSolPnts);

  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    m_contFlxBackupDiff[LEFT][iSolPnt].resize(m_dim+m_ndimplus);
    m_contFlxBackupDiff[RIGHT][iSolPnt].resize(m_dim+m_ndimplus);
    m_epsJacobian[LEFT][iSolPnt].resize(m_nbrEqs);
    m_epsJacobian[RIGHT][iSolPnt].resize(m_nbrEqs);
    m_llavFluxJacobian[LEFT][iSolPnt].resize(m_nbrEqs);
    m_llavFluxJacobian[RIGHT][iSolPnt].resize(m_nbrEqs);
    m_llavRiemannFluxJacobian[LEFT][iSolPnt].resize(m_nbrEqs);
    m_llavRiemannFluxJacobian[RIGHT][iSolPnt].resize(m_nbrEqs);
    
    m_contFlxWoLLAV[iSolPnt].resize(m_dim+m_ndimplus);
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_contFlxWoLLAV[iSolPnt][iDim].resize(m_nbrEqs);
    }
    
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_llavFluxJacobian[LEFT][iSolPnt][iVar].resize(2);
      m_llavFluxJacobian[RIGHT][iSolPnt][iVar].resize(2); 
      m_llavRiemannFluxJacobian[LEFT][iSolPnt][iVar].resize(m_nbrFaceFlxPnts);
      m_llavRiemannFluxJacobian[RIGHT][iSolPnt][iVar].resize(m_nbrFaceFlxPnts);
      
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        m_llavRiemannFluxJacobian[LEFT][iSolPnt][iVar][iFlx].resize(m_nbrEqs);
        m_llavRiemannFluxJacobian[RIGHT][iSolPnt][iVar][iFlx].resize(m_nbrEqs);
      }
      
      for (CFuint iSide = 0; iSide < 2; ++iSide)
      {
        m_llavFluxJacobian[LEFT][iSolPnt][iVar][iSide].resize(m_nbrSolPnts);
        m_llavFluxJacobian[RIGHT][iSolPnt][iVar][iSide].resize(m_nbrSolPnts); 
        
        for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
        {
            
          m_llavFluxJacobian[LEFT][iSolPnt][iVar][iSide][jSol].resize(m_dim+m_ndimplus);
          m_llavFluxJacobian[RIGHT][iSolPnt][iVar][iSide][jSol].resize(m_dim+m_ndimplus);
        
          for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
          {
            m_llavFluxJacobian[LEFT][iSolPnt][iVar][iSide][jSol][iDim].resize(m_nbrEqs);
            m_llavFluxJacobian[RIGHT][iSolPnt][iVar][iSide][jSol][iDim].resize(m_nbrEqs);
          }
        }
      }
    }

    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_contFlxBackupDiff[LEFT][iSolPnt][iDim].resize(m_nbrEqs);
      m_contFlxBackupDiff[RIGHT][iSolPnt][iDim].resize(m_nbrEqs);
    }
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    RealVector temp(m_nbrEqs);
    temp = 0.0;
    m_statesPMinOne.push_back(temp);
    
    m_fluxJacobian[LEFT][iSol].resize(m_nbrEqs);
    m_fluxJacobian[RIGHT][iSol].resize(m_nbrEqs);
    m_gradientFluxJacobian[LEFT][iSol].resize(m_nbrEqs);
    m_gradientFluxJacobian[RIGHT][iSol].resize(m_nbrEqs);
    m_gradientStateJacobian[LEFT][iSol].resize(2);
    m_gradientStateJacobian[RIGHT][iSol].resize(2);
    m_gradVarsToStateJacobian[LEFT][iSol].resize(m_nbrEqs);
    m_gradVarsToStateJacobian[RIGHT][iSol].resize(m_nbrEqs);
    
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_gradientStateJacobian[LEFT][iSol][iSide].resize(m_nbrSolPnts);
      m_gradientStateJacobian[RIGHT][iSol][iSide].resize(m_nbrSolPnts);
      
      for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
      {
        m_gradientStateJacobian[LEFT][iSol][iSide][jSol].resize(m_dim);
        m_gradientStateJacobian[RIGHT][iSol][iSide][jSol].resize(m_dim);
      }
    }

    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_fluxJacobian[LEFT][iSol][iVar].resize(m_dim+m_ndimplus);
      m_fluxJacobian[RIGHT][iSol][iVar].resize(m_dim+m_ndimplus);
      m_gradientFluxJacobian[LEFT][iSol][iVar].resize(m_dim);
      m_gradientFluxJacobian[RIGHT][iSol][iVar].resize(m_dim);
      m_gradVarsToStateJacobian[LEFT][iSol][iVar].resize(m_nbrEqs);
      m_gradVarsToStateJacobian[RIGHT][iSol][iVar].resize(m_nbrEqs);

      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      {
        m_fluxJacobian[LEFT][iSol][iVar][iDim].resize(m_nbrEqs);
        m_fluxJacobian[RIGHT][iSol][iVar][iDim].resize(m_nbrEqs);
      }
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
        m_gradientFluxJacobian[LEFT][iSol][iVar][iDim].resize(m_dim+m_ndimplus);
        m_gradientFluxJacobian[RIGHT][iSol][iVar][iDim].resize(m_dim+m_ndimplus);
        
        for (CFuint jDim = 0; jDim < m_dim+m_ndimplus; ++jDim)
        {
          m_gradientFluxJacobian[LEFT][iSol][iVar][iDim][jDim].resize(m_nbrEqs);
          m_gradientFluxJacobian[RIGHT][iSol][iVar][iDim][jDim].resize(m_nbrEqs); 
        }
      }
    }
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_riemannFluxJacobian[LEFT][iFlx].resize(m_nbrEqs);
    m_riemannFluxJacobian[RIGHT][iFlx].resize(m_nbrEqs);
    m_riemannFluxGradJacobian[iFlx].resize(m_nbrEqs);
    m_flxPntRiemannFluxDiff[iFlx].resize(m_nbrEqs);
    m_flxPntRiemannFluxDiffConv[iFlx].resize(m_nbrEqs);
    m_flxPntRiemannFluxPert[iFlx].resize(m_nbrEqs);
      
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_riemannFluxJacobian[LEFT][iFlx][iVar].resize(m_nbrEqs);
      m_riemannFluxJacobian[RIGHT][iFlx][iVar].resize(m_nbrEqs);
      m_riemannFluxGradJacobian[iFlx][iVar].resize(m_dim);
        
      m_cellGradFlxPntAV[LEFT][iFlx].push_back(new RealVector(m_dim));
      m_cellGradFlxPntAV[RIGHT][iFlx].push_back(new RealVector(m_dim));
      
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
        m_riemannFluxGradJacobian[iFlx][iVar][iDim].resize(m_nbrEqs);
      }
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
  
  m_alphaValues.resize(2);
  m_alphaValues[LEFT].resize(m_nbrSolPnts);
  m_alphaValues[RIGHT].resize(m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstruction::unsetup()
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
  DiffRHSJacobFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
