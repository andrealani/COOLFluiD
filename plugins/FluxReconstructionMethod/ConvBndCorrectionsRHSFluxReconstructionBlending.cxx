#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstructionBlending.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionMethod/QuadFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/HexaFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TriagFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TetraFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/PrismFluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvBndCorrectionsRHSFluxReconstructionBlending, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvBndCorrectionsRHSFluxReconstructionBlendingProvider("ConvBndCorrectionsRHSBlending");

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSFluxReconstructionBlending::ConvBndCorrectionsRHSFluxReconstructionBlending(const std::string& name):
ConvBndCorrectionsRHSFluxReconstruction(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_gradients("gradients"),
  socket_gradientsAV("gradientsAV"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  socket_faceJacobVecSizeFaceFlxPntsP0("faceJacobVecSizeFaceFlxPntsP0"),
  socket_alpha("alpha"),
  m_updateVarSet(CFNULL),
  m_faceBuilder(CFNULL),
  m_face(),
  m_intCell(),
  m_orient(),
  m_dim(),
  m_cellStates(),
  m_waveSpeedUpd(),
  m_waveSpeedUpdP0(),
  m_nbrEqs(),
  m_nbrSolPnts(),
  m_nbrFaceFlxPnts(),
  m_solPntsLocalCoords(CFNULL),
  m_flxPntsLocalCoords(),
  m_allCellFlxPnts(CFNULL),
  m_flxPntCoords(),
  m_faceFlxPntConn(CFNULL),
  m_faceFlxPntConnP0(CFNULL),
  m_faceConnPerOrient(CFNULL),
  m_faceIntegrationCoefs(CFNULL),
  m_faceIntegrationCoefsPerType(CFNULL),
  m_faceJacobVecAbsSizeFlxPnts(),
  m_cellStatesFlxPnt(),
  m_flxPntGhostSol(),
  m_riemannFluxComputer(CFNULL),
  m_flxPntRiemannFlux(CFNULL),
  m_corrFctComputer(CFNULL),
  m_corrFctDiv(),
  m_corrections(),
  m_faceMappedCoordDir(CFNULL),
  m_unitNormalFlxPnts(),
  m_faceJacobVecSizeFlxPnts(),
  m_gradUpdates(),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_flxLocalCoords(CFNULL),
  m_faceFlxPntsLocalCoordsPerType(CFNULL),
  m_flxSolDep(CFNULL),
  m_nbrSolDep(),
  m_faceJacobVecs(),
  m_projectedCorr(),
  m_mappedFaceNormalDir(),
  m_order()
{
}

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSFluxReconstructionBlending::~ConvBndCorrectionsRHSFluxReconstructionBlending()
{
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
ConvBndCorrectionsRHSFluxReconstructionBlending::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_gradients);
  result.push_back(&socket_gradientsAV);
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);
  result.push_back(&socket_faceJacobVecSizeFaceFlxPntsP0);
  result.push_back(&socket_alpha);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::configure ( Config::ConfigArgs& args )
{
  ConvBndCorrectionsRHSFluxReconstruction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::executeOnTrs()
{
  CFAUTOTRACE;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current bnd face TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();
  
  CFLog(VERBOSE,"ConvBndCorrectionsRHSFluxReconstructionBlending::executeOnTRS: " << faceTrs->getName() << "\n");

  // get bndFacesStartIdxs from FluxReconstructionMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];
  CFLog(VERBOSE,"ConvBndCorrectionsRHSFluxReconstructionBlending:: 1 \n");
  // number of face orientations
  cf_assert(bndFacesStartIdxs.size() != 0);
  CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

  // number of TRs
  const CFuint nbTRs = faceTrs->getNbTRs();
  cf_assert(bndFacesStartIdxs.size() == nbTRs);

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;
  CFLog(VERBOSE,"ConvBndCorrectionsRHSFluxReconstructionBlending:: 2 \n");
  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm() || getMethodData().hasArtificialViscosity();
  m_bcStateComputer->preProcess();
  CFLog(VERBOSE,"ConvBndCorrectionsRHSFluxReconstructionBlending:: 3 \n");
  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    nbOrients = bndFacesStartIdxs[iTR].size()-1;

    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      CFLog(VERBOSE,"m_orient: " << m_orient << "\n");
      
      // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
      m_nbrFaceFlxPnts=(*m_faceFlxPntConn)[m_orient].size();
      
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

        // get the neighbouring cell
        m_intCell = m_face->getNeighborGeo(0);

	// get the states in the neighbouring cell
        m_cellStates = m_intCell->getStates();

	CFLog(VERBOSE,"cellID: " << m_intCell->getID() << "\n");
//	if (m_intCell->getID() == 72)
//	{
//	  //CFLog(VERBOSE,"coord state: " << (((*m_cellStates)[0])->getCoordinates()) << "\n");
//	  CFLog(VERBOSE,"face ID: " << m_face->getID() << "\n");
//	}


      /////// NEW HERE --- P0 state ///////
         
      RealMatrix transformationMatrixL;
      RealMatrix filterL;
    
      filterL.resize(m_nbrSolPnts, m_nbrSolPnts);
      transformationMatrixL.resize(m_nbrSolPnts, m_nbrSolPnts);

      // Construct the filter matrix 
      filterL = 0.0;
     filterL(0, 0) = 1.0;

      // Compute the transformation matrix
    
      transformationMatrixL = m_vdm * filterL * m_vdmInv;
    
        // Applying filter
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) 
      {
        for (CFuint i = 0; i < m_nbrSolPnts; ++i) 
        { 
         m_P0State[i][iEq] = 0.;
          for (CFuint j = 0; j < m_nbrSolPnts; ++j) 
          { 
           m_P0State[i][iEq] += transformationMatrixL(i, j) * (*((*(m_cellStates))[j]))[iEq];
          }
        }
      }

      //////////////////////////////////////


        // if cell is parallel updatable or the gradients need to be computed, compute the needed cell data
        if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
        {  
	  // set the bnd face data
	  setBndFaceData(m_face->getID());//faceID 
	  
	  // compute the states and ghost states in the flx pnts
	  computeFlxPntStates();
	
	  // compute the interface flux
          computeInterfaceFlxCorrection();

          // compute the wave speed updates
          computeWaveSpeedUpdates(m_waveSpeedUpd,m_waveSpeedUpdP0);

          // update the wave speeds
          updateWaveSpeed();
	}
	  
	  // if cell is parallel updatable, compute the correction flux
	if ((*m_cellStates)[0]->isParUpdatable())
	{

	  // compute the correction -(FI)divh of the bnd face for each sol pnt
          computeCorrection(m_corrections);

          // update the rhs
          updateRHS();
	  
	  // print out the residual updates for debugging
          if((*m_cellStates)[0]->getLocalID() == 704) //m_intCell->getID() == 35)
          {
	    CFLog(VERBOSE, "ID  = " << (*m_cellStates)[0]->getLocalID() << "\n");
            CFLog(VERBOSE, "UpdateBnd = \n");
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
	      CFLog(VERBOSE, "state " << iState << ": " << *(((*m_cellStates)[iState])->getData()) << "\n");
            }
          }
        } 
        
        // if there is a diffusive term, compute the gradients
        if (hasDiffTerm)
        {
          computeGradientBndFaceCorrections();
        }
        
        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::computeInterfaceFlxCorrection()
{ 
  // compute the riemann flux in the flx pnts
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    m_flxPntRiemannFlux[iFlxPnt] = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[iFlxPnt]),
								      *(m_flxPntGhostSol[iFlxPnt]),
								      m_unitNormalFlxPnts[iFlxPnt]);
    
    //CFLog(INFO, "iFlx: " << iFlxPnt << ", ghost: " << *(m_flxPntGhostSol[iFlxPnt]) << ", inner: " << *(m_cellStatesFlxPnt[iFlxPnt]) << ", normal: " << m_unitNormalFlxPnts[iFlxPnt] << "\n");
    
    // store the local Riemann flux, scaled with geometrical Jacobian
    m_flxPntRiemannFlux[iFlxPnt] *= m_faceJacobVecSizeFlxPnts[iFlxPnt];
    
    if (m_intCell->getID() == 0) CFLog(VERBOSE, "flx: " << m_flxPntRiemannFlux[iFlxPnt] << ", normal: " << m_unitNormalFlxPnts[iFlxPnt] << ", jacob: " << m_faceJacobVecSizeFlxPnts[iFlxPnt] << "\n");
  }

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPntsP0; ++iFlxPnt)
  {
    m_flxPntRiemannFluxP0[iFlxPnt] = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPntP0[iFlxPnt]),
								      *(m_flxPntGhostSolP0[iFlxPnt]),
								      m_unitNormalFlxPntsP0[iFlxPnt]);
    
    //CFLog(INFO, "iFlx: " << iFlxPnt << ", ghost: " << *(m_flxPntGhostSol[iFlxPnt]) << ", inner: " << *(m_cellStatesFlxPnt[iFlxPnt]) << ", normal: " << m_unitNormalFlxPnts[iFlxPnt] << "\n");
    
    // store the local Riemann flux, scaled with geometrical Jacobian
    m_flxPntRiemannFluxP0[iFlxPnt] *= m_faceJacobVecSizeFlxPntsP0[iFlxPnt];
    
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::computeFlxPntStates()
{ 
  // Loop over flux points to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // reset the extrapolated states
    *(m_cellStatesFlxPnt[iFlxPnt]) = 0.0;
    
    // get current flx pnt idx
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    
    // extrapolate the states to current flx pnt
    m_nbrSolDep = ((*m_flxSolDep)[currFlxIdx]).size();
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];

      *(m_cellStatesFlxPnt[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(*((*m_cellStates)[solIdx]));
      if (m_intCell->getID() == 783) CFLog(VERBOSE, "sol: " << *((*m_cellStates)[solIdx])  << "\n");
    }
    
    if (m_intCell->getID() == 783) CFLog(VERBOSE, "inner: " << *(m_cellStatesFlxPnt[iFlxPnt])  << "\n");
  }
  
  // compute ghost states
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);

  // Loop over flux points to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPntsP0; ++iFlxPnt)
  {
    // reset the extrapolated states
    *(m_cellStatesFlxPntP0[iFlxPnt]) = 0.0;
    
    // get current flx pnt idx
    const CFuint currFlxIdx = (*m_faceFlxPntConnP0)[m_orient][iFlxPnt];
    
    // extrapolate the states to current flx pnt
    m_nbrSolDep = ((*m_flxSolDepP0)[currFlxIdx]).size();
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDepP0)[currFlxIdx][iSol];
      *(m_cellStatesFlxPntP0[iFlxPnt]) += (*m_solPolyValsAtFlxPntsP0)[currFlxIdx][solIdx]*(((m_P0State)[solIdx]));
    }
  }

  // compute ghost states
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPntP0,m_flxPntGhostSolP0,m_unitNormalFlxPntsP0,m_flxPntCoordsP0);

  if (m_intCell->getID() == 0)
  {
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
     CFLog(VERBOSE, "flx pnt state: " << *(m_cellStatesFlxPnt[iFlxPnt]) << ", ghost: " << *(m_flxPntGhostSol[iFlxPnt]) << "\n") ;
  }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::setBndFaceData(CFuint faceID)
{  
  // get the correct flxPntsLocalCoords depending on the face type (only applied for Prism for now @todo but also needed if hybrid grids)
  if (m_dim>2)
  {
    // get face geo
    const CFGeoShape::Type geo = m_face->getShape(); 

    if (geo == CFGeoShape::TRIAG) // triag face
    {
      (*m_flxLocalCoords) = (*m_faceFlxPntsLocalCoordsPerType)[0];
      (*m_flxLocalCoordsP0) = (*m_faceFlxPntsLocalCoordsPerTypeP0)[0];
    }
    else  // quad face
    {
      (*m_flxLocalCoords) = (*m_faceFlxPntsLocalCoordsPerType)[1];
      (*m_flxLocalCoordsP0) = (*m_faceFlxPntsLocalCoordsPerTypeP0)[1];
    } 
  }

  // compute flux point coordinates if needed for the BC
  if (m_bcStateComputer->needsSpatialCoordinates())
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      m_flxPntCoords[iFlx] = m_face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
    }

    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPntsP0; ++iFlx)
    {
      m_flxPntCoordsP0[iFlx] = m_face->computeCoordFromMappedCoord((*m_flxLocalCoordsP0)[iFlx]);
    }
  }
          
  // compute face Jacobian vectors
  m_faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
  
  m_faceJacobVecsP0 = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoordsP0);

  // communicate the face to the BC class
  m_bcStateComputer->setFace(m_face);
  
  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  // Loop over flux points to set the unit normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
       // CFLog(INFO,"size of m_faceJacobVecAbsSizeFlxPnts = "<< m_faceJacobVecAbsSizeFlxPnts.size()<<" size of faceJacobVecSizeFaceFlxPnts = "<<faceJacobVecSizeFaceFlxPnts.size()<<" faceGlobalID = "<<faceID<<" iFlx = "<<iFlxPnt<<"\n");

    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[faceID][iFlxPnt];
    
    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPnts[iFlxPnt] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient];
    
    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = m_mappedFaceNormalDir*m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
  }

  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPntsP0 = socket_faceJacobVecSizeFaceFlxPntsP0.getDataHandle();

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPntsP0; ++iFlxPnt)
  {
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPntsP0[iFlxPnt] = faceJacobVecSizeFaceFlxPntsP0[faceID][iFlxPnt];

    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPntsP0[iFlxPnt] = m_faceJacobVecAbsSizeFlxPntsP0[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient];

    // set unit normal vector
    m_unitNormalFlxPntsP0[iFlxPnt] = m_mappedFaceNormalDir*m_faceJacobVecsP0[iFlxPnt]/m_faceJacobVecAbsSizeFlxPntsP0[iFlxPnt];

  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::computeCorrection(vector< RealVector >& corrections)
{ 
  cf_assert(corrections.size() == m_nbrSolPnts);

  DataHandle< CFreal > output = socket_alpha.getDataHandle();
  CFreal alpha = output[((*m_cellStates)[0])->getLocalID()];

  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // reset the corrections
    corrections[iSolPnt] = 0.0;
  }
  
  // loop over flx and sol pnts to compute -divhFI
  /*for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];

    const CFuint flxIdxP0 = (*m_faceFlxPntConnP0)[m_flxPntFaceConn[flxIdx]][0];

    // the current correction factor previously computed
    const RealVector& currentCorrFactor = m_flxPntRiemannFlux[iFlxPnt];

    const RealVector& currentCorrFactorP0 = m_flxPntRiemannFluxP0[0];

    cf_assert(currentCorrFactor.size() == m_nbrEqs);

    // compute the term due to each flx pnt
    m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

      // divergence of the correctionfct
      const CFreal divh = m_corrFctDiv[solIdx][flxIdx];

      const CFreal divhP0 = m_corrFctDivP0[0][flxIdxP0];

      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        corrections[solIdx][iVar] -= currentCorrFactor[iVar] * divh * (1.- alpha) + currentCorrFactorP0[iVar] * divhP0 * alpha; 
      }
    }
  }*/

  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over flux points that affect this solution point
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];

      const CFuint flxIdxP0 = (*m_faceFlxPntConnP0)[m_flxPntFaceConn[flxIdx]][0]; 

      // Check if this flux point affects the current solution point
      m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
      bool flxPntAffectsSolPnt = false;
      for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
      {
        if ((*m_flxSolDep)[flxIdx][jSolPnt] == iSolPnt)
        {
          flxPntAffectsSolPnt = true;
          break;
        }
      }
      
      // Only compute corrections if this flux point affects the current solution point
      if (flxPntAffectsSolPnt)
      {
        // the current correction factor corresponding to the interface flux (stored in cellFlx)
        const RealVector& currentCorrFactor = m_flxPntRiemannFlux[iFlxPnt];
        // divergence of the correction function
        const CFreal divh = m_corrFctDiv[iSolPnt][flxIdx];
        // Fill in the corrections
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          corrections[iSolPnt][iVar] -= currentCorrFactor[iVar] * divh * (1.- alpha);// + currentCorrFactorP0[iVar] * divhP0 * alpha);
        }
      }
    }
    // only one flux point for P0 per face
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][0];
    const CFuint flxIdxP0 = (*m_faceFlxPntConnP0)[m_flxPntFaceConn[flxIdx]][0]; 
    const RealVector& currentCorrFactorP0 = m_flxPntRiemannFluxP0[0];
    const CFreal divhP0 = m_corrFctDivP0[0][flxIdxP0];
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      corrections[iSolPnt][iVar] -= currentCorrFactorP0[iVar] * divhP0 * alpha;// * (*m_cellAvgSolCoefs)[iSolPnt];
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::updateRHS()
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
      rhs[resID+iVar] += resFactor*m_corrections[iState][iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::updateWaveSpeed()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle< CFreal > output = socket_alpha.getDataHandle();
  CFreal alpha = output[((*m_cellStates)[0])->getLocalID()];

  // loop over sol pnts to add the wave speed updates to the socket
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();
    updateCoeff[solID] += (m_waveSpeedUpd*(2.0*m_order+1));//*(1.-alpha) + (m_waveSpeedUpdP0*(1))*alpha;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::computeWaveSpeedUpdates(CFreal& waveSpeedUpd, CFreal& waveSpeedUpdP0)
{
  // reset the wave speed update
  waveSpeedUpd = 0.0;
  waveSpeedUpdP0 = 0.0;
  
  // get the correct m_faceIntegrationCoefs depending on the face type (only applicable for Prism for now) @todo should be updated for hybrid grid
  if (m_dim>2)
  {
    // get face geo
    const CFGeoShape::Type geo = m_face->getShape(); 

    if (geo == CFGeoShape::TRIAG) // triag face
    {
      //(*m_faceIntegrationCoefs).resize(m_nbrFaceFlxPnts);
      (m_faceIntegrationCoefs) = &(*m_faceIntegrationCoefsPerType)[0];
      (m_faceIntegrationCoefsP0) = &(*m_faceIntegrationCoefsPerTypeP0)[0];
    }
    else  // quad face
    {
      //(*m_faceIntegrationCoefs).resize(m_nbrFaceFlxPnts);
      (m_faceIntegrationCoefs) = &(*m_faceIntegrationCoefsPerType)[1];
      (m_faceIntegrationCoefsP0) = &(*m_faceIntegrationCoefsPerTypeP0)[1];
    } 
  }

  // loop over flux pnts to compute the updates
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {   
    const CFreal jacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceIntegrationCoefs)[iFlx];

    // transform update states to physical data to calculate eigenvalues
    m_updateVarSet->computePhysicalData(*(m_cellStatesFlxPnt[iFlx]), m_pData);
    waveSpeedUpd += jacobXIntCoef*m_updateVarSet->getMaxAbsEigenValue(m_pData,m_unitNormalFlxPnts[iFlx]);
    if (waveSpeedUpd != waveSpeedUpd) CFLog(INFO, "nan pData: " << m_pData << "\n");
  }
  
  const CFreal jacobXIntCoef = m_faceJacobVecAbsSizeFlxPntsP0[0]*(*m_faceIntegrationCoefsP0)[0];

    // transform update states to physical data to calculate eigenvalues
    m_updateVarSet->computePhysicalData(*(m_cellStatesFlxPntP0[0]), m_pData);
    waveSpeedUpdP0 += jacobXIntCoef*m_updateVarSet->getMaxAbsEigenValue(m_pData,m_unitNormalFlxPntsP0[0]);
    if (waveSpeedUpdP0 != waveSpeedUpdP0) CFLog(INFO, "nan pData: " << m_pData << "\n");

}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::computeGradientBndFaceCorrections()
{ 
  // Loop over solution pnts to reset the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[iSolPnt][iEq] = 0.0;
    }
  }
      
  // compute the face corrections to the gradients
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlx];

    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFreal avgSol = ((*m_cellStatesFlxPnt[iFlx])[iEq]+(*(m_flxPntGhostSol[iFlx]))[iEq])/2.0;
      m_projectedCorr = (avgSol-(*m_cellStatesFlxPnt[iFlx])[iEq])*m_faceJacobVecSizeFlxPnts[iFlx]*m_unitNormalFlxPnts[iFlx];

      // Loop over solution pnts to calculate the grad updates
      m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
      {
        const CFuint iSolIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

	/// @todo Check if this is also OK for triangles!!
	m_gradUpdates[iSolIdx][iEq] += m_projectedCorr*m_corrFctDiv[iSolIdx][flxIdx];
      }
    }
  }
  
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // get state ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradients[solID][iGrad] += m_gradUpdates[iSol][iGrad];
     
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

//void ConvBndCorrectionsRHSFluxReconstructionBlending::preProcess()
//{

//}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::setup()
{
  CFAUTOTRACE;

  // setup parent class
  FluxReconstructionSolverCom::setup();
  
  // get cell builder
  m_faceBuilder = getMethodData().getFaceBuilder();
  
  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();
  
  // get the Riemann flux
  m_riemannFluxComputer = getMethodData().getRiemannFlux();
  
  // get the correction function computer
  m_corrFctComputer = getMethodData().getCorrectionFunction();
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  
  // compute flux point coordinates
  SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  m_nbrFaceFlxPnts = flxLocalCoords->size();

  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_order = static_cast<CFuint>(order);
  
  const CFGeoShape::Type elemShape = frLocalData[0]->getShape();
  
  if (elemShape == CFGeoShape::TETRA)  // numbering convention in tetra requires face->computeFaceJacobDetVectorAtMappedCoords with -1 factor
    {
      m_mappedFaceNormalDir= -1.;
    }
  else
    {
      m_mappedFaceNormalDir= 1.;
    }

    if (elemShape == CFGeoShape::PRISM)  // (Max number of face flx pnts)
    {
      m_nbrFaceFlxPnts=(order+1)*(order+1);
    }

  // number of sol points
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();
  
  // get solution point local coordinates
  m_solPntsLocalCoords = frLocalData[0]->getSolPntsLocalCoords();
   
  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConn();
	  
  // get the face connectivity per orientation
  m_faceConnPerOrient = frLocalData[0]->getFaceConnPerOrient();
  
  // get the face integration coefficient
  m_faceIntegrationCoefs = frLocalData[0]->getFaceIntegrationCoefs();
  
  // get flux point mapped coordinate directions
  m_faceMappedCoordDir = frLocalData[0]->getFaceMappedCoordDir();
  
  // get all flux points of a cell
  m_allCellFlxPnts = frLocalData[0]->getFlxPntsLocalCoords();
  
  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();
  
  // get the face local coords of the flux points on one face
  m_flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();

  // get the face local coords of the flux points on one face depending on the face type
  m_faceFlxPntsLocalCoordsPerType = frLocalData[0]->getFaceFlxPntsLocalCoordsPerType();

  // get the face integration coefficient depending on the face type
  m_faceIntegrationCoefsPerType = frLocalData[0]->getFaceIntegrationCoefsPerType();

  m_cellAvgSolCoefs = frLocalData[0]->getCellAvgSolCoefs();

  m_flxSolDep = frLocalData[0]->getFlxPntSolDependency();

  m_nbrSolDep = ((*m_flxSolDep)[0]).size();

  // create internal and ghost states
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntGhostSol.push_back(new State());
    m_cellStatesFlxPnt.push_back(new State());
  }
  
  // dimensionality and number of equations
  m_dim   = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  
  RealVector dummyCoord;
  dummyCoord.resize(m_dim);
  dummyCoord = 0.0;

  // set an ID as initialization
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntGhostSol[iFlx]->setLocalID(iFlx);
    m_cellStatesFlxPnt[iFlx]->setLocalID(iFlx);
    
    m_flxPntGhostSol[iFlx]->setSpaceCoordinates(new Node(dummyCoord,false));
    m_cellStatesFlxPnt[iFlx]->setSpaceCoordinates(new Node(dummyCoord,false));
  }

  // resize the physical data temporary vector
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  convTerm->resizePhysicalData(m_pData);
  
  // resize vectors
  m_unitNormalFlxPnts.resize(m_nbrFaceFlxPnts);
  m_faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_flxPntsLocalCoords.resize(m_nbrFaceFlxPnts);
  m_faceJacobVecAbsSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_flxPntCoords.resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFlux.resize(m_nbrFaceFlxPnts);
  m_corrections.resize(m_nbrSolPnts);
  m_corrFctDiv.resize(m_nbrSolPnts);
  m_corrFctDivP0.resize(1);
  m_faceJacobVecs.resize(m_nbrFaceFlxPnts);
  m_projectedCorr.resize(m_dim);
  
  m_Filtered_flxPntRiemannFlux.resize(m_nbrFaceFlxPnts);

  m_vdm = *(frLocalData[0]->getVandermondeMatrix());
  
  m_vdmInv = *(frLocalData[0]->getVandermondeMatrixInv());


  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  { 
    m_flxPntsLocalCoords[iFlx].resize(m_dim);
    m_flxPntCoords[iFlx].resize(m_dim);
    m_unitNormalFlxPnts[iFlx].resize(m_dim);
    m_flxPntRiemannFlux[iFlx].resize(m_nbrEqs);
    m_faceJacobVecs[iFlx].resize(m_dim);
    m_Filtered_flxPntRiemannFlux[iFlx].resize(m_nbrEqs);
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_corrections[iSol].resize(m_nbrEqs);
    m_corrFctDiv[iSol].resize(m_allCellFlxPnts->size());
  }
  
  // resize m_gradUpdates
  m_gradUpdates.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_gradUpdates[iSol].resize(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradUpdates[iSol][iEq].resize(m_dim);
    }
  }
  
  // compute the divergence of the correction function
  m_corrFctComputer->computeDivCorrectionFunction(frLocalData[0],m_corrFctDiv);


  // get the local flx pnts - face connectivity
  m_flxPntFaceConn = *(frLocalData[0]->getFlxPntFaceConn());

  // get the local FR data
  FluxReconstructionElementData* frLocalDataP0;

  // get the order of the polynomial interpolation
  const CFPolyOrder::Type polyOrder = CFPolyOrder::ORDER0;

  // Face type index (default 0 for 2D)
  m_faceType = 0;

  switch(elemShape)
  {
    case CFGeoShape::LINE:
    {
      throw Common::NotImplementedException (FromHere(),"FR has not been implemented for 1D");
    } break;
    case CFGeoShape::TRIAG:
    {
      frLocalDataP0 = new TriagFluxReconstructionElementData(polyOrder);
    } break;
    case CFGeoShape::QUAD:
    {
      frLocalDataP0 = new QuadFluxReconstructionElementData(polyOrder);
    } break;
    case CFGeoShape::TETRA:
    {
      frLocalDataP0 = new TetraFluxReconstructionElementData(polyOrder);      
    } break;
    case CFGeoShape::HEXA:
    {
      frLocalDataP0 = new HexaFluxReconstructionElementData(polyOrder);
    } break;
    case CFGeoShape::PRISM:
    {
      frLocalDataP0 = new PrismFluxReconstructionElementData(polyOrder);      
    } break;      
    default:
    {
      throw Common::NotImplementedException (FromHere(),"FR method not implemented for elements of type "
                                    + StringOps::to_str(elemShape) + ".");
    }
  }
  m_corrFctDivP0[0].resize((frLocalDataP0->getFlxPntsLocalCoords())->size());

  // get the face - flx pnt connectivity for P0
  m_faceFlxPntConnP0 = frLocalDataP0->getFaceFlxPntConn();  

  // compute the divergence of the correction function
  m_corrFctComputer->computeDivCorrectionFunction(frLocalDataP0,m_corrFctDivP0);

  m_solPolyValsAtFlxPntsP0 = frLocalDataP0->getCoefSolPolyInFlxPnts();

  // get the face local coords of the flux points on one face
  m_flxLocalCoordsP0 = frLocalDataP0->getFaceFlxPntsFaceLocalCoords();

  // get the face integration coefficient
  m_faceIntegrationCoefsP0 = frLocalDataP0->getFaceIntegrationCoefs();

  // get the face integration coefficient depending on the face type
  m_faceIntegrationCoefsPerTypeP0 = frLocalDataP0->getFaceIntegrationCoefsPerType();
  
  m_faceFlxPntsLocalCoordsPerTypeP0 = frLocalDataP0->getFaceFlxPntsLocalCoordsPerType();

  // compute flux point coordinates
  m_nbrFaceFlxPntsP0 = m_flxLocalCoordsP0->size();

  if (elemShape == CFGeoShape::PRISM)  // (Max number of face flx pnts)
    {     
      m_nbrFaceFlxPntsP0=(0+1)*(0+1);
    }
    
  m_flxSolDepP0 = frLocalDataP0->getFlxPntSolDependency();

  m_P0State.resize(m_nbrSolPnts);

  m_flxPntRiemannFluxP0.resize(m_nbrFaceFlxPntsP0);

  m_unitNormalFlxPntsP0.resize(m_nbrFaceFlxPntsP0);
  m_faceJacobVecSizeFlxPntsP0.resize(m_nbrFaceFlxPntsP0);
  m_faceJacobVecAbsSizeFlxPntsP0.resize(m_nbrFaceFlxPntsP0);
  m_flxPntCoordsP0.resize(m_nbrFaceFlxPntsP0);
  m_faceJacobVecsP0.resize(m_nbrFaceFlxPntsP0);
  
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {

    m_P0State[iSolPnt].resize(m_nbrEqs);

  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPntsP0; ++iFlx)
  { 
    m_flxPntCoordsP0[iFlx].resize(m_dim);
    m_unitNormalFlxPntsP0[iFlx].resize(m_dim);
    m_flxPntRiemannFluxP0[iFlx].resize(m_nbrEqs);
    m_faceJacobVecsP0[iFlx].resize(m_dim);
  }

  // create internal and ghost states
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPntsP0; ++iFlx)
  {
    m_flxPntGhostSolP0.push_back(new State());
    m_cellStatesFlxPntP0.push_back(new State());
  }
  
  // set an ID as initialization
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPntsP0; ++iFlx)
  {
    m_flxPntGhostSolP0[iFlx]->setLocalID(iFlx);
    m_cellStatesFlxPntP0[iFlx]->setLocalID(iFlx);
    
    m_flxPntGhostSolP0[iFlx]->setSpaceCoordinates(new Node(dummyCoord,false));
    m_cellStatesFlxPntP0[iFlx]->setSpaceCoordinates(new Node(dummyCoord,false));
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstructionBlending::unsetup()
{
  CFAUTOTRACE;
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    deletePtr(m_cellStatesFlxPnt[iFlx]);
    deletePtr(m_flxPntGhostSol[iFlx]);
  }
  m_cellStatesFlxPnt.clear();
  m_flxPntGhostSol.clear();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPntsP0; ++iFlx)
  {
    deletePtr(m_cellStatesFlxPntP0[iFlx]);
    deletePtr(m_flxPntGhostSolP0[iFlx]);
  }
  m_cellStatesFlxPntP0.clear();
  m_flxPntGhostSolP0.clear();

  // unsetup parent class
  FluxReconstructionSolverCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
