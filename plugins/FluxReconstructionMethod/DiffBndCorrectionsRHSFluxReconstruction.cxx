#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/DiffBndCorrectionsRHSFluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/GradientVariables.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {
      
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSFluxReconstruction, 
		       FluxReconstructionSolverData, 
		       FluxReconstructionModule >
DiffBndCorrectionsRHSFluxReconstructionProvider("DiffBndCorrectionsRHS");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSFluxReconstruction::DiffBndCorrectionsRHSFluxReconstruction(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_gradients("gradients"),
  socket_gradientsAV("gradientsAV"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  m_faceBuilder(CFNULL),
  m_bcStateComputer(CFNULL),
  m_face(),
  m_intCell(),
  m_orient(),
  m_dim(),
  m_cellStates(),
  m_waveSpeedUpd(),
  m_cellGrads(),
  m_cellGradFlxPnt(),
  m_flxPntGhostGrads(),
  m_cellVolume(),
  m_cflConvDiffRatio(),
  m_nbrEqs(),
  m_nbrSolPnts(),
  m_nbrFaceFlxPnts(),
  m_solPntsLocalCoords(CFNULL),
  m_allCellFlxPnts(CFNULL),
  m_flxPntsLocalCoords(),
  m_flxPntCoords(),
  m_faceFlxPntConn(CFNULL),
  m_faceConnPerOrient(CFNULL),
  m_faceIntegrationCoefs(CFNULL),
  m_faceIntegrationCoefsPerType(CFNULL),
  m_faceJacobVecAbsSizeFlxPnts(),
  m_cellStatesFlxPnt(),
  m_cellFlx(),
  m_flxPntGhostSol(),
  m_riemannFluxComputer(CFNULL),
  m_flxPntRiemannFlux(CFNULL),
  m_corrFctComputer(CFNULL),
  m_corrFctDiv(),
  m_corrections(),
  m_faceMappedCoordDir(CFNULL),
  m_unitNormalFlxPnts(),
  m_faceJacobVecSizeFlxPnts(),
  m_diffusiveVarSet(CFNULL),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_avgSol(),
  m_avgGrad(),
  m_flxLocalCoords(CFNULL),
  m_faceFlxPntsLocalCoordsPerType(CFNULL),
  m_flxSolDep(CFNULL),
  m_nbrSolDep(),
  m_faceJacobVecs(),
  m_mappedFaceNormalDir(),
  m_order(),
  m_br2Eta(),
  m_ndimplus(),
  m_dimList(),
  m_cellFluxProjVects(),
  m_solJacobDet(),
  m_solPolyDerivAtSolPnts(CFNULL),
  m_solSolDep(CFNULL),
  m_nbrSolSolDep(),
  m_projectedCorr(),
  m_gradVarsSolPnts(),
  m_gradVarStatePtrs(),
  m_compactGradsSolPnts(),
  m_flxPntGradVarsStore(),
  m_flxPntGradVars(),
  m_bndGradVarsStore(),
  m_bndGradVars(),
  m_bndGradFlxPnt()
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSFluxReconstruction::~DiffBndCorrectionsRHSFluxReconstruction()
{
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
DiffBndCorrectionsRHSFluxReconstruction::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_gradients);
  result.push_back(&socket_gradientsAV);
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::executeOnTrs()
{
  CFAUTOTRACE;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current QuadFreeBCFluxReconstruction TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();
  
  CFLog(VERBOSE,"DiffBndCorrectionsRHSFluxReconstruction::executeOnTRS: " << faceTrs->getName() << "\n");

  // get bndFacesStartIdxs from FluxReconstructionMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

  // number of face orientations (should be the same for all TRs)
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
	
	// compute volume
        m_cellVolume = m_intCell->computeVolume();
	
	//cf_assert(m_cellVolume > 0.0);
        if (m_cellVolume<0.0) m_cellVolume = -m_cellVolume;
        
	// set the bnd face data
	setBndFaceData(m_face->getID());//faceID

	// compute the states and ghost states in the flx pnts
	computeFlxPntStates();

	// the compact BR2 face gradient of this boundary face
	setIntCellMetrics();
	computeCompactBR2BndFaceGradient();

	// compute FI
        computeInterfaceFlxCorrection();

        // compute the wave speed updates
        computeWaveSpeedUpdates(m_waveSpeedUpd);
      
        // update the wave speeds
        updateWaveSpeed();
	  
	// if cell is parallel updatable, compute the correction flux
        if ((*m_cellStates)[0]->isParUpdatable())
        {

	  // compute the correction -(FI)divh of the bnd face for each sol pnt
          computeCorrection(m_corrections);

	  // update the rhs
          updateRHS();
        } 
        
        // print out the residual updates for debugging
        if(m_intCell->getID() == 35)//true) //
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
        
        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::computeBndFaceDiffData(GeometricEntity* face,
                                                                     const CFuint orient,
                                                                     BndFaceDiffData& data)
{
  m_face = face;
  m_orient = orient;
  m_intCell = face->getNeighborGeo(0);
  m_cellStates = m_intCell->getStates();
  m_cellVolume = std::abs(m_intCell->computeVolume());

  // flux points of this face
  m_nbrFaceFlxPnts = (*m_faceFlxPntConn)[m_orient].size();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntsLocalCoords[iFlx] = (*m_allCellFlxPnts)[(*m_faceFlxPntConn)[m_orient][iFlx]];
  }

  // diffusive boundary flux of this face, as in executeOnTrs
  setBndFaceData(face->getID());
  computeFlxPntStates();
  setIntCellMetrics();
  computeCompactBR2BndFaceGradient();
  computeInterfaceFlxCorrection();

  resizeBndFaceDiffData(data);

  // face type for the integration coefficients of a 3D face
  const CFuint faceType = face->getShape() == CFGeoShape::TRIAG ? 0 : 1;

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    data.intStates[iFlx] = *m_cellStatesFlxPnt[iFlx];
    data.ghostStates[iFlx] = *m_flxPntGhostSol[iFlx];
    data.bndStates[iFlx] = *m_bndStateFlux[iFlx];
    data.gradVarsFlxPnt[iFlx] = *m_flxPntGradVars[iFlx];
    data.bndGradVars[iFlx] = *m_bndGradVars[iFlx];
    data.diffFluxes[iFlx] = m_flxPntRiemannFlux[iFlx];
    data.unitNormals[iFlx] = m_unitNormalFlxPnts[iFlx];
    data.coords[iFlx] = m_flxPntCoords[iFlx];
    data.faceJacobVecAbsSizes[iFlx] = m_faceJacobVecAbsSizeFlxPnts[iFlx];

    if (m_dim == 3)
    {
      data.faceIntegrationCoefs[iFlx] = (*m_faceIntegrationCoefsPerType)[faceType][iFlx];
    }
    else
    {
      data.faceIntegrationCoefs[iFlx] = (*m_faceIntegrationCoefs)[iFlx];
    }

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      data.bndGrads[iFlx][iEq] = *m_bndGradFlxPnt[iFlx][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::resizeBndFaceDiffData(BndFaceDiffData& data)
{
  data.intStates.resize(m_nbrFaceFlxPnts);
  data.ghostStates.resize(m_nbrFaceFlxPnts);
  data.bndStates.resize(m_nbrFaceFlxPnts);
  data.gradVarsFlxPnt.resize(m_nbrFaceFlxPnts);
  data.bndGradVars.resize(m_nbrFaceFlxPnts);
  data.bndGrads.resize(m_nbrFaceFlxPnts);
  data.diffFluxes.resize(m_nbrFaceFlxPnts);
  data.unitNormals.resize(m_nbrFaceFlxPnts);
  data.coords.resize(m_nbrFaceFlxPnts);
  data.faceJacobVecAbsSizes.resize(m_nbrFaceFlxPnts);
  data.faceIntegrationCoefs.resize(m_nbrFaceFlxPnts);

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    data.intStates[iFlx].resize(m_nbrEqs);
    data.ghostStates[iFlx].resize(m_nbrEqs);
    data.bndStates[iFlx].resize(m_nbrEqs);
    data.gradVarsFlxPnt[iFlx].resize(m_nbrEqs);
    data.bndGradVars[iFlx].resize(m_nbrEqs);
    data.diffFluxes[iFlx].resize(m_nbrEqs);
    data.unitNormals[iFlx].resize(m_dim);
    data.coords[iFlx].resize(m_dim);
    data.bndGrads[iFlx].resize(m_nbrEqs);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      data.bndGrads[iFlx][iEq].resize(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::computeFlxPntStates()
{
  // Loop over flux points to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // reset the states in flx pnts
    *(m_cellStatesFlxPnt[iFlxPnt]) = 0.0;

    // index of current flx pnt
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];

    // Loop over sol points to add the contributions to each sol pnt
    m_nbrSolDep = ((*m_flxSolDep)[currFlxIdx]).size();
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];

      *(m_cellStatesFlxPnt[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(*((*m_cellStates)[solIdx]));
    }
  }

  // compute ghost states
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::setIntCellMetrics()
{
  // mapped coordinate plane normals in the solution points
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_cellFluxProjVects[iDim] = m_intCell->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
  }

  // solution point Jacobian determinants
  m_solJacobDet = m_intCell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::computeCompactBR2BndFaceGradient()
{
  // gradient variables at the solution points of the interior cell
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_gradVarStatePtrs[iSol] = (*m_cellStates)[iSol]->getData();
  }

  m_diffusiveVarSet->setGradientVars(m_gradVarStatePtrs,m_gradVarsSolPnts,m_nbrSolPnts);

  // extrapolated to the flux points of this boundary face
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlx];

    extrapolateGradVarsToFlxPnt(m_gradVarsSolPnts,(*m_flxSolDep)[flxIdx],(*m_solPolyValsAtFlxPnts)[flxIdx],m_nbrEqs,*m_flxPntGradVars[iFlx]);
  }

  // the value the boundary condition lifts the extrapolated gradient variables to
  m_bcStateComputer->computeBndGradVars(m_flxPntGradVars,m_cellStatesFlxPnt,m_flxPntGhostSol,
                                        m_unitNormalFlxPnts,m_flxPntCoords,m_bndGradVars);

  // reset
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_compactGradsSolPnts[iSol][iEq] = 0.0;
    }
  }

  // volume term of the gradient of the gradient variables
  addGradVarsVolumeTerm(m_gradVarsSolPnts,m_cellFluxProjVects,*m_solSolDep,*m_solPolyDerivAtSolPnts,m_nbrSolSolDep,
                       m_nbrSolPnts,m_nbrEqs,m_dim,m_projectedCorr,m_compactGradsSolPnts);

  // eta times the lifting of this face: jump = g_b - a
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlx];

    m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFreal gradVarsJump = (*m_bndGradVars[iFlx])[iEq] - (*m_flxPntGradVars[iFlx])[iEq];
      const CFreal etaGradVarsJump = m_br2Eta*gradVarsJump;

      addGradVarsLifting(etaGradVarsJump,m_faceJacobVecSizeFlxPnts[iFlx],m_unitNormalFlxPnts[iFlx],
                        (*m_flxSolDep)[flxIdx],m_nbrSolDep,m_corrFctDiv,flxIdx,iEq,m_projectedCorr,m_compactGradsSolPnts);
    }
  }

  // division by the Jacobian determinant
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const CFreal invJacobDet = 1.0/m_solJacobDet[iSol];

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_compactGradsSolPnts[iSol][iEq] *= invJacobDet;
    }
  }

  // extrapolate to the flux points of this face
  extrapolateGradsToFlxPnts(m_compactGradsSolPnts,(*m_faceFlxPntConn)[m_orient],m_nbrFaceFlxPnts,
                        *m_flxSolDep,*m_solPolyValsAtFlxPnts,m_nbrEqs,m_cellGradFlxPnt);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::computeInterfaceFlxCorrection()
{
  // the boundary states the diffusive flux and the transport properties are evaluated at
  m_bcStateComputer->computeBndStates(m_cellStatesFlxPnt,m_flxPntGhostSol,
                                      m_unitNormalFlxPnts,m_flxPntCoords,m_bndStateFlux);

  // the boundary gradients: the compact face gradients with the rules of the
  // boundary condition on the normal derivatives
  m_bcStateComputer->computeBndGrads(m_cellGradFlxPnt,m_bndGradFlxPnt,m_bndStateFlux,
                                     m_unitNormalFlxPnts,m_flxPntCoords);

  // compute the diffusive flux in the flx pnts
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = *(m_bndGradFlxPnt[iFlxPnt][iVar]);

      m_avgSol[iVar] = (*m_bndStateFlux[iFlxPnt])[iVar];
    }

    // prepare the flux computation
    prepareFlxPntFluxComputation(iFlxPnt);

    // rules of the boundary condition that need the transport properties
    m_bcStateComputer->constrainBndGrads(m_avgSol,m_avgGrad,
                                         m_unitNormalFlxPnts[iFlxPnt],m_flxPntCoords[iFlxPnt]);

    // store the final boundary gradients
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *m_bndGradFlxPnt[iFlxPnt][iVar] = *m_avgGrad[iVar];
    }

    // compute FI
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFlux[iFlxPnt]);

    // compute FI in the local frame
    m_cellFlx[iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::computeFlux(const RealVector& values, const std::vector< RealVector* >& gradients, const RealVector& normal, const CFreal& radius, RealVector& flux)
{
  flux = m_diffusiveVarSet->getFlux(values,gradients,normal,0);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::setBndFaceData(CFuint faceID)
{ 
  // get the correct flxPntsLocalCoords depending on the face type (only applied for Prism for now @todo but also needed if hybrid grids)
  if (m_dim>2)
  {
    // get face geo
    const CFGeoShape::Type geo = m_face->getShape(); 

    if (geo == CFGeoShape::TRIAG) // triag face
    {
      (*m_flxLocalCoords) = (*m_faceFlxPntsLocalCoordsPerType)[0];
    }
    else  // quad face
    {
      (*m_flxLocalCoords) = (*m_faceFlxPntsLocalCoordsPerType)[1];
    } 

  }

  // compute flux point coordinates
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx] = m_face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);	
  }
          
  // compute face Jacobian vectors
  m_faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
  
  // communicate the face to the BC class
  m_bcStateComputer->setFace(m_face);
  
  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
  // Loop over flux points to compute the unit normals
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[faceID][iFlxPnt];
    
    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPnts[iFlxPnt] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient];
    
    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = m_mappedFaceNormalDir*m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
  }
  
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // set gradients
  const CFuint nbrStates = m_cellStates->size();
  m_cellGrads.resize(nbrStates);
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const CFuint stateID = (*m_cellStates)[iState]->getLocalID();
    m_cellGrads[iState] = &gradients[stateID];
  }
  
//  // compute Jacobian determinants
//  m_jacobDets = m_intCell->computeGeometricShapeFunctionJacobianDeterminant((*m_faceFlxPntCellMappedCoords)[m_orient][LEFT]);
//
//  // compute inverse characteristic lengths
//  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
//  {
//    m_faceInvCharLengths[iFlx] = m_faceJacobVecAbsSizeFlxPnts[iFlx]/m_jacobDets[iFlx];
//  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::computeCorrection(vector< RealVector >& corrections)
{ 
  cf_assert(corrections.size() == m_nbrSolPnts);
  
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // reset the corrections
    corrections[iSolPnt] = 0.0;
  }

  // compute the term due to each flx pnt
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];

    // the current correction factor
    const RealVector& currentCorrFactor = m_cellFlx[iFlxPnt];

    m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

      // divergence of the correctionfct
      const CFreal divh = m_corrFctDiv[solIdx][flxIdx];
    
      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        corrections[solIdx][iVar] += currentCorrFactor[iVar] * divh;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::updateRHS()
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
      //CFLog(INFO, "upres " << resID+iVar << ": " << resFactor*m_corrections[iState][iVar] << "\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::updateWaveSpeed()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // add the updates to the wave speeds
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();
    updateCoeff[solID] += m_waveSpeedUpd*(2.0*m_order+1);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::computeWaveSpeedUpdates(CFreal& waveSpeedUpd)
{
  CFreal visc = 1.0;
  
  // get the correct m_faceIntegrationCoefs depending on the face type (only applicable for Prism for now) @todo should be updated for hybrid grid
  if (m_dim>2)
  {
    // get face geo
    const CFGeoShape::Type geo = m_face->getShape(); 

    if (geo == CFGeoShape::TRIAG) // triag face
    {
      //(*m_faceIntegrationCoefs).resize(m_nbrFaceFlxPnts);
      (m_faceIntegrationCoefs) = &(*m_faceIntegrationCoefsPerType)[0];
    }
    else  // quad face
    {
      //(*m_faceIntegrationCoefs).resize(m_nbrFaceFlxPnts);
      (m_faceIntegrationCoefs) = &(*m_faceIntegrationCoefsPerType)[1];
    } 
  }

  waveSpeedUpd = 0.0;
  //for (CFuint iFlx = 0; iFlx < m_cellFlx.size(); ++iFlx)
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
   
    // transform update states to physical data to calculate eigenvalues
    waveSpeedUpd += visc*jacobXJacobXIntCoef/m_cellVolume;
  }
  //if (waveSpeedUpd > 10.0) CFLog(INFO, "wvspDiffBnd: " << waveSpeedUpd << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::setup()
{
  CFAUTOTRACE;

  // setup parent class
  FluxReconstructionSolverCom::setup();
  
  // get cell builder
  m_faceBuilder = getMethodData().getFaceBuilder();
  
  // get the diffusive varset, the Null one gives the states themselves
  m_diffusiveVarSet = getMethodData().getDiffusiveVar();

  // multiplier of the face lifting in the compact BR2 face gradient
  m_br2Eta = getMethodData().getBR2Eta();
  
  // get the Riemann flux
  m_riemannFluxComputer = getMethodData().getRiemannFlux();
  
  // get the correction function computer
  m_corrFctComputer = getMethodData().getCorrectionFunction();
  
  // get the local spectral FD data
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

  //Setting ndimplus, needed for Triag (and tetra, prism)
  if (elemShape == CFGeoShape::TRIAG || elemShape == CFGeoShape::TETRA || elemShape == CFGeoShape::PRISM)
  {
    m_ndimplus=1;
  }
  else
  {
    m_ndimplus=0;
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

  // get the face integration coefficient depending on the face type
  m_faceIntegrationCoefsPerType = frLocalData[0]->getFaceIntegrationCoefsPerType();
  
  // get flux point mapped coordinate directions
  m_faceMappedCoordDir = frLocalData[0]->getFaceMappedCoordDir();
  
  // get all flux points of a cell
  m_allCellFlxPnts = frLocalData[0]->getFlxPntsLocalCoords();
  
  // get convective/diffusive CFL ratio
  m_cflConvDiffRatio = frLocalData[0]->getCFLConvDiffRatio();
  
  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();

  // get the face local coords of the flux points on one face
  m_flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  
  // get the face local coords of the flux points on one face depending on the face type
  m_faceFlxPntsLocalCoordsPerType = frLocalData[0]->getFaceFlxPntsLocalCoordsPerType();

  m_flxSolDep = frLocalData[0]->getFlxPntSolDependency();

  m_nbrSolDep = ((*m_flxSolDep)[0]).size();

  // get the coefs for derivation of the states in the sol pnts
  m_solPolyDerivAtSolPnts = frLocalData[0]->getCoefSolPolyDerivInSolPnts();

  m_solSolDep = frLocalData[0]->getSolPntSolDependency();

  m_nbrSolSolDep = ((*m_solSolDep)[0]).size();

  // create internal and ghost states
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntGhostSol.push_back(new State());
    m_cellStatesFlxPnt.push_back(new State());
  }
  
  // dimensionality and number of equations
  m_dim = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  
  RealVector dummyCoord;
  dummyCoord.resize(m_dim);
  dummyCoord = 0.0;

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntGhostSol[iFlx]->setLocalID(iFlx);
    m_cellStatesFlxPnt[iFlx]->setLocalID(iFlx);
    
    m_flxPntGhostSol[iFlx]->setSpaceCoordinates(new Node(dummyCoord,false));
    m_cellStatesFlxPnt[iFlx]->setSpaceCoordinates(new Node(dummyCoord,false));
  }
  
  // resize m_faceJacobVecSizeFlxPnts
  m_faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  
  // resize m_unitNormalFlxPnts
  m_unitNormalFlxPnts.resize(m_nbrFaceFlxPnts);
  
  // resize vectors
  m_flxPntsLocalCoords.resize(m_nbrFaceFlxPnts);
  m_faceJacobVecAbsSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_cellFlx.resize(m_nbrFaceFlxPnts);
  m_flxPntCoords.resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFlux.resize(m_nbrFaceFlxPnts);
  m_corrections.resize(m_nbrSolPnts);
  m_cellGradFlxPnt.resize(m_nbrFaceFlxPnts);
  m_flxPntGhostGrads.resize(m_nbrFaceFlxPnts);
  m_avgSol.resize(m_nbrEqs);
  m_bndStateFluxStore.resize(m_nbrFaceFlxPnts);
  m_bndStateFlux.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_bndStateFluxStore[iFlx].resize(m_nbrEqs);
    m_bndStateFlux[iFlx] = &m_bndStateFluxStore[iFlx];
  }
  m_avgGrad.resize(m_nbrEqs);
  m_corrFctDiv.resize(m_nbrSolPnts);
  m_faceJacobVecs.resize(m_nbrFaceFlxPnts);
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntsLocalCoords[iFlx].resize(m_dim);
    m_flxPntCoords[iFlx].resize(m_dim);
    m_unitNormalFlxPnts[iFlx].resize(m_dim);
    m_cellFlx[iFlx].resize(m_nbrEqs);
    m_flxPntRiemannFlux[iFlx].resize(m_nbrEqs);
    m_faceJacobVecs[iFlx].resize(m_dim);
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_flxPntGhostGrads[iFlx].push_back(new RealVector(m_dim));
      m_cellGradFlxPnt[iFlx].push_back(new RealVector(m_dim));
    }
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_corrections[iSol].resize(m_nbrEqs);
    m_corrFctDiv[iSol].resize(m_allCellFlxPnts->size());
  }

  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
  {
    m_avgGrad[iVar] = new RealVector(m_dim);
  }

  // compute the divergence of the correction function
  m_corrFctComputer->computeDivCorrectionFunction(frLocalData[0],m_corrFctDiv);

  // data of the compact BR2 boundary face gradient
  m_dimList.resize(m_dim+m_ndimplus);
  m_cellFluxProjVects.resize(m_dim+m_ndimplus);
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_dimList[iDim].resize(m_nbrSolPnts);
    m_cellFluxProjVects[iDim].resize(m_nbrSolPnts);
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      m_dimList[iDim][iSolPnt] = iDim;
      m_cellFluxProjVects[iDim][iSolPnt].resize(m_dim);
    }
  }
  m_solJacobDet.resize(m_nbrSolPnts);
  m_projectedCorr.resize(m_dim);
  m_gradVarsSolPnts.resize(m_nbrEqs,m_nbrSolPnts);
  m_gradVarStatePtrs.resize(m_nbrSolPnts);

  m_compactGradsSolPnts.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_compactGradsSolPnts[iSol].resize(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_compactGradsSolPnts[iSol][iEq].resize(m_dim);
    }
  }

  m_flxPntGradVarsStore.resize(m_nbrFaceFlxPnts);
  m_bndGradVarsStore.resize(m_nbrFaceFlxPnts);
  m_flxPntGradVars.resize(m_nbrFaceFlxPnts);
  m_bndGradVars.resize(m_nbrFaceFlxPnts);
  m_bndGradFlxPnt.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntGradVarsStore[iFlx].resize(m_nbrEqs);
    m_bndGradVarsStore[iFlx].resize(m_nbrEqs);
    m_flxPntGradVars[iFlx] = &m_flxPntGradVarsStore[iFlx];
    m_bndGradVars[iFlx] = &m_bndGradVarsStore[iFlx];
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_bndGradFlxPnt[iFlx].push_back(new RealVector(m_dim));
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
 
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    deletePtr(m_cellStatesFlxPnt[iFlx]);
    deletePtr(m_flxPntGhostSol[iFlx]);
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      deletePtr(m_flxPntGhostGrads[iFlx][iGrad]);
      deletePtr(m_cellGradFlxPnt[iFlx][iGrad]); 
    }
    m_cellGradFlxPnt[iFlx].clear();
    m_flxPntGhostGrads[iFlx].clear();
  }
  for (CFuint iFlx = 0; iFlx < m_bndGradFlxPnt.size(); ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_bndGradFlxPnt[iFlx].size(); ++iGrad)
    {
      deletePtr(m_bndGradFlxPnt[iFlx][iGrad]);
    }
    m_bndGradFlxPnt[iFlx].clear();
  }
  m_bndGradFlxPnt.clear();
  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
  {
    deletePtr(m_avgGrad[iVar]); 
  }
  m_avgGrad.clear();
  m_cellStatesFlxPnt.clear();
  m_flxPntGhostSol.clear();
  m_cellGradFlxPnt.clear();
  m_flxPntGhostGrads.clear();

  // unsetup parent class
  FluxReconstructionSolverCom::unsetup();

}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
