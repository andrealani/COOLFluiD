#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSFluxReconstruction.hh"
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

MethodCommandProvider< ConvBndCorrectionsRHSFluxReconstruction, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvBndCorrectionsRHSFluxReconstructionProvider("ConvBndCorrectionsRHS");

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSFluxReconstruction::ConvBndCorrectionsRHSFluxReconstruction(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_gradients("gradients"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  m_updateVarSet(CFNULL),
  m_faceBuilder(CFNULL),
  m_bcStateComputer(CFNULL),
  m_corrFctComputer(CFNULL),
  m_riemannFluxComputer(CFNULL),
  m_flxPntRiemannFlux(CFNULL),
  m_faceMappedCoordDir(CFNULL),
  m_flxPntCoords(),
  m_face(),
  m_intCell(),
  m_orient(),
  m_dim(),
  m_corrFctDiv(),
  m_cellStates(),
  m_waveSpeedUpd(),
  m_gradUpdates(),
  m_faceJacobVecAbsSizeFlxPnts(),
  m_faceJacobVecSizeFlxPnts(),
  m_cellStatesFlxPnt(),
  m_allCellFlxPnts(),
  m_cellFlx(),
  m_nbrEqs(),
  m_flxPntsLocalCoords(),
  m_unitNormalFlxPnts(),
  m_corrections(),
  m_nbrSolPnts(),
  m_nbrFaceFlxPnts()
{
}

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSFluxReconstruction::~ConvBndCorrectionsRHSFluxReconstruction()
{
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
ConvBndCorrectionsRHSFluxReconstruction::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_gradients);
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::executeOnTrs()
{
  CFAUTOTRACE;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current QuadFreeBCFluxReconstruction TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();
  
  CFLog(VERBOSE,"ConvBndCorrectionRHSFluxReconstruction::executeOnTRS: " << faceTrs->getName() << "\n");

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

        // get the neighbouring cell
        m_intCell = m_face->getNeighborGeo(0);
	
	// get the states in the neighbouring cell
        m_cellStates = m_intCell->getStates();
	
	//CFLog(VERBOSE,"faceID: " << faceID << ", real face ID: " << m_face->getID() << "\n");
	//CFLog(VERBOSE,"cellID: " << m_intCell->getID() << "\n");
	//CFLog(VERBOSE,"coord state 0: " << (((*m_cellStates)[0])->getCoordinates()) << "\n");
	//CFLog(VERBOSE,"state 0: " << *(((*m_cellStates)[0])->getData()) << "\n");
	
	setBndFaceData(m_face->getID());//faceID

        // if cell is parallel updatable, compute the correction flux
        if ((*m_cellStates)[0]->isParUpdatable())
        {
	  // compute the discontinuous flx in the flx pnts
	  computeDiscontinuousFlx();
	  
	  // set the face ID in the BCStateComputer
	  m_bcStateComputer->setFaceID(m_face->getID());
      
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
        } 
        
        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::computeDiscontinuousFlx()
{
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  // get solution polynomial values at nodes
  vector< vector< CFreal > > solPolyValsAtNodes
	= frLocalData[0]->getSolPolyValsAtNode(m_flxPntsLocalCoords);
	
  // Loop over flux points to extrapolate the states to the flux points and get the 
  // discontinuous normal flux in the flux points and the Riemann flux
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    State* tempVector = new State(solPolyValsAtNodes[iFlxPnt][0]*(*((*(m_cellStates))[0]->getData())));

    m_cellStatesFlxPnt[iFlxPnt] = tempVector;
  
    for (CFuint iSol = 1; iSol < m_nbrSolPnts; ++iSol)
    {
    CFLog(DEBUG_MIN,"cellStates: " << *((*m_cellStates)[iSol]->getData()) << "\n");
    *(m_cellStatesFlxPnt[iFlxPnt]) = (State) (*(m_cellStatesFlxPnt[iFlxPnt]->getData()) + 
			  solPolyValsAtNodes[iFlxPnt][iSol]*(*((*(m_cellStates))[iSol]->getData())));
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::computeInterfaceFlxCorrection()
{ 
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {     
    // compute the normal flux at the current flux point
    m_updateVarSet->computePhysicalData(*(m_cellStatesFlxPnt[iFlxPnt]), m_pData);
    m_cellFlx[iFlxPnt] = m_updateVarSet->getFlux()(m_pData,m_unitNormalFlxPnts[iFlxPnt]);
   }
      
   // compute ghost states
   m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);

   // compute the riemann flux minus the discontinuous flx in the flx pnts
   for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
   {
     m_flxPntRiemannFlux[iFlxPnt] = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[iFlxPnt]),
								      *(m_flxPntGhostSol[iFlxPnt]),
								      m_unitNormalFlxPnts[iFlxPnt]);
     m_cellFlx[iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt] - m_cellFlx[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt];
     if(m_intCell->getID() == 36)
     {
       CFLog(VERBOSE, "state = " << *(m_cellStatesFlxPnt[iFlxPnt]->getData()) << "\n");
       CFLog(VERBOSE, "Ghost state = " << *(m_flxPntGhostSol[iFlxPnt]->getData()) << "\n");
       CFLog(VERBOSE, "RiemannFlux = " << m_flxPntRiemannFlux[iFlxPnt] << "\n");
       CFLog(VERBOSE, "Flux in flx pnt = " << m_cellFlx[iFlxPnt] << "\n");
       CFLog(VERBOSE, "unit vector = " << m_unitNormalFlxPnts[iFlxPnt] << "\n");
       CFLog(VERBOSE, "delta flux = " << m_cellFlx[iFlxPnt] << "\n");
     }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::setBndFaceData(CFuint faceID)
{
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
    
  // compute flux point coordinates
  SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  
  // compute flux point coordinates
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx] = m_face->computeCoordFromMappedCoord((*flxLocalCoords)[iFlx]);	
  }
          
  // compute face Jacobian vectors
  vector< RealVector > faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*flxLocalCoords);
  
  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > >
  faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
  // Loop over flux points to extrapolate the states to the flux points and get the 
  // discontinuous normal flux in the flux points and the Riemann flux
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[faceID][iFlxPnt];
    
    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPnts[iFlxPnt] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient];
    
    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
    
    if(m_intCell->getID() == 36)
     {
       CFLog(VERBOSE, "jacob vec = " << faceJacobVecs[iFlxPnt] << "\n");
       CFLog(VERBOSE, "abs vec = " << m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] << "\n");
     }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::computeCorrection(vector< RealVector >& corrections)
{ 
  cf_assert(corrections.size() == m_nbrSolPnts);
  
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // reset the corrections
    corrections[iSolPnt] = 0.0;
    
    cf_assert(corrections[iSolPnt].size() == m_nbrEqs);

    // compute the term due to each flx pnt
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // the current correction factor
      RealVector currentCorrFactor = m_cellFlx[iFlxPnt];
      cf_assert(currentCorrFactor.size() == m_nbrEqs);
    
      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        corrections[iSolPnt][iVar] -= currentCorrFactor[iVar] * m_corrFctDiv[iSolPnt][(*m_faceFlxPntConn)[m_orient][iFlxPnt]]; 
      }
      if(m_intCell->getID() == 36)
      {
	CFLog(VERBOSE, "FI-FD = " << currentCorrFactor << "\n");
	CFLog(VERBOSE, "iSol: " << iSolPnt << ", flxID = " << (*m_faceFlxPntConn)[m_orient][iFlxPnt] << "\n");
	CFLog(VERBOSE, "div h = " << m_corrFctDiv[iSolPnt][(*m_faceFlxPntConn)[m_orient][iFlxPnt]] << "\n");
      }
    }
    if(m_intCell->getID() == 161)
    {
      CFLog(VERBOSE, "State: " << iSolPnt << ": " << (*m_cellStates)[iSolPnt]->getCoordinates() << "\n");
      CFLog(VERBOSE, "State global ID: " << (*m_cellStates)[iSolPnt]->getGlobalID() << "\n");
      CFLog(VERBOSE, "correction = " << corrections[iSolPnt] << "\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::updateRHS()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();
  
  // get nb of states
  const CFuint nbrStates = m_cellStates->size();

  // update rhs
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      rhs[resID+iVar] += resFactor*m_corrections[iState][iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::updateWaveSpeed()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbrSolPnts = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();
    updateCoeff[solID] += m_waveSpeedUpd;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::computeGradientCorrections()
{
//   // compute the face term contribution to the gradients
//   m_bndCorrectionsComputer->computeGradientCorrections(m_gradUpdates);
// 
//   // add updates to gradients
//   addGradBCTerms();
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::addGradBCTerms()
{
//   // get the gradients
//   DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
// 
//   // update the gradients
//   const CFuint nbrSolPnts = m_gradUpdates.size();
//   for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
//   {
//       // get state ID
//     const CFuint solID = (*m_cellStates)[iSol]->getLocalID();
// 
//       // update gradients
//     for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
//     {
//       gradients[solID][iGrad] += m_gradUpdates[iSol][iGrad];
//     }
//   }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::computeWaveSpeedUpdates(CFreal& waveSpeedUpd)
{
  waveSpeedUpd = 0.0;
  for (CFuint iFlx = 0; iFlx < m_cellFlx.size(); ++iFlx)
  {
    const CFreal jacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx];
				   
    // transform update states to physical data to calculate eigenvalues
    m_updateVarSet->computePhysicalData(*(m_cellStatesFlxPnt[iFlx]), m_pData);
    waveSpeedUpd += jacobXIntCoef*
                    m_updateVarSet->getMaxAbsEigenValue(m_pData,m_unitNormalFlxPnts[iFlx]);
    }
}


//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::setup()
{
  CFAUTOTRACE;

  FluxReconstructionSolverCom::setup();
  
  // get cell builder
  m_faceBuilder = getMethodData().getFaceBuilder();
  
  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();
  
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
  
  // number of sol points
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();
  
  // get solution point local coordinates
  m_solPntsLocalCoords = frLocalData[0]->getSolPntsLocalCoords();
    
  // get flux point local coordinates
  m_flxPntsLocalCoords1D = frLocalData[0]->getFlxPntsLocalCoord1D();
   
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

  // create internal and ghost states and extra variables
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntGhostSol.push_back(new State());
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntGhostSol[iFlx]->setLocalID(iFlx);
  }
  
  // resize m_faceJacobVecSizeFlxPnts
  m_faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  
  // resize m_unitNormalFlxPnts
  m_unitNormalFlxPnts.resize(m_nbrFaceFlxPnts);

  // dimensionality and number of equations
  m_dim   = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // resize the physical data temporary vector
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  convTerm->resizePhysicalData(m_pData);
  
  // resize vectors
  m_flxPntsLocalCoords.resize(m_nbrFaceFlxPnts);
  m_faceJacobVecAbsSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_cellStatesFlxPnt.resize(m_nbrFaceFlxPnts);
  m_cellFlx.resize(m_nbrFaceFlxPnts);
  m_flxPntCoords.resize(m_nbrFaceFlxPnts);
  m_flxPntCoords.resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFlux.resize(m_nbrFaceFlxPnts);
  m_corrections.resize(m_nbrSolPnts);
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntsLocalCoords[iFlx].resize(m_dim);
    m_flxPntCoords[iFlx].resize(m_dim);
    m_unitNormalFlxPnts[iFlx].resize(m_dim);
    m_cellFlx[iFlx].resize(m_nbrEqs);
    m_flxPntRiemannFlux[iFlx].resize(m_nbrEqs);
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_corrections[iSol].resize(m_nbrEqs);
  }
  
  // compute the divergence of the correction function
  m_corrFctComputer->computeDivCorrectionFunction(frLocalData[0],m_corrFctDiv);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::unsetup()
{
  CFAUTOTRACE;

  FluxReconstructionSolverCom::unsetup();

}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
