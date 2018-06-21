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
  socket_gradientsAV("gradientsAV"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  m_updateVarSet(CFNULL),
  m_faceBuilder(CFNULL),
  m_bcStateComputer(CFNULL),
  m_face(),
  m_intCell(),
  m_orient(),
  m_dim(),
  m_cellStates(),
  m_waveSpeedUpd(),
  m_nbrEqs(),
  m_nbrSolPnts(),
  m_nbrFaceFlxPnts(),
  m_solPntsLocalCoords(CFNULL),
  m_flxPntsLocalCoords(),
  m_allCellFlxPnts(CFNULL),
  m_flxPntCoords(),
  m_faceFlxPntConn(CFNULL),
  m_faceConnPerOrient(CFNULL),
  m_faceIntegrationCoefs(CFNULL),
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
  m_flxSolDep(CFNULL),
  m_nbrSolDep()
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
  result.push_back(&socket_gradientsAV);
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

  // get current bnd face TRS
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
  
  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm() || getMethodData().hasArtificialViscosity();
  
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

	CFLog(VERBOSE,"cellID: " << m_intCell->getID() << "\n");
//	if (m_intCell->getID() == 72)
//	{
//	  //CFLog(VERBOSE,"coord state: " << (((*m_cellStates)[0])->getCoordinates()) << "\n");
//	  CFLog(VERBOSE,"face ID: " << m_face->getID() << "\n");
//	}

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
          computeWaveSpeedUpdates(m_waveSpeedUpd);

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
          if(m_intCell->getID() == 35)
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

void ConvBndCorrectionsRHSFluxReconstruction::computeInterfaceFlxCorrection()
{ 
  // compute the riemann flux in the flx pnts
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    m_flxPntRiemannFlux[iFlxPnt] = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[iFlxPnt]),
								      *(m_flxPntGhostSol[iFlxPnt]),
								      m_unitNormalFlxPnts[iFlxPnt]);
    
    // store the local Riemann flux, scaled with geometrical Jacobian
    m_flxPntRiemannFlux[iFlxPnt] *= m_faceJacobVecSizeFlxPnts[iFlxPnt];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::computeFlxPntStates()
{ 
  // Loop over flux points to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // reset the extrapolated states
    *(m_cellStatesFlxPnt[iFlxPnt]) = 0.0;
    
    // get current flx pnt idx
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    
    // extrapolate the states to current flx pnt
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

void ConvBndCorrectionsRHSFluxReconstruction::setBndFaceData(CFuint faceID)
{  
  // compute flux point coordinates if needed for the BC
  if (m_bcStateComputer->needsSpatialCoordinates())
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      m_flxPntCoords[iFlx] = m_face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
    }
  }
          
  // compute face Jacobian vectors
  vector< RealVector > faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
  
  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
  // Loop over flux points to set the unit normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[faceID][iFlxPnt];
    
    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPnts[iFlxPnt] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient];
    
    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
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
  }
  
  // loop over flx and sol pnts to compute -divhFI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];

    // the current correction factor previously computed
    const RealVector currentCorrFactor = m_flxPntRiemannFlux[iFlxPnt];

    cf_assert(currentCorrFactor.size() == m_nbrEqs);

    // compute the term due to each flx pnt
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

      // divergence of the correctionfct
      const CFreal divh = m_corrFctDiv[solIdx][flxIdx];

      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        corrections[solIdx][iVar] -= currentCorrFactor[iVar] * divh; 
      }
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

void ConvBndCorrectionsRHSFluxReconstruction::updateWaveSpeed()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // loop over sol pnts to add the wave speed updates to the socket
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();
    updateCoeff[solID] += m_waveSpeedUpd;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::computeWaveSpeedUpdates(CFreal& waveSpeedUpd)
{
  // reset the wave speed update
  waveSpeedUpd = 0.0;
  
  // loop over flux pnts to compute the updates
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFreal jacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx];
				   
				   
    // transform update states to physical data to calculate eigenvalues
    m_updateVarSet->computePhysicalData(*(m_cellStatesFlxPnt[iFlx]), m_pData);
    waveSpeedUpd += jacobXIntCoef*m_updateVarSet->getMaxAbsEigenValue(m_pData,m_unitNormalFlxPnts[iFlx]);
    //if (waveSpeedUpd > 10.0) CFLog(INFO, "wvspConvBnd: " << waveSpeedUpd << "\n");
  }
  
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::computeGradientBndFaceCorrections()
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
      const RealVector projectedCorr = (avgSol-(*m_cellStatesFlxPnt[iFlx])[iEq])*m_faceJacobVecSizeFlxPnts[iFlx]*m_unitNormalFlxPnts[iFlx];

      // Loop over solution pnts to calculate the grad updates
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
      {
        const CFuint iSolIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

	/// @todo Check if this is also OK for triangles!!
	m_gradUpdates[iSolIdx][iEq] += projectedCorr*m_corrFctDiv[iSolIdx][flxIdx];
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

void ConvBndCorrectionsRHSFluxReconstruction::setup()
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

  m_flxSolDep = frLocalData[0]->getFlxPntSolDependency();

  m_nbrSolDep = ((*m_flxSolDep)[0]).size();

  // create internal and ghost states
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntGhostSol.push_back(new State());
    m_cellStatesFlxPnt.push_back(new State());
  }

  // set an ID as initialization
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntGhostSol[iFlx]->setLocalID(iFlx);
    m_cellStatesFlxPnt[iFlx]->setLocalID(iFlx);
  }

  // dimensionality and number of equations
  m_dim   = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

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
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  { 
    m_flxPntsLocalCoords[iFlx].resize(m_dim);
    m_flxPntCoords[iFlx].resize(m_dim);
    m_unitNormalFlxPnts[iFlx].resize(m_dim);
    m_flxPntRiemannFlux[iFlx].resize(m_nbrEqs);
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
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    deletePtr(m_cellStatesFlxPnt[iFlx]);
    deletePtr(m_flxPntGhostSol[iFlx]);
  }
  m_cellStatesFlxPnt.clear();
  m_flxPntGhostSol.clear();

  // unsetup parent class
  FluxReconstructionSolverCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
