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
  ConvBndFaceTermRHSFluxReconstructionProvider("ConvBndCorrectionsRHS");

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
  m_resUpdates(),
  m_waveSpeedUpd(),
  m_gradUpdates(),
  m_faceJacobVecAbsSizeFlxPnts(),
  m_faceJacobVecSizeFlxPnts(),
  m_cellStatesFlxPnt(),
  m_cellFlx(),
  m_nbrEqs(),
  m_flxPntsLocalCoords(),
  m_unitNormalFlxPnts(),
  m_corrections()
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

  // set the data needed to compute the face terms;
  setCorrectionsData();

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
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  // get number of solution points
  const CFuint nbrSolPnt = frLocalData[0]->getNbrOfSolPnts();
    
  // get number of flux points
  const CFuint nbrFlxPnt1D = (frLocalData[0]->getFlxPntsLocalCoord1D())->size();
  
  // get solution point local coordinates
  m_solPntsLocalCoords = frLocalData[0]->getSolPntsLocalCoords();
    
  // get flux point local coordinates
  m_flxPntsLocalCoords1D = frLocalData[0]->getFlxPntsLocalCoord1D();
    
  // get the face normals
  m_faceNormals = frLocalData[0]->getFaceNormals();
   
  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConn();
	  
  // get the face connectivity per orientation
  m_faceConnPerOrient = frLocalData[0]->getFaceConnPerOrient();
  
  // get the face integration coefficient
  m_faceIntegrationCoefs = frLocalData[0]->getFaceIntegrationCoefs();
  
  // get flux point mapped coordinate directions
  m_faceMappedCoordDir = frLocalData[0]->getFaceMappedCoordDir();
  
  SafePtr<vector< RealVector > > allFlxPnts = frLocalData[0]->getFlxPntsLocalCoords();
  
  m_flxPntsLocalCoords.resize(nbrFlxPnt1D);
  m_faceJacobVecAbsSizeFlxPnts.resize(nbrFlxPnt1D);
  

  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      //CFLog(VERBOSE,"start Idx: " << bndFacesStartIdxs[iTR][m_orient] << "\n");
      
        for (CFuint iFlx = 0; iFlx < nbrFlxPnt1D; ++iFlx)
        {
          m_flxPntsLocalCoords[iFlx].resize(m_dim);
          m_flxPntsLocalCoords[iFlx] = (*allFlxPnts)[(*m_faceFlxPntConn)[m_orient][iFlx]];
        }
      
      //CFLog(VERBOSE,"m_orient: " << m_orient << "\n");
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // loop over faces with this orientation
      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
	//CFLog(VERBOSE,"faceID: " << faceID << "\n");
        // build the face GeometricEntity
        geoData.idx = faceID;
        m_face = m_faceBuilder->buildGE();

        // get the neighbouring cell
        m_intCell = m_face->getNeighborGeo(0);
	
	CFLog(VERBOSE,"cellID: " << m_intCell->getID() << "\n");

        // get the states in the neighbouring cell
        m_cellStates = m_intCell->getStates();
	
	

        // if cell is parallel updatable, compute the correction flux
        if ((*m_cellStates)[0]->isParUpdatable())
        {
	  // set the face ID in the BCStateComputer
	  m_bcStateComputer->setFaceID(m_face->getID());
	  
	  // Resize vectors
          m_cellStatesFlxPnt.resize(0);
          m_cellStatesFlxPnt.resize(nbrFlxPnt1D);
          m_cellFlx.resize(0);
          m_cellFlx.resize(nbrFlxPnt1D);
          m_flxPntCoords.resize(0);
          m_flxPntCoords.resize(nbrFlxPnt1D);
      
          // get solution polynomial values at nodes
          vector< vector< CFreal > > solPolyValsAtNodes
	      = frLocalData[0]->getSolPolyValsAtNode(m_flxPntsLocalCoords);
	  
          RealVector normal = (*m_faceNormals)[m_orient];
	  //CFLog(VERBOSE,"normal: " << normal << "\n");
	  
	  vector<RealVector> normalList;
          normalList.resize(nbrFlxPnt1D);
	  vector<RealVector> flxCoords1D;
	  flxCoords1D.resize(nbrFlxPnt1D);
	  
          // compute flux point coordinates
          for (CFuint iFlx = 0; iFlx < nbrFlxPnt1D; ++iFlx)
          {
	    //CFLog(VERBOSE, "state in flx pnt = " << *(m_cellStatesFlxPnt[iFlx]->getData()) << "\n");
	    flxCoords1D[iFlx].resize(1);
	    m_flxPntCoords[iFlx].resize(m_dim);
	    flxCoords1D[iFlx][0] = (*m_flxPntsLocalCoords1D)[iFlx];
	    m_flxPntCoords[iFlx] = m_face->computeCoordFromMappedCoord(flxCoords1D[iFlx]);	
	    //CFLog(VERBOSE, "flx pnt coord = " << m_flxPntCoords[iFlx] << "\n");
	
	    normalList[iFlx].resize(normal.size());
	    normalList[iFlx] = normal;
          }
          
          // compute face Jacobian vectors
          vector< RealVector > faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(flxCoords1D);
	  
          // Loop over flux points to extrapolate the states to the flux points and get the 
          // discontinuous normal flux in the flux points and the Riemann flux
          for (CFuint iFlxPnt = 0; iFlxPnt < nbrFlxPnt1D; ++iFlxPnt)
          {
	    // get face Jacobian vector sizes in the flux points
	    //CFLog(VERBOSE,"Before Jacob vec size\n");
	    DataHandle< vector< CFreal > >
	      faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
            // get face Jacobian vector size
            m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[faceID][iFlxPnt];
	    //CFLog(VERBOSE,"After Jacob vec size\n");
	    
            // set face Jacobian vector size with sign depending on mapped coordinate direction
            m_faceJacobVecSizeFlxPnts[iFlxPnt] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient];

	    m_unitNormalFlxPnts[iFlxPnt].resize(m_dim);
	    
            // set unit normal vector
            m_unitNormalFlxPnts[iFlxPnt] = faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
	
	    m_cellFlx[iFlxPnt].resize(0);
	    m_cellFlx[iFlxPnt].resize(m_nbrEqs);
	    
	    State* tempVector = new State(solPolyValsAtNodes[iFlxPnt][0]*(*((*(m_cellStates))[0]->getData())));
	
	    m_cellStatesFlxPnt[iFlxPnt] = tempVector;
  
            for (CFuint iSol = 1; iSol < nbrSolPnt; ++iSol)
            {
	      //CFLog(VERBOSE,"cellStates: " << *((*m_cellStates)[iSol]->getData()) << "\n");
              *(m_cellStatesFlxPnt[iFlxPnt]) = (State) (*(m_cellStatesFlxPnt[iFlxPnt]->getData()) + 
							  solPolyValsAtNodes[iFlxPnt][iSol]*(*((*(m_cellStates))[iSol]->getData())));
            }
            
	    // compute the normal flux at the current flux point
	    m_updateVarSet->computePhysicalData(*(m_cellStatesFlxPnt[iFlxPnt]), m_pData);
	    m_cellFlx[iFlxPnt] = m_updateVarSet->getFlux()(m_pData,m_unitNormalFlxPnts[iFlxPnt]);
	    //CFLog(VERBOSE, "Flux in flx pnt = " << m_cellFlx[iFlxPnt] << "\n");
	
          }
      
          // compute ghost states
          m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
      
          m_flxPntRiemannFlux.resize(nbrFlxPnt1D);

          for (CFuint iFlxPnt = 0; iFlxPnt < nbrFlxPnt1D; ++iFlxPnt)
          {
	    //CFLog(VERBOSE, "Ghost state = " << *(m_flxPntGhostSol[iFlxPnt]->getData()) << "\n");
	    m_flxPntRiemannFlux[iFlxPnt].resize(m_nbrEqs);
	    m_flxPntRiemannFlux[iFlxPnt] = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[iFlxPnt]),
									      *(m_flxPntGhostSol[iFlxPnt]),
									      m_unitNormalFlxPnts[iFlxPnt]);
	    //CFLog(VERBOSE, "RiemannFlux = " << m_flxPntRiemannFlux[iFlxPnt] << "\n");
          }
      

          for (CFuint jFlxPnt = 0; jFlxPnt < nbrFlxPnt1D; ++jFlxPnt)
          {
            m_cellFlx[jFlxPnt] = (m_flxPntRiemannFlux[jFlxPnt] - m_cellFlx[jFlxPnt])*m_faceJacobVecSizeFlxPnts[jFlxPnt]/10;
	    //CFLog(VERBOSE,"delta flux = " << m_cellFlx[jFlxPnt] << "\n");
          }
      
          computeWaveSpeedUpdates(m_waveSpeedUpd);
      
          // update the wave speeds
          updateWaveSpeed();
	  
          //m_corrFctComputer->computeCorrectionFunction(frLocalData[m_iElemType],m_corrFct);
          m_corrFctComputer->computeDivCorrectionFunction(frLocalData[0],m_corrFctDiv);
      
          m_corrections.resize(0);
          m_corrections.resize(nbrSolPnt);
       
          for (CFuint iSolPnt = 0; iSolPnt < nbrSolPnt; ++iSolPnt)
          {
	    m_corrections[iSolPnt].resize(0);
	    m_corrections[iSolPnt].resize(m_nbrEqs);

	    for (CFuint iFlxPnt1D = 0; iFlxPnt1D < nbrFlxPnt1D; ++iFlxPnt1D)
	    {
	      RealVector currentCorrFactor = m_cellFlx[iFlxPnt1D];
	      cf_assert(currentCorrFactor.size() == m_nbrEqs);
	    
	      // Fill in the matrix
	      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
	      {
	        m_corrections[iSolPnt][iVar] -= currentCorrFactor[iVar] * m_corrFctDiv[iSolPnt][(*m_faceFlxPntConn)[m_orient][iFlxPnt1D]]; 
	      }
	      if(m_intCell->getID() == 323 || m_intCell->getID() == 72)
	      {
		//CFLog(VERBOSE, "FI-FD = " << currentCorrFactor << "\n");
		//CFLog(VERBOSE, "div h = " << m_corrFctDiv[iSolPnt][(*m_faceFlxPntConn)[m_orient][iFlxPnt1D]] << "\n");
		//CFLog(VERBOSE, "correction = " << m_corrections[iSolPnt] << "\n\n\n");
	      }
	    }
	  }
	  
	  // update the rhs
          updateRHS();
        }

        // release the face
        m_faceBuilder->releaseGE();
      
      }
    }
  }
// CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::setCorrectionsData()
{
//   // set the face term data in the boundary face term computer
//   m_bndCorrectionsComputer->setCorrectionsData();
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSFluxReconstruction::updateRHS()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();
  
  const CFuint nbrStates = m_cellStates->size();

  // update rhs
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      //CFLog(VERBOSE, "Update = " << m_divContFlx[iState][iVar] << "\n");
      rhs[resID+iVar] += resFactor*m_corrections[iState][iVar];
      //CFLog(VERBOSE,"-div F = " << m_divContFlx[iState][iVar] << "\n");
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

// cout << m_flxPntSol[iSide][iFlx]->getLocalID() << " " << flush;

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
  
  // number of flux points
  const CFuint nbrFlxPnts = frLocalData[0]->getFlxPntsLocalCoord1D()->size();

  // create internal and ghost states and extra variables
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    //m_cellStatesFlxPnt.push_back(new State());
    m_flxPntGhostSol.push_back(new State());
  }

  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    //m_cellStatesFlxPnt[iFlx]->setLocalID(iFlx);
    m_flxPntGhostSol[iFlx]->setLocalID(iFlx);
  }
  
  // resize m_faceJacobVecSizeFlxPnts
  m_faceJacobVecSizeFlxPnts.resize(nbrFlxPnts);
  
  // resize m_unitNormalFlxPnts
  m_unitNormalFlxPnts.resize(nbrFlxPnts);

  // dimensionality and number of equations
  m_dim   = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // resize the physical data temporary vector
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  convTerm->resizePhysicalData(m_pData);
  
  m_flxPntCoords.resize(nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx].resize(m_dim);
    
  }
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
