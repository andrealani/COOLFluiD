//Points to remember

//To print the value of Bx in flux points, use the line CFLog(VERBOSE,"Bx = "<<(*(ghostStates[iState]))[0]<<"\n"); in the function ComputeGhostStates(). This is done to check whether the flux points are being matched correctly.
//To check whether faces are being matched correctly, print the coordinates of the faces that are being matched. For the faces that being matched in the y direction, the x-coordinate should be the same and the y coordinates should have the opposite signs. This is done in the line

//Added from FV




#include "Framework/BadFormatException.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/BCPeriodic.hh"
//Added from ConvBnd
//#include "Framework/BaseTerm.hh"
//Added from FV
#include "Common/MPI/MPIStructDef.hh"
#include "MathTools/MathChecks.hh"
#include <iostream>
#include <map>
#include <algorithm>
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCPeriodic,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionModule >
  BCPeriodicProvider("Periodic");

//////////////////////////////////////////////////////////////////////////////

BCPeriodic::BCPeriodic(const std::string& name) :
  BCStateComputer(name),
        //from ConvBnd to include in preProcess done for parallelization (Copying from FV)
        m_flxPntsLocalCoords(),
        m_allCellFlxPnts(),
        m_cellStatesFlxPnt(),
        m_faceBuilder(),
        m_faceMappedCoordDir(),
        m_faceJacobVecs(),
        m_faceJacobVecAbsSizeFlxPnts(),
        m_faceJacobVecSizeFlxPnts(),
        m_unitNormalFlxPnts(),
        //till here
  m_thisTRS(),
  m_flxLocalCoords(), 
  m_flxPntCoords(),
  m_nbrFaceFlxPnts(),
  m_dim(),
  m_orient(),      
  m_nbrSolDep(),
  m_solPolyValsAtFlxPnts(),
  m_otherFace(),
  m_intCell(),
  m_cellStates(),
  m_cellGrads(),
  _faceConnectivityMap(),
  _fluxPointConnectivityMap(), 
  _orientMap(),      
  m_flxSolDep(),
  m_faceFlxPntConn(),
  _faceBuilder(),      
  socket_gradients("gradients"), 
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),      
  _globalToLocalTRSFaceID(),
  _localWestEastMap(),
//  socket_states("states"),
//  socket_gstates("gstates"),
//  socket_nodes("nodes"),
  _translationVector(),
  _localConnectivityMap()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);


  _threshold = 1e-5; 
  setParameter("Threshold",&_threshold);

 // _translationVector = std::vector<CFreal>();
  setParameter("TranslationVector",&_translationVector);

}

//////////////////////////////////////////////////////////////////////////////

BCPeriodic::~BCPeriodic()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::computeGhostStates(const vector< State* >& intStates,
					  vector< State* >& ghostStates,
					  const std::vector< RealVector >& normals,
					  const std::vector< RealVector >& coords)
{



//CFLog(VERBOSE,"here2\n");
  //CFLog(VERBOSE, "BCPeriodic - Computing ghost states.. \n");
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbrStates = intStates.size();
  cf_assert(nbrStates == ghostStates.size());
  cf_assert(nbrStates == normals.size());

  // we prepare the face builder
  FaceToCellGEBuilder::GeoData& faceData = _faceBuilder->getDataGE();
  faceData.cellsTRS = cellTrs;
  faceData.facesTRS = m_thisTRS; 
  faceData.isBoundary = true;
  //CFLog(VERBOSE, "BCPeriodic - Preparing face builder.. \n");
  //cout <<"HERE 1 face id: "<< m_face->getID()<<endl;

  // Get the localFaceID from the map, knowing the faceGlobalID
  const CFuint faceLocalID = _globalToLocalTRSFaceID.find(m_face->getID());
//CFLog(VERBOSE,"here3\n");
  // If we assume that all flux point of this face are connected to same face at the other side, we
  // can do this outside the flux points loop
  //cout <<"HERE 2 "<<endl;

  const CFuint otherFaceLocalID = _faceConnectivityMap[faceLocalID*nbrStates];
  //cout <<"HERE 3 "<<endl;

 //CFLog(VERBOSE,"here3.1\n");
  //CFLog(VERBOSE, "BCPeriodic - We are studying face "<<m_face->getID()<<" with faceLocalID "<<faceLocalID<<" linked to "<<otherFaceLocalID<<"\n");
  // If we take this assumption further, we could store the _faceConnectivityMap only for the faces,
  // instead that for each flux point!  
//CFLog(VERBOSE,"here3.2\n");
  // We build the geometric entity of that face
  faceData.idx = otherFaceLocalID;
  //CFLog(VERBOSE,"here3.3\n");
  m_otherFace = _faceBuilder->buildGE();
  //const CFuint otherFaceGlobalID = m_otherFace->getID();
//CFLog(VERBOSE,"here4\n");
  if(_nbProcesses == 1)
{
  // loop over the flux points
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
     // We compute the local index of the flux point
     CFuint localFlxPntID = faceLocalID*nbrStates + iState; 

     // get the neighbouring cell
     m_intCell = m_otherFace->getNeighborGeo(0);

     // get the states in the neighbouring cell
     m_cellStates = m_intCell->getStates();

     // We get assigned flux point
     CFuint iFlxPnt = _fluxPointConnectivityMap[localFlxPntID];
     //CFLog(VERBOSE, "BCPeriodic - Computing BC at FlxPnt "<<iState<< ", linked to "<<iFlxPnt<<"\n");

     // reset the extrapolated states
     *(ghostStates[iState]) = 0.0;
    
     // get current flx pnt idx
     const CFuint currFlxIdx = (*m_faceFlxPntConn)[_orientMap[localFlxPntID]][iFlxPnt];
    
     // extrapolate the states to current flx pnt
     for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
     {
        const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];

        *(ghostStates[iState]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(*((*m_cellStates)[solIdx]));
     }
     
	//CFLog(VERBOSE,"Bx = "<<(*(ghostStates[iState]))[0]<<"\n");
  } 
//CFLog(VERBOSE,"here5\n");

  // release the face
  _faceBuilder->releaseGE();
  //CFLog(VERBOSE, "BCPeriodic - Computing BC at FlxPnt END \n");

/*
  
  CFreal vn = 0.;
  CFreal area2 = 0.;
  
  cf_assert(m_velocityIDs.size() == dim);
  for (CFuint iState = 0; iState < intStates.size() ; ++iState){
   for (CFuint i = 0; i < m_velocityIDs.size(); ++i) { 
     const CFreal nxComp = normals[iState][i];
     vn += nxComp*(*intStates[iState])[m_velocityIDs[i]];
     area2 += nxComp*nxComp;
   }
   
   CFuint jxx = 0;
   for (CFuint i = 0; i < m_isVelocityComp.size(); ++i) {
     if (!m_isVelocityComp[i]) {
       (*ghostStates[iState])[i] = (*intStates[iState])[i];
     }
     else {
       (*ghostStates[iState])[i]= (*intStates[iState])[i] - 2.0*vn*normals[iState][jxx]/area2;
       jxx++;
     }
     // if (i == 2 &&  (*ghostStates[iState])[i] < -0.000001) CFLog(VERBOSE, "intState: " << *intStates[iState] << ", ghost: " << *ghostStates[iState] << "\n");
   }
   
   //CFLog(DEBUG_MAX, "Periodic::setGhostState() => ghostState = " << *ghostState << "\n"); 
  }
*/
}//end of nbProcesses == 1 loop

    else{
  for(CFuint iState = 0;iState < nbrStates; ++iState){
       CFuint localFlxPntID = faceLocalID*nbrStates + iState;
       CFLog(VERBOSE,"faceLocalID = "<<faceLocalID<<"\n");
       CFLog(VERBOSE,"localFlxPntID = "<<localFlxPntID<<"\n");
       //CFLog(VERBOSE,"else here 1\n");
         const CFuint j = _localConnectivityMap[localFlxPntID].size() - 1;
         CFLog(VERBOSE,"j = "<<j<<"\n");
         //CFLog(VERBOSE,"else here 1.1\n");
  const CFuint periodicFaceID = _localConnectivityMap[localFlxPntID][j].getFaceID();
  CFLog(VERBOSE,"otherfacelocalID = "<<otherFaceLocalID+iState<<" periodicFaceID = "<<periodicFaceID<<"\n");
  //CFLog(VERBOSE,"else here 1.2\n");
  cf_assert(localFlxPntID < _localConnectivityMap.size());
  //CFLog(VERBOSE,"else here 1.3\n");
    cf_assert(j < _localConnectivityMap[localFlxPntID].size()); 
    
    //CFLog(VERBOSE,"else here 2\n");
    
    CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicState() => 1 end\n");
    const CFuint periodicProcess = _localConnectivityMap[localFlxPntID][j].getProcess();
    for(CFuint h=0; h<NbEqs; h++){
      cf_assert(periodicProcess < _recvdispls.size());
      CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicState() => 2 start\n");
      const CFuint idx = _recvdispls[periodicProcess] + NbEqs*periodicFaceID + h;
      cf_assert(idx < _recvbuf.size());
      (*_periodicState)[h] = _recvbuf[idx];
      CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicState() => 2 end\n");
      CFLog(VERBOSE,"recvbuf = "<<_recvbuf[idx]<<"\n");
    }//CFLog(VERBOSE,"else here 3\n");
    CFLog(DEBUG_MIN, "BCPeriodic::computePeriodicState() => end\n");
     CFLog(DEBUG_MIN, "BCPeriodic::setGhostState() => start\n");
  *(ghostStates[iState]) = 0.0;
  //CFLog(VERBOSE,"else here 4\n");
  State *const periodicState = _periodicState;
  //CFLog(VERBOSE,"else here 5\n");
    *(ghostStates[iState]) = (*periodicState);
    CFLog(VERBOSE,"ghost state = "<<*(ghostStates[iState])<<"\n");
 CFLog(VERBOSE,"else here 6\n");
  CFLog(DEBUG_MIN, "BCPeriodic::setGhostState() => end\n");
    
      }
  // release the face
  _faceBuilder->releaseGE();
  
  }
//  CFLog(VERBOSE,"End here\n");


}

//if (nbProcesses > 1)
//{

   

//}
//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::computeGhostGradients
(const std::vector< std::vector< RealVector* > >& intGrads,
 std::vector< std::vector< RealVector* > >& ghostGrads,
 const std::vector< RealVector >& normals,
 const std::vector< RealVector >& coords)
{

  //CFLog(VERBOSE, "BCPeriodic - Computing ghost gradients.. CAUTION! This has not been tested yet!\n");
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  const CFuint nbrGradVars = intGrads[0].size();
  const CFuint nbrStateGrads = intGrads.size();
  const CFuint nbrStates = m_cellStates->size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // we prepare the face builder
  FaceToCellGEBuilder::GeoData& faceData = _faceBuilder->getDataGE();
  faceData.cellsTRS = cellTrs;
  faceData.facesTRS = m_thisTRS; //store both references as protecteed values
  faceData.isBoundary = true;
  //CFLog(VERBOSE, "BCPeriodic - (grads) Preparing face builder.. \n");
//cout <<"HERE 1 "<<endl;
  // Get the localFaceID from the map, knowing the faceGlobalID
  const CFuint faceLocalID = _globalToLocalTRSFaceID.find(m_face->getID());
//cout <<"HERE 2 "<<endl;

  // If we assume that all flux point of this face are connected to same face at the other side, we
  // can do this outside the flux points loop
  const CFuint otherFaceLocalID = _faceConnectivityMap[faceLocalID*nbrStateGrads];
 
 // CFLog(VERBOSE, "BCPeriodic - (grads) We are studying face "<<m_face->getID()<<" with faceLocalID "<<faceLocalID<<" linked to "<<otherFaceLocalID<<"\n");
  // If we take this assumption further, we could store the _faceConnectivityMap only for the faces,
  // instead that for each flux point!  

  // We build the geometric entity of that face
  faceData.idx = otherFaceLocalID;
  m_otherFace = _faceBuilder->buildGE();
  //const CFuint otherFaceGlobalID = m_otherFace->getID();

  // loop over the flux points
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
     // We compute the local index of the flux point
     CFuint localFlxPntID = faceLocalID*nbrStateGrads + iState; 

     // get the neighbouring cell
     m_intCell = m_otherFace->getNeighborGeo(0);

     // get the states in the neighbouring cell
     m_cellStates = m_intCell->getStates();

     // We set the gradients (the resize could be done in setup? In DiffBnd is done at setBndFaceData)
     m_cellGrads.resize(nbrStates);
     for (CFuint iState = 0; iState < nbrStates; ++iState)
     {
        const CFuint stateID = (*m_cellStates)[iState]->getLocalID();
        m_cellGrads[iState] = &gradients[stateID];
     }

     // We get assigned flux point
     CFuint iFlxPnt = _fluxPointConnectivityMap[localFlxPntID];
   //  CFLog(VERBOSE, "BCPeriodic - (grads) Computing BC at FlxPnt "<<iState<< ", linked to "<<iFlxPnt<<"\n");

     // reset the grads in flx pnts
     for (CFuint iVar = 0; iVar < nbrGradVars; ++iVar)
     {
        *(ghostGrads[iFlxPnt][iVar]) = 0.0;
     }
  
     // index of current flx pnt
     const CFuint currFlxIdx = (*m_faceFlxPntConn)[_orientMap[localFlxPntID]][iFlxPnt];
    
     // Loop over sol points to add the contributions to each sol pnt
     for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
     {
        const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];
      
        for (CFuint iVar = 0; iVar < nbrGradVars; ++iVar)
        {
           *(ghostGrads[iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*((*(m_cellGrads[solIdx]))[iVar]);
        }
     }

  } 
  // release the face
  _faceBuilder->releaseGE();
 // CFLog(VERBOSE, "BCPeriodic - (grads) Computing BC at FlxPnt END \n");

//cout <<"HERE 3 "<<endl;


/*
  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();

  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // number of gradient variables
  cf_assert(nbrStateGrads > 0);

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState) {
    // normal
    const RealVector& normal = normals[iState];

    // tangential unit vector
    m_tangent[XX] = -normal[YY];
    m_tangent[YY] =  normal[XX];

    const vector< RealVector* >& velocityGradI = intGrads[iState];
    m_velocityNGradI = 0.;
    m_velocityTGradI = 0.;
    CFuint jxx = 0;
    for (CFuint i = 0; i < m_isVelocityComp.size(); ++i) {
      if (!m_isVelocityComp[i]) {
	const RealVector& varGradI =  *intGrads[iState][i];
	RealVector& varGradG =  *ghostGrads[iState][i];
	const CFreal nVarGrad = MathTools::MathFunctions::innerProd(varGradI, normal);
	varGradG = varGradI - 2.0*nVarGrad*normal;
      }
      else {
	// internal normal and tangential component
	m_velocityNGradI +=  *velocityGradI[i]*normal[jxx];
	m_velocityTGradI +=  *velocityGradI[i]*m_tangent[jxx];
	++jxx;
      }
    }

    m_velocityNGradG = m_velocityNGradI;

    jxx = 0;
    CFreal nGradUT = 0.;
    for (CFuint i = 0; i < m_isVelocityComp.size(); ++i) {
      if (m_isVelocityComp[i]){
        // ghost normal and tangential component
        nGradUT += m_velocityTGradI[jxx]*normal[jxx];  
	++jxx;
      }
    }

    m_velocityTGradG = m_velocityTGradI - 2.0*nGradUT*normal;

    // project onto x- and y-axis
    jxx = 0;
    for (CFuint i = 0; i < m_isVelocityComp.size(); ++i) {
      if (m_isVelocityComp[i]) {
        *ghostGrads[iState][i] = m_velocityNGradG*normal[jxx] + m_velocityTGradG*m_tangent[jxx];
	++jxx;
      }
    }
  }
*/
}

//////////////////////////////////////////////////////////////////////////////


void BCPeriodic::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Threshold","Tolerance for considering coordinates matching.");
  options.addConfigOption< std::vector<CFreal> >("TranslationVector","TranslationVector between planes");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > BCPeriodic::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  // Add the needed DataSocketSinks
//  result.push_back(&socket_normals);
//  result.push_back(&socket_states);
//  result.push_back(&socket_gstates);
//  result.push_back(&socket_nodes);
  result.push_back(&socket_gradients);
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);
  return result;
}
             
//////////////////////////////////////////////////////////////////////////////


void BCPeriodic::createFaceOrientationStartIndexes()
{
  CFAUTOTRACE;
   MPI_Barrier(_comm);//the processing of a an individual processor is paused temporarily till all the processors compute till some certain values
  //CFLog(VERBOSE,"here\n");
  MPI_Allgather(&nbGeoEnts, 1, MPIStructDef::getMPIType(&nbGeoEnts), 
	&_nbFacesPerProcess[0], 1, MPIStructDef::getMPIType(&nbGeoEnts), _comm);

  CFLog(VERBOSE, "\n\n createFaceOrientationStartIndexes() \n\n");

  // START INDEXES FOR INNER CFGeoEnt::CELLS
  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get connectivity where the face start indexes are stored
  SafePtr< ConnectivityTable< CFuint > > innerFacesStartIdxsConnTable
    = MeshDataStack::getActive()->getConnectivity("innerFacesStartIdxs");

  // get number of orientations
  const CFuint nbrFaceOrientsPlus1 = innerFacesStartIdxsConnTable->nbRows();

  // resize innerFacesStartIdxs
  innerFacesStartIdxs.resize(nbrFaceOrientsPlus1);

  // put indexes in innerFacesStartIdxs
  for (CFuint iIdx = 0; iIdx < nbrFaceOrientsPlus1; ++iIdx)
  {
    innerFacesStartIdxs[iIdx] = (*innerFacesStartIdxsConnTable)(iIdx,0);
  }

  // START INDEXES FOR BOUNDARY CFGeoEnt::CELLS
  // get bndFacesStartIdxs from FluxReconstructionSolverData
  map< std::string , vector< vector< CFuint > > >& bndFacesStartIdxs = getMethodData().getBndFacesStartIdxs();

  // get TRS list
  vector<Common::SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();

  // number of TRSs
  const CFuint nbTRSs = trsList.size();

  // loop over boundary TRSs
  CFuint iBndTRS = 0;
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    const std::string trsName = trsList[iTRS]->getName();

    if (trsList[iTRS]->hasTag("boundary"))
    {

      // get connectivity where the face start indexes are stored
      SafePtr< ConnectivityTable< CFuint > > boundaryFacesStartIdxsConnTable
        = MeshDataStack::getActive()->getConnectivity(trsName+"boundaryFacesStartIdxs");

      // get number of TRs in this TRS
      const CFuint nbTRs = trsList[iTRS]->getNbTRs();
      // number of boundary face orientations + 1
      const CFuint nbBndFaceOrientP1 = boundaryFacesStartIdxsConnTable->nbRows()/nbTRs;

      // array over TRs and startIdxs
      vector< vector< CFuint > > trStartIdxs(nbTRs);

      // loop over TRs
      CFuint iStartIdx = 0;
      for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
      {
        // resize trStartIdxs[iTR]
        trStartIdxs[iTR].resize(nbBndFaceOrientP1);

        // loop over start indexes
        for (CFuint iIdx = 0; iIdx < nbBndFaceOrientP1; ++iIdx, ++iStartIdx)
        {
          trStartIdxs[iTR][iIdx] = (*boundaryFacesStartIdxsConnTable)(iStartIdx,0);
        }
      }
      cf_assert(iStartIdx == boundaryFacesStartIdxsConnTable->nbRows());

      // store trStartIdxs in boundary TRS map
      bndFacesStartIdxs[trsName] = trStartIdxs;

      // increment boundary TRS counter
      ++iBndTRS;
    } 
  }

}



//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::setup()
{
  CFAUTOTRACE;
  CFLog(INFO, "BCPeriodic::setup \n");

/*
  addConfigOptionsTo(this);

  _threshold = 1e-5;
  setParameter("Threshold",&_threshold);
  setParameter("TranslationVector",&_translationVector);
*/


  // setup of the parent class
  BCStateComputer::setup();
  
  // no flux point coordinates required
  // Set this to false once the BC is working
  m_needsSpatCoord = false;// true; //
  
  _faceBuilder = getMethodData().getSecondFaceBuilder();
  


  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  CFLog(VERBOSE, "BCPeriodic::setup 2 \n");

  // get the face local coords of the flux points on one face
  m_flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  m_nbrFaceFlxPnts = m_flxLocalCoords->size();
  m_dim   = PhysicalModelStack::getActive()->getDim();
CFLog(VERBOSE, "BCPeriodic::here 1 \n");

  // get the face - flx pnt connectivity per orient (ConvBnd) 
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConn();
  m_flxSolDep = frLocalData[0]->getFlxPntSolDependency();
  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();

  m_nbrSolDep = ((*m_flxSolDep)[0]).size();
CFLog(VERBOSE, "BCPeriodic::here 2 \n");
  m_flxPntCoords.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx].resize(m_dim);
  }
  CFLog(VERBOSE, "BCPeriodic::setup 3 - "<<m_dim<<" "<<m_nbrFaceFlxPnts<<"\n");

  //from ConvBnd to include in preProcess done for parallelization (Copying from FV)
  
  m_faceBuilder = getMethodData().getFaceBuilder();
  
  //CFLog(VERBOSE,"here 1");
  // create internal and ghost states
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    
    m_cellStatesFlxPnt.push_back(new State());
  }

  // set an ID as initialization
  //CFLog(VERBOSE,"here 2");
   
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    
    m_cellStatesFlxPnt[iFlx]->setLocalID(iFlx);
  }
  //CFLog(VERBOSE,"here 3");
  
  m_flxPntsLocalCoords.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  { 
    m_flxPntsLocalCoords[iFlx].resize(m_dim);
    
  }
  
  //CFLog(VERBOSE,"here 4");
  // get all flux points of a cell
  m_allCellFlxPnts = frLocalData[0]->getFlxPntsLocalCoords();
  
  // get flux point mapped coordinate directions
  m_faceMappedCoordDir = frLocalData[0]->getFaceMappedCoordDir();
  
  // resize vectors
  m_unitNormalFlxPnts.resize(m_nbrFaceFlxPnts);
  m_faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_faceJacobVecAbsSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_faceJacobVecs.resize(m_nbrFaceFlxPnts);
  
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  { 
   // m_flxPntsLocalCoords[iFlx].resize(m_dim);
    m_unitNormalFlxPnts[iFlx].resize(m_dim);
    m_faceJacobVecs[iFlx].resize(m_dim);
  }
//till here
  
  //CFLog(VERBOSE,"here 5");
  /* Check if translation vector has right dimension */
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_assert(_translationVector.size()==dim);

  /* Convert std::vector<CFreal> to RealVector version */
  RealVector translationVector(dim);
  RealVector backTranslationVector(dim);
  CFLog(VERBOSE, "BCPeriodic::setup 4 dim "<<dim<<" \n");
  for(CFuint iDim=0; iDim<dim; ++iDim) {
    translationVector[iDim] = _translationVector[iDim];
    backTranslationVector[iDim] = -translationVector[iDim];
  }
  CFLog(VERBOSE, " ==> translationVector = " << translationVector << " \n");

  // get TRS list
  vector< SafePtr< TopologicalRegionSet > > trsList = MeshDataStack::getActive()->getTrsList();

//  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

//  IA: The BC object is created for each BC, so no need to do this!!
  // number of boundary TRSs with this BC
  //const CFuint nbrBCTRSs = m_trsNames.size();
  //_localConnectivityMap.resize(nbrBCTRSs);
  // get boundary TRSs
  //vector< SafePtr< TopologicalRegionSet > > bcTRSs(nbrBCTRSs);
  const CFuint nbTRSs = trsList.size();
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    CFLog(VERBOSE, "iTRS "<< iTRS << " \n");
    if (m_trsNames[0]==trsList[iTRS]->getName() ){
       m_thisTRS = trsList[iTRS];
       CFLog(VERBOSE, "Matching BC "<<m_trsNames[0]<<" with "<<m_thisTRS->getName() << "\n");
    }
   /*
    for (CFuint iBCTRS = 0; iBCTRS < nbrBCTRSs; ++iBCTRS)
    {
      if (m_trsNames[iBCTRS] == trsList[iTRS]->getName())
      {
        if (bcTRSs[iBCTRS].isNull())
        {
          bcTRSs[iBCTRS] = trsList[iTRS];
        }
        else
        {
          throw BadFormatException (FromHere(),"Two TRSs with the same name found!");
        }
      }
    }
   */
  }
  CFLog(VERBOSE,"End of TRS loop \n");

  // Loop over the boundary TRS
  //for (CFuint iBCTRS=0; iBCTRS < nbrBCTRSs; iBCTRS++){

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& faceData = _faceBuilder->getDataGE();
  faceData.cellsTRS = cellTrs;
  faceData.facesTRS = m_thisTRS;
  faceData.isBoundary = true;
//
//     _faceBuilder.setup();
//     _faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
//     FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
//     faceData.isBFace = true;
//     faceData.trs = bcTRSs[iBCTRS];


     // We get the number of faces in this TRS
  nbGeoEnts = m_thisTRS->getLocalNbGeoEnts();  //getNbTRs();
  CFLog(VERBOSE, "BCPeriodic::setup 5 nbGeoEnts - "<<nbGeoEnts<<"\n");
     // We get the number of equations
  NbEqs = PhysicalModelStack::getActive()->getNbEq();

  /* Setup MPI-related parameters */
  setupMPI();

  /* Check for every local face if it is an east or a west face */
  std::vector<FlxPntStruct> westFaces;
  westFaces.reserve(nbGeoEnts); // oversized to avoid frequent memory reallocation

  std::vector<FlxPntStruct> eastFaces;
  eastFaces.reserve(nbGeoEnts); // oversized to avoid frequent memory reallocation

  CFLog(VERBOSE, "BCPeriodic::setup 5.1\n");
  //_localWestEastMap.reserve(nbGeoEnts);

  RealVector faceNormal(dim);
  RealVector centreNode(dim);
  FlxPntStruct thisFlxPnt;

// This has not been yet been initialized!
//void StdSetup::createFaceOrientationStartIndexes()
  
createFaceOrientationStartIndexes();
    map< std::string , vector< vector< CFuint > > >&
      bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
    vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[m_thisTRS->getName()];

  //CFLog(VERBOSE, "BCPeriodic::setup 5.2 - "<<bndFacesStartIdxsPerTRS<<"\n");
  cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;//4;// 
  CFLog(VERBOSE, "BCPeriodic::setup 5.3\n");
  const CFuint nbTRs = m_thisTRS->getNbTRs();
  CFLog(VERBOSE, "BCPeriodic::setup 5.4\n");
  // loop over TRs
  CFuint localFaceID = 0;
  bool isFirstOrientation = true;
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    CFLog(VERBOSE, "BCPeriodic::setup 6 - nbTRs "<< nbTRs << "\n");
    for (m_orient = 0; m_orient < nbOrients; ++m_orient){
       // start and stop index of the faces with this orientation
       const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
       const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

       CFLog(VERBOSE, "startFaceIdx "<<startFaceIdx<<"/ stopFaceIdx "<<stopFaceIdx<<" with orientation "<<m_orient<<"\n");
       // loop over faces with this orientation
       for (CFuint iFace = startFaceIdx; iFace < stopFaceIdx; ++iFace){

  //   for (CFuint iFace = 0; iFace<nbTRs; iFace++){
    
          //CFLog(VERBOSE, "BCPeriodic::setup 6.0 - "<<iFace<<"\n");
          /* Face setup */
          faceData.idx = iFace;
          GeometricEntity *const face = _faceBuilder->buildGE();
          const CFuint faceGlobalID = face->getID();
          _globalToLocalTRSFaceID.insert(faceGlobalID,localFaceID);
    
    
          /* Calculate coordinate */
/*  FV
       const std::vector<Node*>& nodes = *face->getNodes();
       const CFuint nbNodesInFace = nodes.size();
       centreNode = 0.;
       for(CFuint iNode=0; iNode< nbNodesInFace; ++iNode) {
         centreNode += *(nodes[iNode]);
       }
       centreNode /= nbNodesInFace;
*/
       //IA;  From ConvBndCorrectionsRHSFluxReconstruction::setBndFaceData(CFuint faceID)
       
    
       //IA: In FV, the matching is done comparing faces. In this case, we will work directly
       // with the flux points
          
          for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      // Compute coordinates
             m_flxPntCoords[iFlx] = face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
    }
         // compute face Jacobian vectors
         m_faceJacobVecs = face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
         
         
         //CFLog(VERBOSE,"here I am\n");
         // get face Jacobian vector sizes in the flux points
         DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
         //CFLog(VERBOSE," size of faceJacobVecSizeFaceFlxPnts = "<<faceJacobVecSizeFaceFlxPnts.size()<<"\n");
          for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
          {
           //  CFLog(VERBOSE,"here I am 2\n");

             /* Load this flux point in a structure */
             thisFlxPnt.setCentreCoordinates(m_flxPntCoords[iFlx]);
             thisFlxPnt.setGlobalFaceID(faceGlobalID);
             thisFlxPnt.setLocalFaceID(localFaceID);
             thisFlxPnt.setLocalFluxID(iFlx);
             thisFlxPnt.setOrient(m_orient);
             //CFLog(VERBOSE,"here I am 3\n");
             //IA: How to get access to normals in setup of the BC with FR? Can I also use sockets?
             //  In that case, the normals are stored for each face or flux point?
             //
             // Next day:
             //
             // In FR socket_normals is not initialized!!!! So we can do the match depending upon each
             // orientation. We will assume that inside each TRS we will have two boundary walls, each one with
             // diferent orientation. I think that it is safe to assume this.
             //
             // In this way, we only need to detect which two orientation are active and then match them! For sake
             // of clarity, the same notation as in FV (west and east) will be keep, but this can be changed (and optimized!!)
             /*
             DataHandle<CFreal> normals = socket_normals.getDataHandle();
             RealVector faceNormal(dim);
             for(CFuint iDim=0; iDim<dim; ++iDim) {
                faceNormal[iDim]=normals[startID+iDim];
             }
             CFLog(VERBOSE, "normal "<<faceNormal<<" --- Vector "<<translationVector<<"\n");
             CFreal proj = MathTools::MathFunctions::innerProd(translationVector, faceNormal);
             if (proj > 0) {
                eastFaces.push_back(thisFlxPnt);
                //_localWestEastMap.insert(iFace,1);
                //std::cout<<"EAST \t rank = "<<_rank<<" faceGlobalID["<<iFace<<"] = "<< faceGlobalID<<"\n";
             }
             else if (proj < 0){
                westFaces.push_back(thisFlxPnt);
                //_localWestEastMap.insert(iFace,0);
             } 
             */
             if(_nbProcesses==1){
             if (!isFirstOrientation){
                eastFaces.push_back(thisFlxPnt);
             }else{
                westFaces.push_back(thisFlxPnt);
             }
             }
             else{
             
                 // get face Jacobian vector size
               //  CFLog(VERBOSE,"here I am 2\n");
                 //CFLog(VERBOSE,"size of m_faceJacobVecAbsSizeFlxPnts = "<< m_faceJacobVecAbsSizeFlxPnts.size()<<" size of faceJacobVecSizeFaceFlxPnts = "<<faceJacobVecSizeFaceFlxPnts.size()<<" faceGlobalID = "<<faceGlobalID<<" iFlx = "<<iFlx<<"\n");
    m_faceJacobVecAbsSizeFlxPnts[iFlx] = faceJacobVecSizeFaceFlxPnts[faceGlobalID][iFlx];
    //CFLog(VERBOSE,"here I am 4\n");
    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPnts[iFlx] = m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[m_orient];
    //CFLog(VERBOSE,"here I am 5\n");
    // set unit normal vector
    m_unitNormalFlxPnts[iFlx] = m_faceJacobVecs[iFlx]/m_faceJacobVecAbsSizeFlxPnts[iFlx];
      //       CFLog(VERBOSE,"here I am 6\n");
    CFreal proj = MathTools::MathFunctions::innerProd(translationVector, m_unitNormalFlxPnts[iFlx]);
    //CFLog(VERBOSE,"here I am 7\n");
             if (proj>0){
                eastFaces.push_back(thisFlxPnt);
             }else if(proj < 0){
                westFaces.push_back(thisFlxPnt);
             }
             }
             /* release geometry */
             _faceBuilder->releaseGE();
          } // End of the flux points loop
          localFaceID++;
       } // End of the faces loop
       if (startFaceIdx!=stopFaceIdx){
          isFirstOrientation=false;
       }
    } // End of the orientation loop
  } // End of the TR loop
  //_localWestEastMap.sortKeys();
      
  CFLog(VERBOSE, "BCPeriodic::setup 6.2\n");

    //IA: Now we should have classified the flux points of the TRS in west and east types
    // Now we have to match them using the user-defined translation vector.
    //
    // We aim to have two arrays. The first of those will contain the globalFaceID of the face
    // to which the flux point is linked to. The second one, will contain the local flux point
    // index at that face. They are, respectively:
    //   * _faceConnectivityMap
    //   * _fluxPointConnectivityMap
    //
    // But as we are using now the orientation, we do not know if we are treating actual west or
    // east faces. This can cause that we try to look for the translated flux points outside of the
    // actual domain. An easy fix to this is to do the match between west and east faces. If there is
    // no match, interchange them. This is not the least computationally intensive solution, but as it is
    // only done in the setup will not have huge impact on the overall runtime. 

    


 CFLog(VERBOSE,"here6.2.1");
  CFuint nbWestFaces = westFaces.size();
  CFuint nbEastFaces = eastFaces.size();
  CFLog(VERBOSE,"here6.2.1");
  std::vector<CFuint> nbWestFacesPerProcess(_nbProcesses,0);
  std::vector<CFuint> nbEastFacesPerProcess(_nbProcesses,0);
  CFLog(VERBOSE,"here6.2.2");
  MPI_Allgather(&nbWestFaces, 1, MPIStructDef::getMPIType(&nbWestFaces), 
	&nbWestFacesPerProcess[0], 1, MPIStructDef::getMPIType(&nbWestFaces), _comm);
  MPI_Allgather(&nbEastFaces, 1, MPIStructDef::getMPIType(&nbEastFaces), 
	&nbEastFacesPerProcess[0], 1, MPIStructDef::getMPIType(&nbEastFaces), _comm);
  
  CFuint totalNbWestFaces = 0;
  CFuint totalNbEastFaces = 0;
  CFLog(VERBOSE,"here6.2.3");
  if (_nbProcesses > 1) {
    CFLog(VERBOSE, " ==> number of west faces on each processor: \n");
    for(CFuint iP=0; iP<_nbProcesses; ++iP) {
      CFLog(VERBOSE, "     "<<iP<<" --> " << nbWestFacesPerProcess[iP] << " \n");
      totalNbWestFaces += nbWestFacesPerProcess[iP];
    }
    CFLog(VERBOSE, "     sum = " << totalNbWestFaces << " \n");
    CFLog(VERBOSE, " ==> number of east faces on each processor: \n");
    for(CFuint iP=0; iP<_nbProcesses; ++iP) {
      CFLog(VERBOSE, "     "<<iP<<" --> " << nbEastFacesPerProcess[iP] << " \n");
      totalNbEastFaces += nbEastFacesPerProcess[iP];
    }
    CFLog(VERBOSE, "     sum = " << totalNbEastFaces << " \n");
  }
  else{
    CFLog(VERBOSE, " ==> number of west faces: " << nbWestFaces << " \n");
    CFLog(VERBOSE, " ==> number of east faces: " << nbEastFaces << " \n");
  }
/* Definition of MPI struct */
  FaceMPIStruct faceMPIStruct;
 // CFLog(VERBOSE, "nbWestFaces " << nbWestFaces << " nbEastFaces " << nbEastFaces << "\n");
  _faceConnectivityMap.resize(nbGeoEnts*m_nbrFaceFlxPnts);
  _fluxPointConnectivityMap.resize(nbGeoEnts*m_nbrFaceFlxPnts);
  _orientMap.resize(nbGeoEnts*m_nbrFaceFlxPnts);
  
   /* 
   * Assemble local connectivity map using MPI. It stores for each local face
   0* on which processor the periodic face is, and its local and global ID
   */
  _localConnectivityMap.resize(nbGeoEnts*m_nbrFaceFlxPnts);

  CFuint matches = 0;
  if(_nbProcesses == 1){  
  /*  for(CFuint iP=0; iP<_nbProcesses; ++iP) { */

  /* Find WestFace connectivity */
  const CFuint nbW = nbWestFaces; //nbWestFacesPerProcess[iP];
  for(CFuint iFlx=0; iFlx<nbW; iFlx++) {

     /* fill the mpi_struct */
     //if(iP == _rank) faceMPIStruct.copy(westFaces[iFace]);

     /* MPI Broadcast struct */
     //faceMPIStruct.broadcast(iP);

     /* Unload the buffer */
     FlxPntStruct& westFace = westFaces[iFlx]; //= static_cast<FlxPntStruct>(faceMPIStruct);

     /* Find eastFace corresponding to this westFace */
     std::vector<FlxPntStruct>::iterator found = std::find_if(eastFaces.begin(), eastFaces.end(), findTranslated(westFace,translationVector,_threshold));
     if (found != eastFaces.end()) {
        //if (_localConnectivityMap[found->getLocalFaceID()].size() == 0)   matches++;
        //_localConnectivityMap[found->getLocalFaceID()].push_back(PairStruct(iP,westFace.getLocalFaceID()));
        CFLog(VERBOSE, "a - westFace "<<westFace.getGlobalFaceID()<<" has been matched with east "<<found->getGlobalFaceID()<<"\n");
        _faceConnectivityMap[westFace.getLocalFaceID()*m_nbrFaceFlxPnts + westFace.getLocalFluxID()] = found->getLocalFaceID();
        _fluxPointConnectivityMap[westFace.getLocalFaceID()*m_nbrFaceFlxPnts + westFace.getLocalFluxID()] = found->getLocalFluxID();
        _orientMap[westFace.getLocalFaceID()*m_nbrFaceFlxPnts + westFace.getLocalFluxID()] = found->getOrient();
        matches++;
     }else{
        std::vector<FlxPntStruct>::iterator found = std::find_if(eastFaces.begin(), eastFaces.end(), findTranslated(westFace,backTranslationVector,_threshold));
        if (found != eastFaces.end()) {
           CFLog(VERBOSE, "b - westFace "<<westFace.getGlobalFaceID()<<" has been matched with east "<<found->getGlobalFaceID()<<"\n");
           _faceConnectivityMap[westFace.getLocalFaceID()*m_nbrFaceFlxPnts + westFace.getLocalFluxID()] = found->getLocalFaceID();
           _fluxPointConnectivityMap[westFace.getLocalFaceID()*m_nbrFaceFlxPnts + westFace.getLocalFluxID()] = found->getLocalFluxID();
           _orientMap[westFace.getLocalFaceID()*m_nbrFaceFlxPnts + westFace.getLocalFluxID()] = found->getOrient();
           matches++;
        } 
     }

  } // end westFace connectivity

  const CFuint nbE = nbEastFaces; //nbEastFacesPerProcess[iP];
  for(CFuint iFlx=0; iFlx<nbE; iFlx++) {

     /* fill the mpi_struct */
     //if(iP == _rank) faceMPIStruct.copy(eastFaces[iFace]);

     /* MPI Broadcast struct */
     //faceMPIStruct.broadcast(iP);

     /* Unload the buffer */
     FlxPntStruct& eastFace = eastFaces[iFlx]; //static_cast<FlxPntStruct>(faceMPIStruct);

     /* Find westFace corresponding to this eastFace */
     std::vector<FlxPntStruct>::iterator found = std::find_if(westFaces.begin(), westFaces.end(), findTranslated(eastFace,backTranslationVector,_threshold));
     if (found != westFaces.end()) {
        //if (_localConnectivityMap[found->getLocalFaceID()].size() == 0)   matches++;
        //_localConnectivityMap[found->getLocalFaceID()].push_back(PairStruct(iP,eastFace.getLocalFaceID()));
        CFLog(VERBOSE, "a - eastFace "<<eastFace.getGlobalFaceID()<<" has been matched with west "<<found->getGlobalFaceID()<<"\n");
        _faceConnectivityMap[eastFace.getLocalFaceID()*m_nbrFaceFlxPnts + eastFace.getLocalFluxID()] = found->getLocalFaceID();
        _fluxPointConnectivityMap[eastFace.getLocalFaceID()*m_nbrFaceFlxPnts + eastFace.getLocalFluxID()] = found->getLocalFluxID();
        _orientMap[eastFace.getLocalFaceID()*m_nbrFaceFlxPnts + eastFace.getLocalFluxID()] = found->getOrient();
        matches++;
     }else{
        std::vector<FlxPntStruct>::iterator found = std::find_if(westFaces.begin(), westFaces.end(), findTranslated(eastFace,translationVector,_threshold));
        if (found != westFaces.end()) {
           //if (_localConnectivityMap[found->getLocalFaceID()].size() == 0)   matches++;
           //_localConnectivityMap[found->getLocalFaceID()].push_back(PairStruct(iP,eastFace.getLocalFaceID()));
           CFLog(VERBOSE, "b - eastFace "<<eastFace.getGlobalFaceID()<<" has been matched with west "<<found->getGlobalFaceID()<<"\n");
           _faceConnectivityMap[eastFace.getLocalFaceID()*m_nbrFaceFlxPnts + eastFace.getLocalFluxID()] = found->getLocalFaceID();
           _fluxPointConnectivityMap[eastFace.getLocalFaceID()*m_nbrFaceFlxPnts + eastFace.getLocalFluxID()] = found->getLocalFluxID();
           _orientMap[eastFace.getLocalFaceID()*m_nbrFaceFlxPnts + eastFace.getLocalFluxID()] = found->getOrient();
           matches++;
        }
     }
  } // end eastFace connectivity

  //} // End of nbProcesses loop

  } // End if nbProcesses=1
  
    else { 
    std::vector<CFuint> sendGlobalIDW(nbWestFaces, 0);
    std::vector<CFuint> sendLocalIDW(nbWestFaces, 0);
    std::vector<CFuint> sendLocalFluxIDW(nbWestFaces, 0); 
    std::vector<CFreal> sendCoordsW(nbWestFaces*dim, 0);
    std::vector<CFuint> recvGlobalIDW(totalNbWestFaces, 0);
    std::vector<CFuint> recvLocalIDW(totalNbWestFaces, 0);
    std::vector<CFuint> recvLocalFluxIDW(totalNbWestFaces, 0);
    std::vector<CFreal> recvCoordsW(totalNbWestFaces*dim, 0);

    std::vector<int> displsW(_nbProcesses, 0);
    std::vector<int> rcountsW(_nbProcesses, 0);
    std::vector<int> displsWDim(_nbProcesses, 0);
    std::vector<int> rcountsWDim(_nbProcesses, 0);

    int localLengthW = nbWestFaces;
    int localLengthWDim = nbWestFaces*dim;

    std::vector<CFuint> sendGlobalIDE(nbEastFaces, 0);
    std::vector<CFuint> sendLocalIDE(nbEastFaces, 0);
    std::vector<CFuint> sendLocalFluxIDE(nbEastFaces, 0);
    std::vector<CFreal> sendCoordsE(nbEastFaces*dim, 0);
    std::vector<CFuint> recvGlobalIDE(totalNbEastFaces, 0);
    std::vector<CFuint> recvLocalIDE(totalNbEastFaces, 0);
    std::vector<CFuint> recvLocalFluxIDE(totalNbEastFaces, 0);
    std::vector<CFreal> recvCoordsE(totalNbEastFaces*dim, 0);

    std::vector<int> displsE(_nbProcesses, 0);
    std::vector<int> rcountsE(_nbProcesses, 0);
    std::vector<int> displsEDim(_nbProcesses, 0);
    std::vector<int> rcountsEDim(_nbProcesses, 0);

    int localLengthE = nbEastFaces;
    int localLengthEDim = nbEastFaces*dim;
    
    displsW[0] = 0;
    displsE[0] = 0;
    for (CFuint iP=1;iP<_nbProcesses;iP++) {
      displsW[iP] = displsW[iP-1] + nbWestFacesPerProcess[iP-1];
      displsWDim[iP] = displsW[iP]*dim;
      displsE[iP] = displsE[iP-1] + nbEastFacesPerProcess[iP-1];
      displsEDim[iP] = displsE[iP]*dim;
    }

    for(CFuint iP=0; iP<_nbProcesses; ++iP) {
      rcountsW[iP] = nbWestFacesPerProcess[iP];
      rcountsWDim[iP] = rcountsW[iP]*dim;
      rcountsE[iP] = nbEastFacesPerProcess[iP];
      rcountsEDim[iP] = rcountsE[iP]*dim;
    }
      
    /* Fill the send vector */
    for(CFuint iFace=0; iFace<nbWestFaces; iFace++) {

      /* fill the sentd info */
      sendGlobalIDW[iFace] = westFaces[iFace].getGlobalFaceID();
      sendLocalIDW[iFace] = westFaces[iFace].getLocalFaceID();//+westFaces[iFace].getLocalFluxID();  //Changed from FVMPI
      sendLocalFluxIDW[iFace] = westFaces[iFace].getLocalFluxID();
      for (CFuint iDim = 0; iDim<dim; iDim++){
	sendCoordsW[iFace*dim + iDim] = westFaces[iFace].getCentreCoordinates()[iDim];
      }
    }
    for(CFuint iFace=0; iFace<nbEastFaces; iFace++) {

      /* fill the sentd info */
      sendGlobalIDE[iFace] = eastFaces[iFace].getGlobalFaceID();
      sendLocalIDE[iFace] = eastFaces[iFace].getLocalFaceID();//+eastFaces[iFace].getLocalFluxID(); //Changed from FVMPI
      sendLocalFluxIDE[iFace] = eastFaces[iFace].getLocalFluxID();
      for (CFuint iDim = 0; iDim<dim; iDim++){
	sendCoordsE[iFace*dim + iDim] = eastFaces[iFace].getCentreCoordinates()[iDim];
      }
    }
    MPI_Barrier(_comm);

    MPI_Datatype MPI_CFUINT = Common::MPIStructDef::getMPIType(&sendGlobalIDW[0]);
    MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&sendCoordsW[0]);
    MPI_Allgatherv(&sendGlobalIDW[0], localLengthW, MPI_CFUINT,
		    &recvGlobalIDW[0], &rcountsW[0], &displsW[0], MPI_CFUINT, _comm);
    MPI_Allgatherv(&sendLocalIDW[0], localLengthW, MPI_CFUINT,
		    &recvLocalIDW[0], &rcountsW[0], &displsW[0], MPI_CFUINT, _comm);
    MPI_Allgatherv(&sendLocalFluxIDW[0], localLengthW, MPI_CFUINT,
		    &recvLocalFluxIDW[0], &rcountsW[0], &displsW[0], MPI_CFUINT, _comm);
    MPI_Allgatherv(&sendCoordsW[0], localLengthWDim, MPI_CFREAL,
		    &recvCoordsW[0], &rcountsWDim[0], &displsWDim[0], MPI_CFREAL, _comm);
    MPI_Allgatherv(&sendGlobalIDE[0], localLengthE, MPI_CFUINT,
		    &recvGlobalIDE[0], &rcountsE[0], &displsE[0], MPI_CFUINT, _comm);
    MPI_Allgatherv(&sendLocalIDE[0], localLengthE, MPI_CFUINT,
		    &recvLocalIDE[0], &rcountsE[0], &displsE[0], MPI_CFUINT, _comm);
    MPI_Allgatherv(&sendLocalFluxIDE[0], localLengthE, MPI_CFUINT,
		    &recvLocalFluxIDE[0], &rcountsE[0], &displsE[0], MPI_CFUINT, _comm);
    MPI_Allgatherv(&sendCoordsE[0], localLengthEDim, MPI_CFREAL,
		    &recvCoordsE[0], &rcountsEDim[0], &displsEDim[0], MPI_CFREAL, _comm);
    MPI_Barrier(_comm);

    // Constructing all the faces
    /// WEST
    std::vector<FlxPntStruct> totalWestFaces;
    RealVector coordinates(dim);
    totalWestFaces.resize(totalNbWestFaces);
    for (CFuint iFace=0; iFace<totalNbWestFaces; iFace++) {
      for (CFuint iDim=0; iDim<dim; iDim++) {
	coordinates[iDim] = recvCoordsW[iFace*dim + iDim];
      }
      /* Load this face in a structure */
      totalWestFaces[iFace].setCentreCoordinates(coordinates);
      totalWestFaces[iFace].setGlobalFaceID(recvGlobalIDW[iFace]);
      totalWestFaces[iFace].setLocalFaceID(recvLocalIDW[iFace]);
      totalWestFaces[iFace].setLocalFluxID(recvLocalFluxIDW[iFace]);
    }
    ///EAST
    std::vector<FlxPntStruct> totalEastFaces(totalNbEastFaces);
    for (CFuint iFace=0; iFace<totalNbEastFaces; iFace++) {
      for (CFuint iDim=0; iDim<dim; iDim++) {
	coordinates[iDim] = recvCoordsE[iFace*dim + iDim]; 
      }
      /* Load this face in a structure */
      totalEastFaces[iFace].setCentreCoordinates(coordinates);
      totalEastFaces[iFace].setGlobalFaceID(recvGlobalIDE[iFace]);
      totalEastFaces[iFace].setLocalFaceID(recvLocalIDE[iFace]);
      totalEastFaces[iFace].setLocalFluxID(recvLocalFluxIDE[iFace]);
    }
   // FlxPntStruct& westFace = westFaces[iFace];
    //FlxPntStruct& eastFace = eastFaces[iFace];
    /* Find WestFace connectivity */
    for(CFuint iP=0; iP<_nbProcesses; iP++){
      const CFuint nbW = nbWestFacesPerProcess[iP]; 
      for(CFuint iFace=0; iFace<nbW; iFace++) {
	FlxPntStruct westFace = totalWestFaces[displsW[iP] + iFace];
	/* Find eastFace corresponding to this westFace */
	std::vector<FlxPntStruct>::iterator found = std::find_if(eastFaces.begin(), eastFaces.end(), findTranslated(westFace,translationVector,_threshold));
	if (found != eastFaces.end()) {
          CFLog(VERBOSE, "a - westFace "<<westFace.getGlobalFaceID()<<" has been matched with east "<<found->getGlobalFaceID()<<"\n");
	  if (_localConnectivityMap[found->getLocalFaceID()*m_nbrFaceFlxPnts+found->getLocalFluxID()].size() == 0)   matches++;
	  _localConnectivityMap[found->getLocalFaceID()*m_nbrFaceFlxPnts+found->getLocalFluxID()].push_back(PairStruct(iP,westFace.getLocalFaceID()*m_nbrFaceFlxPnts+westFace.getLocalFluxID()));
          CFLog(VERBOSE,"found local face id = "<<found->getLocalFaceID()<<" westFace local Face id = "<<westFace.getLocalFaceID()<<"\n");
          //CFLog(VERBOSE,"coordintes of found = ")
	}else{
             std::vector<FlxPntStruct>::iterator found = std::find_if(eastFaces.begin(), eastFaces.end(), findTranslated(westFace,backTranslationVector,_threshold));
             if (found != eastFaces.end()) {
          CFLog(VERBOSE, "b - westFace "<<westFace.getGlobalFaceID()<<" has been matched with east "<<found->getGlobalFaceID()<<"\n");
	  if (_localConnectivityMap[found->getLocalFaceID()*m_nbrFaceFlxPnts+found->getLocalFluxID()].size() == 0)   matches++;
	  _localConnectivityMap[found->getLocalFaceID()*m_nbrFaceFlxPnts+found->getLocalFluxID()].push_back(PairStruct(iP,westFace.getLocalFaceID()*m_nbrFaceFlxPnts+westFace.getLocalFluxID()));
	}
        }
      } 
      
      /* Find eastFace connectivity */
      const CFuint nbE = nbEastFacesPerProcess[iP];
      for(CFuint iFace=0; iFace<nbE; iFace++) {
	//CFLog(VERBOSE,"iD = "<<displsE[iP] + iFace <<"\n");
	FlxPntStruct eastFace = totalEastFaces[displsE[iP] + iFace]; 
	
	/* Find westFace corresponding to this eastFace */
	std::vector<FlxPntStruct>::iterator found = std::find_if(westFaces.begin(), westFaces.end(), findTranslated(eastFace,backTranslationVector,_threshold));
	if (found != westFaces.end()) {
          CFLog(VERBOSE, "a - eastFace "<<eastFace.getGlobalFaceID()<<" has been matched with west "<<found->getGlobalFaceID()<<"\n");
	  if (_localConnectivityMap[found->getLocalFaceID()*m_nbrFaceFlxPnts+found->getLocalFluxID()].size() == 0)   matches++;
	  _localConnectivityMap[found->getLocalFaceID()*m_nbrFaceFlxPnts+found->getLocalFluxID()].push_back(PairStruct(iP,eastFace.getLocalFaceID()*m_nbrFaceFlxPnts+eastFace.getLocalFluxID()));
          CFLog(VERBOSE,"found local face id = "<<found->getLocalFaceID()<<" eastFace local Face id = "<<eastFace.getLocalFaceID()<<"\n");
	}else{
            std::vector<FlxPntStruct>::iterator found = std::find_if(westFaces.begin(), westFaces.end(), findTranslated(eastFace,translationVector,_threshold));
	if (found != westFaces.end()) {
          CFLog(VERBOSE, "b - eastFace "<<eastFace.getGlobalFaceID()<<" has been matched with west "<<found->getGlobalFaceID()<<"\n");
	  if (_localConnectivityMap[found->getLocalFaceID()*m_nbrFaceFlxPnts+found->getLocalFluxID()].size() == 0)   matches++;
	  _localConnectivityMap[found->getLocalFaceID()*m_nbrFaceFlxPnts+found->getLocalFluxID()].push_back(PairStruct(iP,eastFace.getLocalFaceID()*m_nbrFaceFlxPnts+eastFace.getLocalFluxID()));
	}
        }
      } // end eastFace connectivity
    }
    
    MPI_Barrier(_comm);    
  }
  
 // CFLog(VERBOSE,"BCPeriodic::setup() => step 3 took " << stp.read() << "s\n");
  
  if(matches!=m_nbrFaceFlxPnts*nbGeoEnts)
    CFLog(VERBOSE, "Only " << matches << "/"<< m_nbrFaceFlxPnts*nbGeoEnts << " matches found. Wrong TranslationVector or try increasing Threshold \n");
  
  CFLog(VERBOSE, "BCPeriodic::setup() => end\n");


   CFLog(VERBOSE, "matches "<<matches<<"\n");
   CFLog(VERBOSE, "END BCPeriodic::setup \n");

}

//////////////////////////////////////////////////////////////////////////////
//Added from FV
void BCPeriodic::setupMPI() 
{
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  _nbProcesses = PE::GetPE().GetProcessorCount(nsp); //gets the no. of processors being used
  _rank = PE::GetPE().GetRank(nsp); //returns the number id of the processor
  _comm = PE::GetPE().GetCommunicator(nsp);//sends data from one core to the the other
  _nbFacesPerProcess.resize(_nbProcesses,0);
  _periodicState = new Framework::State;
  _periodicState->resize(NbEqs);
//_periodicState->resize(_nE);
  CFLog(VERBOSE,"_nbFacesPerProcess = "<<_nbFacesPerProcess[1]<<"\n");
  MPI_Barrier(_comm);//the processing of a an individual processor is paused temporarily till all the processors compute till some certain values
  CFLog(VERBOSE,"here\n");
  MPI_Allgather(&nbGeoEnts, 1, MPIStructDef::getMPIType(&nbGeoEnts), 
	&_nbFacesPerProcess[0], 1, MPIStructDef::getMPIType(&nbGeoEnts), _comm);//gathers data from all processors
CFLog(VERBOSE,"_nbFacesPerProcess = "<<_nbFacesPerProcess[1]<<"\n");
  //MPI_Allgather(&_nbTrsFaces, 1, MPIStructDef::getMPIType(&_nbTrsFaces), 
	//&_nbFacesPerProcess[0], 1, MPIStructDef::getMPIType(&_nbTrsFaces), _comm);//gathers data from all processors

  
  // build counts and displacements for preProcess step
  
  _sendcounts.resize(_nbProcesses);
  _recvcounts.resize(_nbProcesses);
  _recvdispls.resize(_nbProcesses);
  _senddispls.resize(_nbProcesses);
  
  _recvdispls[0] = 0;
  _senddispls[0] = 0;
  _LastDisplacement = 0;
  for(CFuint iP=0; iP<_nbProcesses; iP++){
    _sendcounts[iP] = m_nbrFaceFlxPnts*NbEqs*_nbFacesPerProcess[_rank];
    _recvcounts[iP] = m_nbrFaceFlxPnts*NbEqs*_nbFacesPerProcess[iP];
//_sendcounts[iP] = _nE*_nbFacesPerProcess[_rank];
//    _recvcounts[iP] = _nE*_nbFacesPerProcess[iP];
  }
  if(_nbProcesses>1){
    for(CFuint iP=1; iP<_nbProcesses; iP++){
      _recvdispls[iP] = _recvdispls[iP-1] + _recvcounts[iP-1];
      _senddispls[iP] = _senddispls[iP-1] + _sendcounts[iP-1];
    }
  }
  _LastDisplacement = _recvdispls[_nbProcesses-1] + _recvcounts[_nbProcesses-1];
  
  _recvbuf.resize(_LastDisplacement, 0.);
  _sendbuf.reserve(m_nbrFaceFlxPnts*nbGeoEnts*NbEqs*_nbProcesses);
//  _sendbuf.reserve(_nbTrsFaces*_nE*_nbProcesses);
  
  _recvbufLimiter.resize(_LastDisplacement, 0.);
  _sendbufLimiter.reserve(m_nbrFaceFlxPnts*nbGeoEnts*NbEqs*_nbProcesses);
//_sendbufLimiter.reserve(_nbTrsFaces*_nE*_nbProcesses);  

  // allocate data for the transfer of cell-based gradients
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _recvbufGrad.resize(dim);
  for (CFuint i = 0; i < dim; ++i) {
    _recvbufGrad[i].resize(_LastDisplacement, 0.);
  }
  _sendbufGrad.resize(dim);
  for (CFuint i = 0; i < dim; ++i) {
    _sendbufGrad[i].reserve(m_nbrFaceFlxPnts*nbGeoEnts*NbEqs*_nbProcesses);
//_sendbufGrad[i].reserve(_nbTrsFaces*_nE*_nbProcesses);
  }
}
//till here

//////////////////////////////////////////////////////////////////////////////
void BCPeriodic::preProcess()
{

  CFLog(VERBOSE, "BCPeriodic::preProcess() => start\n");
  
  if (_nbProcesses > 1) {
      
      CFAUTOTRACE;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current bnd face TRS
  SafePtr<TopologicalRegionSet> faceTrs = m_thisTRS;
      
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
    
    _sendbuf.clear();

  std::vector<CFreal> boundaryState;
    if (nbGeoEnts > 0) {
      boundaryState.reserve(m_nbrFaceFlxPnts*nbGeoEnts*NbEqs); 
    }
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
	 // setBndFaceData(m_face->getID());//faceID 
	  
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
      if (m_intCell->getID() == 223) CFLog(VERBOSE, "sol: " << *((*m_cellStates)[solIdx])  << "\n");
    }
    
    for(CFuint e=0; e<NbEqs; e++){
        boundaryState.push_back((*(m_cellStatesFlxPnt[iFlxPnt]))[e]);
      }
    
    CFLog(VERBOSE, "inner: " << *(m_cellStatesFlxPnt[iFlxPnt])  << "\n");//if (m_intCell->getID() == 223) 
  }
  
    

  } //if cell is parallel updatable
        
      m_faceBuilder->releaseGE();
      } //end of loop over faces with this orientation
    } // end of orientation loop
  } //end of TR loop
  
  const CFuint nbFluxPoints = m_nbrFaceFlxPnts*nbGeoEnts*NbEqs;//*_nbProcesses
    for(CFuint iP=0; iP<_nbProcesses; iP++){
      for(CFuint i=0; i< nbFluxPoints; i++){
        _sendbuf.push_back(boundaryState[i]);
        CFLog(VERBOSE,"_sendbuf = "<<_sendbuf[i]<<"\n");
      }
    }

    MPI_Barrier(_comm);    
    
    MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&_sendbuf[0]);
    
    MPI_Alltoallv(&_sendbuf[0], &_sendcounts[0], &_senddispls[0], MPI_CFREAL,
		  &_recvbuf[0], &_recvcounts[0], &_recvdispls[0], MPI_CFREAL, _comm);
    MPI_Barrier(_comm);
  } //if nb processes > 1
  CFLog(VERBOSE, "BCPeriodic::preProcess() => end\n");
}

//////////////////////////////////////////////////////////////////////////////
  


//////////////////////////////////////////////////////////////////////////////
void BCPeriodic::unsetup()
{
  CFAUTOTRACE;
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    deletePtr(m_cellStatesFlxPnt[iFlx]);
  }
  m_cellStatesFlxPnt.clear();

  // unsetup of the parent class
  BCStateComputer::unsetup();
}
    
//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

