#include "Framework/BadFormatException.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/BCPeriodic.hh"

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
  m_flxLocalCoords(), 
  m_flxPntCoords(),
  m_nbrFaceFlxPnts(),
  m_dim(),
  m_orient(),
  m_thisTRS(),
  m_nbrSolDep(),
  m_solPolyValsAtFlxPnts(),
  m_otherFace(),
  m_intCell(),
  m_cellStates(),
  m_cellGrads(),
  m_flxSolDep(),
  m_faceFlxPntConn(),
  _globalToLocalTRSFaceID(),
  _localWestEastMap(),
  _orientMap(),
  _faceConnectivityMap(),
  _fluxPointConnectivityMap(),
  _faceBuilder(),
  socket_gradients("gradients"),
//  socket_states("states"),
//  socket_gstates("gstates"),
//  socket_nodes("nodes"),
  _translationVector()
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
{CFLog(VERBOSE,"here2\n");
  CFLog(VERBOSE, "BCPeriodic - Computing ghost states.. \n");
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbrStates = intStates.size();
  cf_assert(nbrStates == ghostStates.size());
  cf_assert(nbrStates == normals.size());

  // we prepare the face builder
  FaceToCellGEBuilder::GeoData& faceData = _faceBuilder->getDataGE();
  faceData.cellsTRS = cellTrs;
  faceData.facesTRS = m_thisTRS; 
  faceData.isBoundary = true;
  CFLog(VERBOSE, "BCPeriodic - Preparing face builder.. \n");

  // Get the localFaceID from the map, knowing the faceGlobalID
  const CFuint faceLocalID = _globalToLocalTRSFaceID.find(m_face->getID());
CFLog(VERBOSE,"here3\n");
  // If we assume that all flux point of this face are connected to same face at the other side, we
  // can do this outside the flux points loop
  const CFuint otherFaceLocalID = _faceConnectivityMap[faceLocalID*nbrStates];
 CFLog(VERBOSE,"here3.1\n");
  CFLog(VERBOSE, "BCPeriodic - We are studying face "<<m_face->getID()<<" with faceLocalID "<<faceLocalID<<" linked to "<<otherFaceLocalID<<"\n");
  // If we take this assumption further, we could store the _faceConnectivityMap only for the faces,
  // instead that for each flux point!  
CFLog(VERBOSE,"here3.2\n");
  // We build the geometric entity of that face
  faceData.idx = otherFaceLocalID;
  CFLog(VERBOSE,"here3.3\n");
  m_otherFace = _faceBuilder->buildGE();
  //const CFuint otherFaceGlobalID = m_otherFace->getID();
CFLog(VERBOSE,"here4\n");
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
     CFLog(VERBOSE, "BCPeriodic - Computing BC at FlxPnt "<<iState<< ", linked to "<<iFlxPnt<<"\n");

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
	CFLog(VERBOSE,"Bx = "<<(*(ghostStates[iState]))[0]<<"\n");
  } 
CFLog(VERBOSE,"here5\n");
  // release the face
  _faceBuilder->releaseGE();
  CFLog(VERBOSE, "BCPeriodic - Computing BC at FlxPnt END \n");

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
     // if (i == 2 &&  (*ghostStates[iState])[i] < -0.000001) CFLog(INFO, "intState: " << *intStates[iState] << ", ghost: " << *ghostStates[iState] << "\n");
   }
   
   //CFLog(DEBUG_MAX, "Periodic::setGhostState() => ghostState = " << *ghostState << "\n"); 
  }
*/
}

//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::computeGhostGradients
(const std::vector< std::vector< RealVector* > >& intGrads,
 std::vector< std::vector< RealVector* > >& ghostGrads,
 const std::vector< RealVector >& normals,
 const std::vector< RealVector >& coords)
{

  CFLog(INFO, "BCPeriodic - Computing ghost gradients.. CAUTION! This has not been tested yet!\n");
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  const CFuint nbrGradVars = intGrads[0].size();
  const CFuint nbrStateGrads = intGrads.size();
  const CFuint nbrStates = m_cellStates->size();
  cf_assert(nbrStates == ghostGrads.size());
  cf_assert(nbrStates == normals.size());

  // we prepare the face builder
  FaceToCellGEBuilder::GeoData& faceData = _faceBuilder->getDataGE();
  faceData.cellsTRS = cellTrs;
  faceData.facesTRS = m_thisTRS; //store both references as protecteed values
  faceData.isBoundary = true;
  CFLog(VERBOSE, "BCPeriodic - (grads) Preparing face builder.. \n");

  // Get the localFaceID from the map, knowing the faceGlobalID
  const CFuint faceLocalID = _globalToLocalTRSFaceID.find(m_face->getID());

  // If we assume that all flux point of this face are connected to same face at the other side, we
  // can do this outside the flux points loop
  const CFuint otherFaceLocalID = _faceConnectivityMap[faceLocalID*nbrStateGrads];
 
  CFLog(VERBOSE, "BCPeriodic - (grads) We are studying face "<<m_face->getID()<<" with faceLocalID "<<faceLocalID<<" linked to "<<otherFaceLocalID<<"\n");
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
     CFLog(VERBOSE, "BCPeriodic - (grads) Computing BC at FlxPnt "<<iState<< ", linked to "<<iFlxPnt<<"\n");

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
  CFLog(VERBOSE, "BCPeriodic - (grads) Computing BC at FlxPnt END \n");



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
  
  return result;
}
             
//////////////////////////////////////////////////////////////////////////////


void BCPeriodic::createFaceOrientationStartIndexes()
{
  CFAUTOTRACE;

  CFLog(INFO, "\n\n createFaceOrientationStartIndexes() \n\n");

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
  CFLog(INFO, "BCPeriodic::setup 1 \n");


  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  _faceBuilder = getMethodData().getSecondFaceBuilder();


  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  CFLog(INFO, "BCPeriodic::setup 2 \n");

  // get the face local coords of the flux points on one face
  m_flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  m_nbrFaceFlxPnts = m_flxLocalCoords->size();
  m_dim   = PhysicalModelStack::getActive()->getDim();


  // get the face - flx pnt connectivity per orient (ConvBnd) 
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConn();
  m_flxSolDep = frLocalData[0]->getFlxPntSolDependency();
  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();

  m_nbrSolDep = ((*m_flxSolDep)[0]).size();

  m_flxPntCoords.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx].resize(m_dim);
  }
  CFLog(INFO, "BCPeriodic::setup 3 - "<<m_dim<<" "<<m_nbrFaceFlxPnts<<"\n");


  /* Check if translation vector has right dimension */
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_assert(_translationVector.size()==dim);

  /* Convert std::vector<CFreal> to RealVector version */
  RealVector translationVector(dim);
  RealVector backTranslationVector(dim);
  CFLog(INFO, "BCPeriodic::setup 4 dim "<<dim<<" \n");
  for(CFuint iDim=0; iDim<dim; ++iDim) {
    translationVector[iDim] = _translationVector[iDim];
    backTranslationVector[iDim] = -translationVector[iDim];
  }
  CFLog(INFO, " ==> translationVector = " << translationVector << " \n");

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
    CFLog(INFO, "iTRS "<< iTRS << " \n");
    if (m_trsNames[0]==trsList[iTRS]->getName() ){
       m_thisTRS = trsList[iTRS];
       CFLog(INFO, "Matching BC "<<m_trsNames[0]<<" with "<<m_thisTRS->getName() << "\n");
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
  CFLog(INFO,"End of TRS loop \n");

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
  const CFuint nbGeoEnts = m_thisTRS->getLocalNbGeoEnts();  //getNbTRs();
  CFLog(INFO, "BCPeriodic::setup 5 nbGeoEnts - "<<nbGeoEnts<<"\n");

  /* Check for every local face if it is an east or a west face */
  std::vector<FlxPntStruct> westFaces;
  westFaces.reserve(nbGeoEnts); // oversized to avoid frequent memory reallocation

  std::vector<FlxPntStruct> eastFaces;
  eastFaces.reserve(nbGeoEnts); // oversized to avoid frequent memory reallocation

  CFLog(INFO, "BCPeriodic::setup 5.1\n");
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

  //CFLog(INFO, "BCPeriodic::setup 5.2 - "<<bndFacesStartIdxsPerTRS<<"\n");
  cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = 4;//bndFacesStartIdxs[0].size()-1; 
  CFLog(INFO, "BCPeriodic::setup 5.3\n");
  const CFuint nbTRs = m_thisTRS->getNbTRs();
  CFLog(INFO, "BCPeriodic::setup 5.4\n");
  // loop over TRs
  CFuint localFaceID = 0;
  bool isFirstOrientation = true;
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    CFLog(INFO, "BCPeriodic::setup 6 - nbTRs "<< nbTRs << "\n");
    for (m_orient = 0; m_orient < nbOrients; ++m_orient){
       // start and stop index of the faces with this orientation
       const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
       const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

       CFLog(INFO, "startFaceIdx "<<startFaceIdx<<"/ stopFaceIdx "<<stopFaceIdx<<" with orientation "<<m_orient<<"\n");
       // loop over faces with this orientation
       for (CFuint iFace = startFaceIdx; iFace < stopFaceIdx; ++iFace){

  //   for (CFuint iFace = 0; iFace<nbTRs; iFace++){
    
          //CFLog(INFO, "BCPeriodic::setup 6.0 - "<<iFace<<"\n");
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

             /* Load this flux point in a structure */
             thisFlxPnt.setCentreCoordinates(m_flxPntCoords[iFlx]);
             thisFlxPnt.setGlobalFaceID(faceGlobalID);
             thisFlxPnt.setLocalFaceID(localFaceID);
             thisFlxPnt.setLocalFluxID(iFlx);
             thisFlxPnt.setOrient(m_orient);

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
             CFLog(INFO, "normal "<<faceNormal<<" --- Vector "<<translationVector<<"\n");
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

             if (!isFirstOrientation){
                eastFaces.push_back(thisFlxPnt);
             }else{
                westFaces.push_back(thisFlxPnt);
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
      
  CFLog(INFO, "BCPeriodic::setup 6.2\n");

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

    

  CFuint nbWestFaces = westFaces.size();
  CFuint nbEastFaces = eastFaces.size();

  CFLog(INFO, "nbWestFaces " << nbWestFaces << " nbEastFaces " << nbEastFaces << "\n");
  _faceConnectivityMap.resize(nbGeoEnts*m_nbrFaceFlxPnts);
  _fluxPointConnectivityMap.resize(nbGeoEnts*m_nbrFaceFlxPnts);
  _orientMap.resize(nbGeoEnts*m_nbrFaceFlxPnts);

  CFuint matches = 0;
  /* if(_nbProcesses == 1){  We assume =1 */
  /*  for(CFuint iP=0; iP<_nbProcesses; ++iP) { */

  /* Find WestFace connectivity */
  const CFuint nbW = nbWestFaces; //nbWestFacesPerProcess[iP];
  for(CFuint iFlx=0; iFlx<nbW; iFlx++) {

     /* fill the mpi_struct */
     //if(iP == _rank) faceMPIStruct.copy(westFaces[iFace]);

     /* MPI Broadcast struct */
     //faceMPIStruct.broadcast(iP);

     /* Unload the buffer */
     FlxPntStruct& westFace = westFaces[iFlx]; //= static_cast<FaceStruct>(faceMPIStruct);

     /* Find eastFace corresponding to this westFace */
     std::vector<FlxPntStruct>::iterator found = std::find_if(eastFaces.begin(), eastFaces.end(), findTranslated(westFace,translationVector,_threshold));
     if (found != eastFaces.end()) {
        //if (_localConnectivityMap[found->getLocalFaceID()].size() == 0)   matches++;
        //_localConnectivityMap[found->getLocalFaceID()].push_back(PairStruct(iP,westFace.getLocalFaceID()));
        CFLog(INFO, "a - westFace "<<westFace.getGlobalFaceID()<<" has been matched with east "<<found->getGlobalFaceID()<<"\n");
        _faceConnectivityMap[westFace.getLocalFaceID()*m_nbrFaceFlxPnts + westFace.getLocalFluxID()] = found->getLocalFaceID();
        _fluxPointConnectivityMap[westFace.getLocalFaceID()*m_nbrFaceFlxPnts + westFace.getLocalFluxID()] = found->getLocalFluxID();
        _orientMap[westFace.getLocalFaceID()*m_nbrFaceFlxPnts + westFace.getLocalFluxID()] = found->getOrient();
        matches++;
     }else{
        std::vector<FlxPntStruct>::iterator found = std::find_if(eastFaces.begin(), eastFaces.end(), findTranslated(westFace,backTranslationVector,_threshold));
        if (found != eastFaces.end()) {
           CFLog(INFO, "b - westFace "<<westFace.getGlobalFaceID()<<" has been matched with east "<<found->getGlobalFaceID()<<"\n");
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
     FlxPntStruct& eastFace = eastFaces[iFlx]; //static_cast<FaceStruct>(faceMPIStruct);

     /* Find westFace corresponding to this eastFace */
     std::vector<FlxPntStruct>::iterator found = std::find_if(westFaces.begin(), westFaces.end(), findTranslated(eastFace,backTranslationVector,_threshold));
     if (found != westFaces.end()) {
        //if (_localConnectivityMap[found->getLocalFaceID()].size() == 0)   matches++;
        //_localConnectivityMap[found->getLocalFaceID()].push_back(PairStruct(iP,eastFace.getLocalFaceID()));
        CFLog(INFO, "a - eastFace "<<eastFace.getGlobalFaceID()<<" has been matched with west "<<found->getGlobalFaceID()<<"\n");
        _faceConnectivityMap[eastFace.getLocalFaceID()*m_nbrFaceFlxPnts + eastFace.getLocalFluxID()] = found->getLocalFaceID();
        _fluxPointConnectivityMap[eastFace.getLocalFaceID()*m_nbrFaceFlxPnts + eastFace.getLocalFluxID()] = found->getLocalFluxID();
        _orientMap[eastFace.getLocalFaceID()*m_nbrFaceFlxPnts + eastFace.getLocalFluxID()] = found->getOrient();
        matches++;
     }else{
        std::vector<FlxPntStruct>::iterator found = std::find_if(westFaces.begin(), westFaces.end(), findTranslated(eastFace,translationVector,_threshold));
        if (found != westFaces.end()) {
           //if (_localConnectivityMap[found->getLocalFaceID()].size() == 0)   matches++;
           //_localConnectivityMap[found->getLocalFaceID()].push_back(PairStruct(iP,eastFace.getLocalFaceID()));
           CFLog(INFO, "b - eastFace "<<eastFace.getGlobalFaceID()<<" has been matched with west "<<found->getGlobalFaceID()<<"\n");
           _faceConnectivityMap[eastFace.getLocalFaceID()*m_nbrFaceFlxPnts + eastFace.getLocalFluxID()] = found->getLocalFaceID();
           _fluxPointConnectivityMap[eastFace.getLocalFaceID()*m_nbrFaceFlxPnts + eastFace.getLocalFluxID()] = found->getLocalFluxID();
           _orientMap[eastFace.getLocalFaceID()*m_nbrFaceFlxPnts + eastFace.getLocalFluxID()] = found->getOrient();
           matches++;
        }
     }
  } // end eastFace connectivity

  //} // End of nbProcesses loop

  //} // End if nbProcesses=1

   CFLog(INFO, "matches "<<matches<<"\n");
   CFLog(INFO, "END BCPeriodic::setup \n");

}

//////////////////////////////////////////////////////////////////////////////

void BCPeriodic::unsetup()
{
  CFAUTOTRACE;

  // unsetup of the parent class
  BCStateComputer::unsetup();
}
    
//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

