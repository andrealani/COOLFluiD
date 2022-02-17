#include "Framework/MethodStrategyProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/BCNeumannFromFile.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"
#include "Common/OldLookupTable.hh"

#include <fstream>

#include "Framework/MeshData.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/TrsNotFoundException.hh"
#include "MathTools/MathChecks.hh"
#include "Common/SwapEmpty.hh"
#include "Common/FilesystemException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCNeumannFromFile,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionModule >
  BCNeumannFromFileProvider("NeumannFromFile");

//////////////////////////////////////////////////////////////////////////////

void BCNeumannFromFile::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<std::string> >
   ("TRSName","Name of the TRSs on which values must be prescribed");

  options.addConfigOption<std::string>
    ("FileNameTw", "Name of the file with the temperature distribution");

  options.addConfigOption<CFreal>
    ("RotationAngle", "Rotation angle (in degrees) to apply to the wall distribution coordinates");

  options.addConfigOption<std::vector<CFuint> >
    ("RotationCoordIDs", "IDs (must be 2) of the coordinates lying in the rotation plane");

  options.addConfigOption<bool>
    ("AddMinusConvention", "Add a minus sign to the read data due to a difference in convention"); 
  
  options.addConfigOption<CFuint>
    ("NbClosestPoints", "Number of closest points for surface interpolation"); 
  
  options.addConfigOption<vector<CFint> >
    ("ExtractCoordXYIDs", "IDs corresponding to the x,y coordinate (z=0) for which plane is extracted");
}

//////////////////////////////////////////////////////////////////////////////

BCNeumannFromFile::BCNeumannFromFile(const std::string& name) :
  BCStateComputer(name),
  m_extractCoordZID(-1),
  m_flxPntsLocalCoords(),
  m_allCellFlxPnts(),
  m_faceBuilder(),
  m_thisTRS(),
  m_flxLocalCoords(), 
  m_flxPntCoords(),
  m_nbrFaceFlxPnts(),
  m_dim(),      
  m_orient(),      
  m_intCell(),     
  m_faceFlxPntConn(),
  _faceBuilder(),          
  m_globalToLocalTRSFaceID(),
  m_flxPntTws(),
  m_nbGeoEnts(),
  m_nbEqs()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_trsName = vector<std::string>();
  this->setParameter("TRSName",&m_trsName);

  m_fileNameTw = "";
  this->setParameter("FileNameTw",&m_fileNameTw);

  m_angle = 0.0;
  this->setParameter("RotationAngle",&m_angle);

  m_xvec = std::vector<CFuint>();
  this->setParameter("RotationCoordIDs",&m_xvec);

  m_addMinus = false;
  this->setParameter("AddMinusConvention",&m_addMinus);
  
  m_nbClosestPoints = 0;
  this->setParameter("NbClosestPoints",&m_nbClosestPoints);
  
  m_extractCoordXYID = vector<CFint>();
  this->setParameter("ExtractCoordXYIDs",&m_extractCoordXYID);
}

//////////////////////////////////////////////////////////////////////////////

BCNeumannFromFile::~BCNeumannFromFile()
{
  CFAUTOTRACE;

}

//////////////////////////////////////////////////////////////////////////////

void BCNeumannFromFile::computeGhostStates(const vector< State* >& intStates,
                                     vector< State* >& ghostStates,
                                     const std::vector< RealVector >& normals,
                                     const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = intStates.size();
  cf_assert(nbrStates == ghostStates.size());
  cf_assert(nbrStates == normals.size());

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    *ghostStates[iState] = *intStates[iState];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCNeumannFromFile::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                        std::vector< std::vector< RealVector* > >& ghostGrads,
                                        const std::vector< RealVector >& normals,
                                        const std::vector< RealVector >& coords)
{
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // we prepare the face builder
  FaceToCellGEBuilder::GeoData& faceData = _faceBuilder->getDataGE();
  faceData.cellsTRS = cellTrs;
  faceData.facesTRS = m_thisTRS; 
  faceData.isBoundary = true;

  // Get the localFaceID from the map, knowing the faceGlobalID
  const CFuint faceLocalID = m_globalToLocalTRSFaceID.find(m_face->getID());

  // number of gradient variables
  cf_assert(nbrStateGrads > 0);
  const CFuint nbrGradVars = intGrads[0].size();

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    // normal
    const RealVector& normal = normals[iState];

    // We compute the local index of the flux point
    const CFuint localFlxPntID = faceLocalID*nbrStateGrads + iState; 

    for (CFuint iGradVar = 0; iGradVar < m_nbEqs; ++iGradVar)
    {
      const RealVector& varGradI =  *intGrads[iState][iGradVar];
      RealVector& varGradG =  *ghostGrads[iState][iGradVar];
      const CFreal nVarGrad = MathFunctions::innerProd(varGradI, normal);
      CFLog(VERBOSE, "face local ID: " << faceLocalID << ", local Flx ID: " << localFlxPntID << "\n");
      
      if (m_addMinus) 
      {
        varGradG = varGradI - 2.0*(m_flxPntTws[faceLocalID*nbrStateGrads + localFlxPntID] + nVarGrad)*normal;
      }
      else
      {
        varGradG = varGradI + 2.0*(m_flxPntTws[faceLocalID*nbrStateGrads + localFlxPntID] - nVarGrad)*normal;
      }
    }
  }

  _faceBuilder->releaseGE();
}

//////////////////////////////////////////////////////////////////////////////

void BCNeumannFromFile::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();
  
  m_dim = PhysicalModelStack::getActive()->getDim();
  if (m_nbClosestPoints == 0) {
    m_nbClosestPoints = m_dim + (m_dim - 2); // 2 in 2D, 4 in 3D
  }
  
  if (m_extractCoordXYID.size() > 0) 
  {
    if ((m_extractCoordXYID[0] == 0 && m_extractCoordXYID[1] == 1) ||
	  (m_extractCoordXYID[0] == 1 && m_extractCoordXYID[1] == 0)) 
    {
      m_extractCoordZID = 2;
    }
      
    if ((m_extractCoordXYID[0] == 0 && m_extractCoordXYID[1] == 2) ||
	  (m_extractCoordXYID[0] == 2 && m_extractCoordXYID[1] == 0)) 
    {
      m_extractCoordZID = 1;
    } 
      
    if ((m_extractCoordXYID[0] == 1 && m_extractCoordXYID[1] == 2) ||
	  (m_extractCoordXYID[0] == 2 && m_extractCoordXYID[1] == 1)) 
    {
      m_extractCoordZID = 0;
    }
  }
    
  if (std::abs(m_angle) > 0.0) 
  {
    // conversion to radiants
    m_angle *= MathConsts::CFrealPi()/180.;
    cf_always_assert(m_xvec.size() == 2);
  }
  
  vector<SurfaceData*> surfaces;
  readSurfaceData(surfaces);
  
  const CFuint nbSurf = surfaces.size();
  cf_assert(nbSurf >= 1);
  
  RealVector tmpNode(m_dim);
  
  m_faceBuilder = getMethodData().getFaceBuilder();
  _faceBuilder = getMethodData().getSecondFaceBuilder();
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);

  // get the face local coords of the flux points on one face
  m_flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  m_nbrFaceFlxPnts = m_flxLocalCoords->size();
  
  m_flxPntCoords.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx].resize(m_dim);
  }
  
  m_flxPntsLocalCoords.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  { 
    m_flxPntsLocalCoords[iFlx].resize(m_dim);
  }
  
  // get all flux points of a cell
  m_allCellFlxPnts = frLocalData[0]->getFlxPntsLocalCoords();
  
  // get TRS list
  vector< SafePtr< TopologicalRegionSet > > trsList = MeshDataStack::getActive()->getTrsList();
  
  const CFuint nbTRSs = trsList.size();
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    CFLog(INFO, "iTRS "<< iTRS << " \n");
    if (m_trsNames[0]==trsList[iTRS]->getName() ){
       m_thisTRS = trsList[iTRS];
       CFLog(INFO, "Matching BC "<<m_trsNames[0]<<" with "<<m_thisTRS->getName() << "\n");
    }
  }
  
  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& faceData = _faceBuilder->getDataGE();
  faceData.cellsTRS = cellTrs;
  faceData.facesTRS = m_thisTRS;
  faceData.isBoundary = true;

  // We get the number of faces in this TRS
  m_nbGeoEnts = m_thisTRS->getLocalNbGeoEnts();  //getNbTRs();

  // We get the number of equations
  m_nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  std::vector<FlxPntStruct> bndFlxPnts;
  bndFlxPnts.reserve(m_nbGeoEnts*m_nbrFaceFlxPnts); // oversized to avoid frequent memory reallocation

  FlxPntStruct thisFlxPnt;
  
  createFaceOrientationStartIndexes();
  
  map< std::string , vector< vector< CFuint > > >&
      bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
    vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[m_thisTRS->getName()];
    
    cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1; 

  const CFuint nbTRs = m_thisTRS->getNbTRs();

  // loop over TRs
  CFuint localFaceID = 0;
  bool isFirstOrientation = true;
  
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    for (m_orient = 0; m_orient < nbOrients; ++m_orient){
       // start and stop index of the faces with this orientation
       const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
       const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

       // loop over faces with this orientation
       for (CFuint iFace = startFaceIdx; iFace < stopFaceIdx; ++iFace){

          /* Face setup */
          faceData.idx = iFace;
          GeometricEntity *const face = _faceBuilder->buildGE();
          const CFuint faceGlobalID = face->getID();
          m_globalToLocalTRSFaceID.insert(faceGlobalID,localFaceID);
          
          for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
          {
            // Compute coordinates
             m_flxPntCoords[iFlx] = face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
          }

          for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
          {
             /* Load this flux point in a structure */
             thisFlxPnt.setCentreCoordinates(m_flxPntCoords[iFlx]);
             thisFlxPnt.setGlobalFaceID(faceGlobalID);
             thisFlxPnt.setLocalFaceID(localFaceID);
             thisFlxPnt.setLocalFluxID(iFlx);
             thisFlxPnt.setOrient(m_orient);
             
             bndFlxPnts.push_back(thisFlxPnt);
             
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
  
  CFuint nbBndFlxPnts = bndFlxPnts.size();
  
  ClosestPointData closestPoint;
  closestPoint.surfaceIDs.resize(m_nbClosestPoints); 
  closestPoint.pointsIDs.resize(m_nbClosestPoints); 
  closestPoint.r.resize(m_nbClosestPoints); 
  
  m_flxPntTws.resize(nbBndFlxPnts);
  
  const CFuint nbBnd = nbBndFlxPnts; 
  for(CFuint iFlx=0; iFlx<nbBnd; iFlx++) 
  {
    FlxPntStruct& currFlxPnt = bndFlxPnts[iFlx];

    // during this preprocessing m_nbClosestPoints closest neighbors are sought
    closestPoint.reset();
      
    bool flagOut = false;
    for (CFuint is = 0; is < nbSurf && (!flagOut); ++is) 
    {
      const SurfaceData& sf = *surfaces[is];
      const CFuint nbPoints = sf.Tw.size();
      for (CFuint ip = 0; ip < nbPoints; ++ip) 
      {
	sf.xyz.putRow(ip,tmpNode);

        const CFreal distance = MathFunctions::getDistance(currFlxPnt.getCentreCoordinates(), tmpNode);

        if (distance < 1.0e-8) 
        {
	  // in this case we assume that the current node coincides with the mapping node
	  // set the matching Twall in the array
	  //cf_assert(iNode < nodalTwall.size());
	  cf_assert(ip < sf.Tw.size());
	  //nodalTwall[iNode] = sf.Tw[ip]; 
          currFlxPnt.setTw(sf.Tw[ip]);
          m_flxPntTws[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPnts + currFlxPnt.getLocalFluxID()] = sf.Tw[ip];
	  flagOut = true;
	  break;
	}
	else 
        {
	  CFint counter = -1;
	  for (CFuint n = 0; n < m_nbClosestPoints; ++n) 
          {
	    cf_assert(n < closestPoint.r.size());
	    if (distance < closestPoint.r[n]) 
            {
	      counter++;
	    }
	  }
	    
	  if (counter >= 0) 
          {
	      // counter == (m_nbClosestPoints-1) corresponds to the point with min distance within the stencil
	      // counter == 0                     corresponds to the point with max distance within the stencil
	    for (CFuint i = 0; i < static_cast<CFuint>(counter); ++i) 
            {
              closestPoint.regressionFromTo(i+1, i);
	    }
	      
	    cf_assert(counter < closestPoint.surfaceIDs.size());

            closestPoint.surfaceIDs[counter] = is;
	    cf_assert(counter < closestPoint.pointsIDs.size());
	    closestPoint.pointsIDs[counter] = ip;
	    cf_assert(counter < closestPoint.r.size());
	    closestPoint.r[counter] = distance;
	  }
	}
      }
    }

    if (!flagOut) 
    {
      CFreal matchingTw = 0.0;
      CFreal sumWeights = 0.0;
      
      for (CFuint n = 0; n < m_nbClosestPoints; ++n) 
      {
	const CFuint idxs = closestPoint.surfaceIDs[n];
	cf_assert(idxs < surfaces.size());
	const SurfaceData& sf = *surfaces[idxs];
	cf_assert(closestPoint.r[n] > 0.);
	const CFreal weight = 1./closestPoint.r[n]; 
	sumWeights += weight;
	const CFuint idxp = closestPoint.pointsIDs[n];
	cf_assert(idxp < sf.Tw.size());
	matchingTw += weight*sf.Tw[idxp]; 
      }
      
      matchingTw /= sumWeights;
	
      // set the matching Twall in the array
      //nodalTwall[iNode] = matchingTw;
      currFlxPnt.setTw(matchingTw);
      m_flxPntTws[currFlxPnt.getLocalFaceID()*m_nbrFaceFlxPnts + currFlxPnt.getLocalFluxID()] = matchingTw;
    }
  }
  
  // cleanup the memory
  for (CFuint is = 0; is < nbSurf; ++is) {
    deletePtr(surfaces[is]);
  }
}

//////////////////////////////////////////////////////////////////////////////


void BCNeumannFromFile::createFaceOrientationStartIndexes()
{
  CFAUTOTRACE;
   //MPI_Barrier(_comm);//the processing of a an individual processor is paused temporarily till all the processors compute till some certain values

   //MPI_Allgather(&m_nbGeoEnts, 1, MPIStructDef::getMPIType(&m_nbGeoEnts), 
	//&_nbFacesPerProcess[0], 1, MPIStructDef::getMPIType(&m_nbGeoEnts), _comm);

  CFLog(INFO, " BCNeumannFromFile::createFaceOrientationStartIndexes()\n");

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

void BCNeumannFromFile::readSurfaceData(std::vector<SurfaceData*>& surfaces)
{ 
  const CFuint dim = (m_extractCoordZID < 0) ? PhysicalModelStack::getActive()->getDim() : 3;

  boost::filesystem::path fname;
  if (Common::StringOps::startsWith(m_fileNameTw,"."))
  {
    fname = boost::filesystem::path(m_fileNameTw);
  }
  else
  {
    fname = Environment::DirPaths::getInstance().getBaseDir() / boost::filesystem::path(m_fileNameTw);
  }

  // check the file type
  std::ifstream fin(fname.string().c_str());
  if(!fin) throw Common::FilesystemException (FromHere(),"Could not open file: " + fname.string());
  // The format is as follows:
  
  // NumberOfSurfaces
  // SURFACE_NAME1 NumberOfPoints
  // x y z T
  // ...
  // SURFACE_NAME2 NumberOfPoints
  // x y z T
  // ...
  
  CFLogInfo("BCNeumannFromFile::readSurfaceData() => START reading file " <<
	    m_fileNameTw << "\n");
  
  CFuint nbSurf = 0;
  fin >> nbSurf;
  surfaces.resize(nbSurf);
  
  // store all the surface data
  for (CFuint is = 0; is < nbSurf; ++is) {
    std::string nameSurf = "";
    fin >> nameSurf;
    CFuint nbPoints = 0;
    fin >> nbPoints;
    
    CFLogInfo("nbSurf = " << is << "/" << nbSurf << ", NameSurf = " <<
	      nameSurf << ", nbPoints = " << nbPoints << "\n");
    SurfaceData* sf = new SurfaceData();
    sf->xyz.resize(nbPoints,dim);
    sf->Tw.resize(nbPoints);
    
    for (CFuint ip = 0; ip < nbPoints; ++ip) {
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
	fin >> sf->xyz(ip,iDim);
      }
      
      if (m_xvec.size() > 0) {
	cf_assert(std::abs(m_angle) > 0.0);
	
	// rotate coordinates
	const CFreal x1 =  sf->xyz(ip,m_xvec[0])*std::cos(m_angle) + sf->xyz(ip,m_xvec[1])*std::sin(m_angle);
	const CFreal y1 = -sf->xyz(ip,m_xvec[0])*std::sin(m_angle) + sf->xyz(ip,m_xvec[1])*std::	cos(m_angle);
	sf->xyz(ip,m_xvec[0]) = x1;				   
	sf->xyz(ip,m_xvec[1]) = y1;				   
      }	
      fin >> sf->Tw[ip];
    }
	
    surfaces[is] = (m_extractCoordZID < 0) ?  sf : extractLineData(sf);
  }
  
  CFLogInfo("BCNeumannFromFile::readSurfaceData() => END reading file " <<
	    m_fileNameTw << "\n"); 
}

//////////////////////////////////////////////////////////////////////////////

typename BCNeumannFromFile::SurfaceData* 
BCNeumannFromFile::extractLineData(SurfaceData* surface)
{
  // if m_extractCoordZID >= 0 a line distribution is extracted from a {x y z T} distribution for z=0
  //                           with z = x[m_extractCoordZID]
  
  CFreal minZ = MathConsts::CFrealMax();
  const CFuint nbPoints = surface->Tw.size();
  for (CFuint ip = 0; ip < nbPoints; ++ip) {
    minZ = std::min(std::abs(surface->xyz(ip, m_extractCoordZID)), minZ);
  }
  
  CFuint nbPointsOnLine = 0;
  for (CFuint ip = 0; ip < nbPoints; ++ip) {
    if (MathChecks::isEqual(std::abs(surface->xyz(ip, m_extractCoordZID)), minZ)) {
      nbPointsOnLine++;
    }
  }
  
  SurfaceData* newData = new SurfaceData(); 
  newData->xyz.resize(nbPointsOnLine,2);
  newData->Tw.resize(nbPointsOnLine);
  
  CFuint counter = 0;
  for (CFuint ip = 0; ip < nbPoints; ++ip) {
    if (MathChecks::isEqual(std::abs(surface->xyz(ip, m_extractCoordZID)), minZ)) {
      cout << surface->Tw[ip] << endl;
      newData->xyz(counter, XX) = surface->xyz(ip, m_extractCoordXYID[XX]);
      newData->xyz(counter, YY) = surface->xyz(ip, m_extractCoordXYID[YY]);
      newData->Tw[counter] = surface->Tw[ip];
      counter++;
    }
  }
  
  cf_always_assert(counter == nbPointsOnLine);
  
  // remove old SurfaceData
  deletePtr(surface);
  
  //replace with new SurfaceData
  return newData;  
}

//////////////////////////////////////////////////////////////////////////////

void BCNeumannFromFile::unsetup()
{
  CFAUTOTRACE;

  // unsetup of the parent class
  BCStateComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
