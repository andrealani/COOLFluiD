#include "LagrangianSolver/ParallelVector/ParallelVector.hh"
#include "ParticleTracking.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

////////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace LagrangianSolver {

//////////////////////////////////////////////////////////////////////////////

ParticleTracking::ParticleTracking(const std::string& name) :
  SocketBundleSetter(),
  m_particleDataType(),
  m_particleCommonData(),
  m_exitFaceID(0),
  m_entryFaceID(0),
  m_exitCellID(0),
  m_entryCellID(0),
  m_cellIdx(0),
  m_faceIdx(0),
  m_dim(2), 
  m_normals(CFNULL),
  m_cartNormal(2)
{
}

//////////////////////////////////////////////////////////////////////////////

ParticleTracking::~ParticleTracking()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParticleTracking::setupAlgorithm()
{
  m_normals = m_sockets.normals.getDataHandle();
  m_dim = Framework::PhysicalModelStack::getActive()->getDim();
}

//////////////////////////////////////////////////////////////////////////////

void ParticleTracking::getAxiNormals(CFuint faceID,
				     RealVector& CartPosition, 
				     RealVector& faceNormal)
{
  faceNormal.resize(3);
  
  //Create the Cylindrical normal vector
  
  //1. Get the 2D cartesian normal
  getCartNormals(faceID, CartPosition, m_cartNormal);
  
  CFLog(DEBUG_MAX, "ParticleTracking::getAxiNormals() => cartNormal = " 
	<< m_cartNormal[0] << " " << m_cartNormal[1]<< "\n");
  
  //2. Convert them in Cilindrical coordinates
  const CFreal c1 = CartPosition[1];
  const CFreal c2 = CartPosition[2];
  const CFreal invNorm = 1./std::sqrt(c1*c1 + c2*c2);
  
  faceNormal[XX] = m_cartNormal[0];
  faceNormal[YY] = m_cartNormal[1]*c1 * invNorm;
  faceNormal[ZZ] = m_cartNormal[1]*c2 * invNorm;
}
  
/////////////////////////////////////////////////////////////////////////////
  
void ParticleTracking::getCartNormals(CFuint faceID, 
				      RealVector& CartPosition, 
				      RealVector& faceNormal)
{
  const CFuint dim = m_dim;
  const CFuint startID = faceID*dim;
  DataHandle<CFreal> normals = m_normals;
  cf_assert(startID < normals.size() );
  // outwards pointing normals (AL: are you sure about this???)
  for(CFuint i=0; i<dim;++i){
    faceNormal[i]= -normals[startID+i];
  }
  faceNormal.normalize();
}

/////////////////////////////////////////////////////////////////////////////

void ParticleTracking::setFaceTypes(MathTools::CFMat<CFint> & wallTypes, 
				    std::vector<std::string>& wallNames,
                                    std::vector<std::string>& boundaryNames)
{
  setupAlgorithm();
  
  FaceTrsGeoBuilder::GeoData& facesData = m_faceBuilder.getDataGE();
  CellTrsGeoBuilder::GeoData& cellsData = m_cellBuilder.getDataGE();
  
  LagrangianSolver::ParallelVector < CFuint > ownerCellRank, ownerCellLocalID;

  ownerCellRank.setDataSockets(m_sockets);
  ownerCellLocalID.setDataSockets(m_sockets);

  ownerCellRank.getSharedEntries();
  ownerCellLocalID.getSharedEntries();
  
  const CFuint nbCells = cellsData.trs->getLocalNbGeoEnts();
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  
  //vector<CFuint> ownerCellRank(nbCells,m_myProcessRank);
  CFuint myProcessRank = PE::GetPE().GetRank(nsp);
  ownerCellRank.resize(nbCells,myProcessRank);
  ownerCellLocalID.resize(nbCells);
  
  //fill with the local iD's
  for(CFuint i=0; i<nbCells; ++i){
    cellsData.idx = i;
    GeometricEntity *const cell = m_cellBuilder.buildGE();
    ownerCellLocalID[i]=cell->getState(0)->getLocalID();
    m_cellBuilder.releaseGE();
  }
  
  ownerCellRank.sincronizeAssign();
  ownerCellLocalID.sincronizeAssign();
  
  //1. Get the total number of faces.
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  //cout<<"number of faces= "<< nbFaces<<endl;

  wallTypes.resize(nbFaces, 4, -1);

  SafePtr<TopologicalRegionSet> ptrTrs;

  std::vector<std::string> mediaNames;
  mediaNames.push_back("InnerCells");

  //Set type and physics ID for the cells
  CFuint nbMediaTrs = mediaNames.size();
  for(CFuint j=0; j<nbMediaTrs; ++j){
    SafePtr<TopologicalRegionSet> MediumCells = MeshDataStack::getActive()->getTrs(mediaNames[j]);
    cellsData.trs = MediumCells;
    const CFuint nbCellsMedia = MediumCells->getLocalNbGeoEnts();
    for(CFuint i=0; i<nbCellsMedia; ++i){
      cellsData.idx = i;
      GeometricEntity *const cell = m_cellBuilder.buildGE();
      if(!cell->getState(0)->isGhost() ){
        CFuint nFaces = cell->nbNeighborGeos();
        for(CFuint f=0; f<nFaces; ++f){
          GeometricEntity* const face = cell->getNeighborGeo(f);
          wallTypes(face->getID(),0) = INTERNAL_FACE;
          wallTypes(face->getID(),1) = cell->getState(0)->getLocalID();
          //wallTypes(face->getID(),2) = cell->getState(1)->getLocalID();
          wallTypes(face->getID(),2) = ownerCellRank[cell->getState(0)->getLocalID() ];
          wallTypes(face->getID(),3) = ownerCellLocalID[cell->getState(0)->getLocalID() ];
        }
      }
      m_cellBuilder.releaseGE();
    }
  }


  //Set type and physics ID for the walls
  m_faceBuilder.getDataGE().isBFace = true;
  CFuint nbWallTrs = wallNames.size();
  for(CFuint j=0; j<nbWallTrs; ++j){
    SafePtr<TopologicalRegionSet> WallFaces = MeshDataStack::getActive()->getTrs(wallNames[j]);
    facesData.trs = WallFaces;
    const CFuint nbFacesWall = WallFaces->getLocalNbGeoEnts();
    for(CFuint i=0; i<nbFacesWall; ++i){
      facesData.idx = i;
      GeometricEntity *const face = m_faceBuilder.buildGE();

      bool isState0G= face->getState(0)->isGhost();
      bool isState1G= face->getState(1)->isGhost();

      if(isState0G ^ isState1G){ //XOR
        //CFuint s = (isState0G ? 1 : 0);
        wallTypes(face->getID(),0) = WALL_FACE;
        wallTypes(face->getID(),1) = face->getState(0)->getLocalID();
        //wallTypes(face->getID(),2) = ownerCellRank[face->getState(s)->getLocalID() ];
        //wallTypes(face->getID(),3) = ownerCellLocalID[face->getState(s)->getLocalID() ];
        wallTypes(face->getID(),2) = face->getState(1)->getLocalID();
      }

//      if(face->getState(1)->isGhost()){
//        wallTypes(face->getID(),1) = face->getState(1)->getLocalID();
//      }

      m_faceBuilder.releaseGE();
    }
  }

  //Set type for the boundaries
  m_faceBuilder.getDataGE().isBFace = true;
  const CFuint nbBoundariesTrs = boundaryNames.size();
  for(CFuint j=0; j<nbBoundariesTrs; ++j){
    SafePtr<TopologicalRegionSet> WallFaces = MeshDataStack::getActive()->getTrs(boundaryNames[j]);
    facesData.trs = WallFaces;
    const CFuint nbFacesWall = WallFaces->getLocalNbGeoEnts();
    for(CFuint i=0; i<nbFacesWall; ++i){
      facesData.idx = i;
      GeometricEntity *const face = m_faceBuilder.buildGE();
      const bool isState0G= face->getState(0)->isGhost();
      const bool isState1G= face->getState(1)->isGhost();
      
      if(isState0G ^ isState1G){ //XOR
        //CFuint s = (isState0G ? 1 : 0);
        wallTypes(face->getID(),0) = BOUNDARY_FACE;
        wallTypes(face->getID(),1) = face->getState(0)->getLocalID();
        //wallTypes(face->getID(),2) = ownerCellRank[face->getState(s)->getLocalID() ];
        //wallTypes(face->getID(),3) = ownerCellLocalID[face->getState(s)->getLocalID() ];
        wallTypes(face->getID(),2) = face->getState(1)->getLocalID();
      }
//      if(face->getState(1)->isGhost()){
//        wallTypes(face->getID(),1) = face->getState(1)->getLocalID();
//      }

      m_faceBuilder.releaseGE();
    }
  }


  //Set type for the Partition Faces

//  CFLog(INFO,"Cell ID, Cell Rank, LocalID: \n");
//  for(CFuint i=0; i<nbCells; ++i){
//      CFLog(INFO,i<<' '<<ownerCellRank[i]<<' '<<ownerCellLocalID[i]<<"\n");
//  }

  ptrTrs = MeshDataStack::getActive()->getTrs("InnerFaces");

  facesData.trs = ptrTrs;
  m_faceBuilder.getDataGE().isBFace = false;
  CFuint nbInnerFaces = ptrTrs->getLocalNbGeoEnts();
  for(CFuint i=0; i< nbInnerFaces; ++i){
    facesData.idx = i;
    GeometricEntity *face = m_faceBuilder.buildGE();

    bool isState0G= face->getState(0)->isParUpdatable();
    bool isState1G= face->getState(1)->isParUpdatable();

    if(isState0G ^ isState1G){ //XOR
      CFuint s = (isState0G ? 1 : 0);
      wallTypes( face->getID(),0) = COMP_DOMAIN_FACE;
      wallTypes( face->getID(),1) = face->getState(0)->getLocalID();
      wallTypes( face->getID(),2) = ownerCellRank[ face->getState(s)->getLocalID() ];
      wallTypes( face->getID(),3) = ownerCellLocalID[ face->getState(s)->getLocalID() ];
    }
    m_faceBuilder.releaseGE();
  }

//  CFLog(INFO,"Matrix: \n");
//  for(CFuint i=0; i<wallTypes.nbRows();++i){
//    CFLog(INFO,wallTypes(i,0)<<' '<<wallTypes(i,1)<<' '<<wallTypes(i,2)<<' '<<wallTypes(i,3)<<'\n');
//  }
}

//////////////////////////////////////////////////////////////////////////////

}
}
