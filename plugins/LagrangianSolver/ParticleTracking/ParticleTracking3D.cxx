#include "ParticleTracking3D.hh"
#include <iostream>
#include <cmath>
#include <algorithm>
#include "Framework/MeshData.hh"

#include <sys/types.h>
#include <sys/wait.h> 
#include <unistd.h>

#define DEBUG 0

namespace COOLFluiD {

namespace LagrangianSolver {


ParticleTracking3D::ParticleTracking3D(const std::string& name):
    ParticleTracking(name),
    m_exitPoint(3),
    m_entryPoint(3),
    m_direction(3)
{

}

ParticleTracking3D::~ParticleTracking3D(){

}

void ParticleTracking3D::getCommonData(CommonData &data){
  static RealVector initialPoint(3);
  getExitPoint(initialPoint);

  data.currentPoint[0]=initialPoint[0];
  data.currentPoint[1]=initialPoint[1];
  data.currentPoint[2]=initialPoint[2];

  data.direction[0] = m_particleCommonData.direction[0];
  data.direction[1] = m_particleCommonData.direction[1];
  data.direction[2] = m_particleCommonData.direction[2];

  data.cellID = m_particleCommonData.cellID;
}

void ParticleTracking3D::myComputeCentroid( std::vector<Framework::Node*>& nodes, Vec3& centroid ){
   //average x y z
   centroid[0] = 0.;
   centroid[1] = 0.;
   centroid[2] = 0.;
   const CFuint nbNodes = nodes.size();
   for(CFuint i=0; i< nbNodes; ++i){
     centroid[0] += (*nodes[i])[0];
     centroid[1] += (*nodes[i])[1];
     centroid[2] += (*nodes[i])[2];
   }

  centroid[0] /= static_cast<CFreal>(nbNodes);
  centroid[1] /= static_cast<CFreal>(nbNodes);
  centroid[2] /= static_cast<CFreal>(nbNodes);
}

void ParticleTracking3D::setupAlgorithm(){
//    m_maxNbFaces = Framework::MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell();

//    CFuint nbFaces=Framework::MeshDataStack::getActive()->Statistics().getNbFaces();
//    m_centroids.reserve(nbFaces*3);

    using namespace Framework;
    using namespace Common;

//    CellTrsGeoBuilder::GeoData& cellsData = m_cellBuilder.getDataGE();
//    SafePtr<TopologicalRegionSet> MediumCells = MeshDataStack::getActive()->getTrs("InnerCells");
//    cellsData.trs = MediumCells;
//    const CFuint nbCellsMedia = MediumCells->getLocalNbGeoEnts();
//    for(CFuint i=0; i<nbCellsMedia; ++i){
//      cellsData.idx = i;
//      GeometricEntity *const cell = m_cellBuilder.buildGE();
//      CFuint nFaces = cell->nbNeighborGeos();
//      for(CFuint f=0; f<nFaces; ++f){
//        GeometricEntity* const face = cell->getNeighborGeo(f);
//        CFuint faceID = face->getID();
//        const RealVector centroid = face->computeCentroid();
//        m_centroids[faceID*3+0] = centroid[0];
//        m_centroids[faceID*3+1] = centroid[1];
//        m_centroids[faceID*3+2] = centroid[2];
//      }
//      m_cellBuilder.releaseGE();
//    }
}



#define EPSILON 1e-12

#define DOT(v1,v2) ( v1[0]*v2[0] +  v1[1]*v2[1] + v1[2]*v2[2] )  

#define CROSS(product, v1, v2) \
  product[0] = v1[1]*v2[2] - v1[2]*v2[1]  ; \
  product[1] = v1[2]*v2[0] - v1[0]*v2[2]  ; \
  product[2] = v1[0]*v2[1] - v1[1]*v2[0]  ;

#define SUB(dest, v1, v2 ) \
  dest[0] = v1[0] - v2[0]; \
  dest[1] = v1[1] - v2[1]; \
  dest[2] = v1[2] - v2[2];

//TODO: orient face 
bool triangle_intersection( const Vec3   V1,  // Triangle vertices
                            const Vec3   V2,
                            const Vec3   V3,
                            const Vec3    O,  //Ray origin
                            const Vec3    D,  //Ray direction
                                 CFreal* out )
{
  Vec3 e1, e2;  //Edge1, Edge2
  Vec3 P, Q, T;
  CFreal det, inv_det, u, v;
  CFreal t;
 
  //Find vectors for two edges sharing V1
  SUB(e1, V2, V1);
  SUB(e2, V3, V1);
  //Begin calculating determinant - also used to calculate u parameter
  CROSS(P, D, e2);
  //if determinant is near zero, ray lies in plane of triangle
  det = DOT(e1, P);
  //NOT CULLING
  if(det > -EPSILON && det < EPSILON) return 0;
  inv_det = 1.f / det;
 
  //calculate distance from V1 to ray origin
  SUB(T, O, V1);
 
  //Calculate u parameter and test bound
  u = DOT(T, P) * inv_det;
  //The intersection lies outside of the triangle
  if(u < 0.f || u > 1.f) return 0;
 
  //Prepare to test v parameter
  CROSS(Q, T, e1);
 
  //Calculate V parameter and test bound
  v = DOT(D, Q) * inv_det;
  //The intersection lies outside of the triangle
  if(v < 0.f || u + v  > 1.f) return 0;
 
  t = DOT(e2, Q) * inv_det;
 
  if(t > EPSILON) { //ray intersection
    *out = t;
    return true;
  }
 
  // No hit, no win
  return false;
}

 

void ParticleTracking3D::trackingStep(){

  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;

  //cout<<"%\n%Start tracking; CellID: "<<m_exitCellID<<endl<<"%\n";
  m_exitFaceID=-1;
  //m_exitCellID=-1;

  const RealVector &initialPoint = m_exitPoint;

  static DataHandle<CFint> faceIsOutwards= m_sockets.isOutward.getDataHandle();
  //static DataHandle<CFreal> normals= m_sockets.normals.getDataHandle();
  //static DataHandle<CFreal> faceCenters= m_sockets.faceCenters.getDataHandle();
  //vector<CFreal> &faceCenters = m_centroids;

  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  m_entryCellID = m_exitCellID;
  cellData.idx = m_entryCellID;

  GeometricEntity *const cell = m_cellBuilder.buildGE();
  CFuint nFaces = cell->getNbFacets();

  const Vec3 rayO = { m_exitPoint[0],
                      m_exitPoint[1],
                      m_exitPoint[2] };
  const Vec3 rayD = { m_direction[0],
                      m_direction[1],
                      m_direction[2] }; 

  #if DEBUG == 1
  stringstream rayData;
  stringstream faceData;
  rayData << "# X Y Z \n";
  rayData << rayO[0] << ' ' << rayO[1] << ' ' <<rayO[2] << '\n'; 
  faceData << "# X Y Z \n";
  #endif

  bool found = false;
  for(CFuint f=0; f<nFaces ; ++f){

    GeometricEntity* const face = cell->getNeighborGeo(f);
    
    vector<Node*>& faceNodes = *face->getNodes();
    RealVector centroid2(3);

    Vec3 centroid;
    myComputeCentroid(faceNodes, centroid);
    const CFuint nbNodes= faceNodes.size();
    const CFuint nbTris = nbNodes;
    const CFuint faceID = face->getID();
    CFreal outT;
    CFreal outTotalT = 0.;
    //find if the face normal point inside the cell
    //if so, we should reverse the circulation of the face points
    //so we have only outwards-facing normals    
    const bool reverseCirculation = ( static_cast<CFint>( faceIsOutwards[faceID] ) == m_entryCellID );
    
    for(CFuint iTri = 0; iTri < nbTris; ++iTri){
      //triangularize the faces using 2 points along the face and the centroid
      //Mega slow!
      CFuint iTri_1 =( iTri == nbTris-1 ) ? 0 : iTri+1;
      
      const Vec3 V1 = { (*faceNodes[iTri  ])[0],
      	                (*faceNodes[iTri  ])[1], 
                        (*faceNodes[iTri  ])[2] };

      const Vec3 V2 = { (*faceNodes[iTri_1])[0], 
        	        (*faceNodes[iTri_1])[1], 
         	        (*faceNodes[iTri_1])[2] }; 
      
      //reverse the triangle vertices to change the circulation and the normal
      const Vec3 &V1_reversed = reverseCirculation ? V2 : V1;
      const Vec3 &V2_reversed = reverseCirculation ? V1 : V2;
      #if DEBUG == 1 
      faceData << V1_reversed[0] << ' ' << V1_reversed[1] << ' ' << V1_reversed[2] << '\n';
      faceData << V2_reversed[0] << ' ' << V2_reversed[1] << ' ' << V2_reversed[2] << '\n';
      faceData <<    centroid[0] << ' ' <<    centroid[1] << ' ' <<    centroid[2] << '\n';
      #endif
         
      if( triangle_intersection(V1_reversed, V2_reversed, centroid, rayO, rayD, &outT)  ){
        outTotalT = found ? outTotalT : outT;
        outTotalT = (outT < outTotalT) ? outT : outTotalT; 
        //outTotalT = outT;
        found = true;
     }
 
    }
    #if DEBUG == 1
    faceData<<"\n \n";
    #endif

    if ( found ){
        m_stepDist = outTotalT;
        m_exitPoint = initialPoint + m_direction * m_stepDist;
     
        m_exitFaceID = faceID;
        m_exitCellID = face->getState(0)->getLocalID();

        m_exitCellID =  (m_exitCellID ==  m_entryCellID && !face->getState(1)->isGhost() ) ?
        face->getState(1)->getLocalID() : m_exitCellID;

        m_cellBuilder.releaseGE();
        #if DEBUG == 1
        rayData<< m_exitPoint[0] << ' ' << m_exitPoint[1] <<' ' << m_exitPoint[2] << '\n';
        break;
        #else
        return;
        #endif
    }
  }

  #if DEBUG == 1
  ofstream rayDataFile("rayData.dat",  std::ofstream::app);
  ofstream faceDataFile("faceData.dat", std::ofstream::app);
  
  rayDataFile << rayData.rdbuf();
  faceDataFile << faceData.rdbuf();
  
  faceDataFile.close();
  rayDataFile.close();
  #endif


  if ( !found ){ 
  #if DEBUG ==1 
  std::FILE* pipehandle=popen("gnuplot","w");
  std::fprintf(pipehandle,"set term x11 reset\n");
  std::fprintf(pipehandle,"splot \"faceData.dat\" with lines lc rgb \'black\' ,  \"rayData.dat\" with linespoints lc rgb \'red\' \n");

  std::fflush(pipehandle);
  std::cin.ignore();
  std::fprintf(pipehandle,"quit");
  std::fflush(pipehandle);
  std::fclose(pipehandle);
  #endif

    CFLog(VERBOSE, "ParticleTracking3D::trackingStep() => Can't find an exit Point!!\n");
  }
 
  m_cellBuilder.releaseGE();
}

void ParticleTracking3D::newParticle(CommonData &particle){

//  std::cout<<"%*******************\n%NEW PARTICLE\n%************************************\n";

  static RealVector buffer(3);
  ParticleTracking::newParticle(particle);

  m_entryCellID = m_particleCommonData.cellID;
  m_exitCellID = m_entryCellID;
  m_exitFaceID=-1;

  buffer[0] = m_particleCommonData.direction[0];
  buffer[1] = m_particleCommonData.direction[1];
  buffer[2] = m_particleCommonData.direction[2];

  m_exitPoint[0] = m_particleCommonData.currentPoint[0];
  m_exitPoint[1] = m_particleCommonData.currentPoint[1];
  m_exitPoint[2] = m_particleCommonData.currentPoint[2];

  newDirection(buffer);

  #if DEBUG == 1
  //clear the files
  std::ofstream rayDataFile("rayData.dat",  std::ofstream::trunc);
  std::ofstream faceDataFile("faceData.dat", std::ofstream::trunc);
  rayDataFile.close();
  faceDataFile.close();
  #endif

}


}
}
