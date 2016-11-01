#include "LagrangianSolver/ParticleTracking/ParticleTracking3D.hh"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sys/types.h>
#include <sys/wait.h> 
#include <unistd.h>

#include "Framework/MeshData.hh"

#define DEBUG 0

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace LagrangianSolver {

//////////////////////////////////////////////////////////////////////////////

ParticleTracking3D::ParticleTracking3D(const std::string& name):
  ParticleTracking(name),
  m_faceIsOutwards(CFNULL),
  m_exitPoint(3),
  m_entryPoint(3),
  m_direction(3),
  m_initialPoint(3),
  m_buffer(3),
  m_centroid(),
  m_v1(),
  m_v2(),
  m_e1(),
  m_e2(),
  m_P(), 
  m_Q(),
  m_T(),
  m_rayO(),
  m_rayD()
{
}

//////////////////////////////////////////////////////////////////////////////

ParticleTracking3D::~ParticleTracking3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParticleTracking3D::getCommonData(CommonData &data)
{
  getExitPoint(m_initialPoint);
  data.currentPoint = m_initialPoint;
  data.direction = m_particleCommonData.direction;
  data.cellID = m_particleCommonData.cellID;
}
  
//////////////////////////////////////////////////////////////////////////////
  
inline void ParticleTracking3D::myComputeCentroid(const std::vector<Framework::Node*>& nodes, 
						  Vec3& centroid)
{
  //average x y z
  centroid = 0.;
  const CFuint nbNodes = nodes.size();
  for (CFuint i=0; i< nbNodes; ++i) {
    centroid += (*nodes[i]);
  }
  centroid /= static_cast<CFreal>(nbNodes);
}
  
//////////////////////////////////////////////////////////////////////////////

void ParticleTracking3D::setupAlgorithm()
{
  m_faceIsOutwards = m_sockets.isOutward.getDataHandle();
}  
  
//////////////////////////////////////////////////////////////////////////////
  
//TODO: orient face 
bool ParticleTracking3D::triangleIntersection(const Vec3& V1, // Triangle vertices
					      const Vec3& V2,
					      const Vec3& V3,
					      const Vec3& O,  //Ray origin
					      const Vec3& D,  //Ray direction
					      CFreal* out)
{
  using namespace MathTools;
  
  const CFreal EPSILON = 1e-12;
  CFreal det = 0.;
  CFreal inv_det = 0.;
  CFreal u = 0.;
  CFreal v = 0.;
  CFreal t = 0.;
  
  //Find vectors for two edges sharing V1
  m_e1 = V2 - V1;
  m_e2 = V3 - V1;
  //Begin calculating determinant - also used to calculate u parameter
  MathFunctions::crossProd(D, m_e2, m_P);
  //if determinant is near zero, ray lies in plane of triangle
  det = MathFunctions::innerProd(m_e1, m_P);
  //NOT CULLING
  if(det > -EPSILON && det < EPSILON) return 0;
  inv_det = 1.f / det;
  
  //calculate distance from V1 to ray origin
  m_T = O - V1;
  
  //Calculate u parameter and test bound
  u = MathFunctions::innerProd(m_T, m_P) * inv_det;
  //The intersection lies outside of the triangle
  if(u < 0.f || u > 1.f) return 0;
  
  //Prepare to test v parameter
  MathFunctions::crossProd(m_T, m_e1, m_Q);
  
  //Calculate V parameter and test bound
  v = MathFunctions::innerProd(D, m_Q) * inv_det;
  //The intersection lies outside of the triangle
  if(v < 0.f || u + v  > 1.f) return 0;
 
  t = MathFunctions::innerProd(m_e2, m_Q) * inv_det;
 
  if(t > EPSILON) { //ray intersection
    *out = t;
    return true;
  }
  
  // No hit, no win
  return false;
}

//////////////////////////////////////////////////////////////////////////////

void ParticleTracking3D::trackingStep()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;
  
  CFLog(DEBUG_MAX, "ParticleTracking3D::trackingStep() => Start tracking CellID: "<< m_exitCellID << "\n");
  
  m_exitFaceID = -1;
  
  const RealVector &initialPoint = m_exitPoint;
  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  m_entryCellID = m_exitCellID;
  cellData.idx = m_entryCellID;
  
  GeometricEntity *const cell = m_cellBuilder.buildGE();
  const CFuint nbFaces = cell->getNbFacets();
  
  m_rayO = m_exitPoint;
  m_rayD = m_direction;
  
#if DEBUG == 1
  stringstream rayData;
  stringstream faceData;
  rayData << "# X Y Z \n";
  rayData << m_rayO << '\n'; 
  faceData << "# X Y Z \n";
#endif

  bool found = false;
  for(CFuint f=0; f< nbFaces ; ++f){
    GeometricEntity* const face = cell->getNeighborGeo(f);
    const vector<Node*>& faceNodes = *face->getNodes();
    myComputeCentroid(faceNodes, m_centroid);
    const CFuint faceID = face->getID();
    CFreal outT = 0.;
    CFreal outTotalT = 0.;
    //find if the face normal point inside the cell
    //if so, we should reverse the circulation of the face points
    //so we have only outwards-facing normals    
    const bool reverseCirculation = ((CFint)(m_faceIsOutwards[faceID]) == m_entryCellID);
    const CFuint nbTris = faceNodes.size();
    for (CFuint iTri = 0; iTri < nbTris; ++iTri) {
      //triangularize the faces using 2 points along the face and the centroid
      //Mega slow!
      const CFuint iTri_1 =( iTri == nbTris-1 ) ? 0 : iTri+1;
      m_v1 = (*faceNodes[iTri]);
      m_v2 = (*faceNodes[iTri_1]);
      
      //reverse the triangle vertices to change the circulation and the normal
      const Vec3& V1_reversed = reverseCirculation ? m_v2 : m_v1;
      const Vec3& V2_reversed = reverseCirculation ? m_v1 : m_v2;
#if DEBUG == 1 
      faceData << V1_reversed << '\n';
      faceData << V2_reversed << '\n';
      faceData << m_centroid  << '\n';
#endif
      
      if (triangleIntersection(V1_reversed, V2_reversed, m_centroid, m_rayO, m_rayD, &outT)) {
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
      rayData << m_exitPoint << '\n';
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
  
//////////////////////////////////////////////////////////////////////////////

void ParticleTracking3D::newParticle(CommonData &particle)
{
  ParticleTracking::newParticle(particle);

  m_entryCellID = m_particleCommonData.cellID;
  m_exitCellID = m_entryCellID;
  m_exitFaceID = -1;
  
  m_buffer = m_particleCommonData.direction;
  m_exitPoint = m_particleCommonData.currentPoint;
  
  newDirection(m_buffer);
  
#if DEBUG == 1
  //clear the files
  std::ofstream rayDataFile("rayData.dat",  std::ofstream::trunc);
  std::ofstream faceDataFile("faceData.dat", std::ofstream::trunc);
  rayDataFile.close();
  faceDataFile.close();
#endif
}

//////////////////////////////////////////////////////////////////////////////

}
}
