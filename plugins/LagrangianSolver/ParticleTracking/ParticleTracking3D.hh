#ifndef COOLFluiD_LagrangianSolver_ParticleTracking3D_hh
#define COOLFluiD_LagrangianSolver_ParticleTracking3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "LagrangianSolver/ParticleTracking/ParticleTracking.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace LagrangianSolver {

//////////////////////////////////////////////////////////////////////////////

typedef MathTools::CFVec<CFreal,3> Vec3;

//////////////////////////////////////////////////////////////////////////////

class ParticleTracking3D : public ParticleTracking {
public:
  ParticleTracking3D(const std::string& name);
  
  ~ParticleTracking3D();
  
  void setupAlgorithm();
  
  void getCommonData(CommonData &data);
  
  void newDirection(RealVector &direction) {m_direction = direction;}
  
  void trackingStep();
  
  void getExitPoint(RealVector &exitPoint) {exitPoint = m_exitPoint;}
  
  void newParticle(CommonData &particle);
  
  void getNormals(CFuint faceID, RealVector &CartPosition, RealVector &faceNormal){
    getCartNormals(faceID, CartPosition, faceNormal);
  }
  
  void myComputeCentroid( const std::vector<Framework::Node*>& nodes, Vec3& centroid );
  
  CFreal getStepDistance() {return m_stepDist;}
  
  bool triangleIntersection(const Vec3&   V1,  // Triangle vertices
			    const Vec3&   V2,
			    const Vec3&   V3,
			    const Vec3&    O,  //Ray origin
			    const Vec3&    D,  //Ray direction
			    CFreal* out);
  
private:
  
  /// array telling for each faceID the neighbor cellID  for which the normal is outward
  Framework::DataHandle<CFint> m_faceIsOutwards;
  
  std::vector<CFreal> m_centroids;
  CFuint m_maxNbFaces;
  RealVector m_exitPoint;
  RealVector m_entryPoint;
  RealVector m_direction;
  RealVector m_initialPoint;
  RealVector m_buffer;
  Vec3 m_centroid;
  Vec3 m_v1;
  Vec3 m_v2;
  Vec3 m_rayO;
  Vec3 m_rayD;
  Vec3 m_e1;
  Vec3 m_e2;
  Vec3 m_P; 
  Vec3 m_Q;
  Vec3 m_T;
  CFreal m_stepDist;
};
  
//////////////////////////////////////////////////////////////////////////////

}// namespace LagrangianSolver

}// namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
