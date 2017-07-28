#ifndef COOLFluiD_RadiativeTransfer_ParticleTracking3D_hh
#define COOLFluiD_RadiativeTransfer_ParticleTracking3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "ParticleTracking.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace LagrangianSolver {

//////////////////////////////////////////////////////////////////////////////

typedef CFreal Vec3[3];
  
class ParticleTracking3D : public ParticleTracking
{
public:
  ParticleTracking3D(const std::string& name);
  
  ~ParticleTracking3D();
  
  void setupAlgorithm();
  
  void getCommonData(CommonData &data);
  
  void newDirection(RealVector &direction){
    m_direction = direction;
  }
  
  void trackingStep();
  
  void getExitPoint(RealVector &exitPoint) {exitPoint = m_exitPoint;}
  
  void newParticle(CommonData &particle);
  
  void getNormals(CFuint faceID, RealVector &CartPosition, RealVector &faceNormal)
  {
    getCartNormals(faceID, CartPosition, faceNormal);
  }
  
  CFreal getStepDistance() {return m_stepDist;}

private:
  
  /// compute the triangle intersection
  /// @param V1   triangle vertex 1
  /// @param V2   triangle vertex 2
  /// @param V3   triangle vertex 3 (centroid)
  /// @param O    ray origin
  /// @param D    ray direction
  /// @param out  ray intersection
  bool triangle_intersection(const Vec3& V1, const Vec3& V2, 
			     const Vec3& V3, const Vec3& O, 
			     const Vec3& D, CFreal* out);
  
private:
  
  Framework::DataHandle<CFint> m_isOutward; 
  Framework::DataHandle<CFreal> m_faceCenters; 
  std::vector<CFreal> m_centroids;
  CFuint m_maxNbFaces;
  RealVector m_exitPoint;
  RealVector m_entryPoint;
  RealVector m_direction;
  RealVector m_initialPoint;
  CFreal m_stepDist;
};

//////////////////////////////////////////////////////////////////////////////

}
  
}

//////////////////////////////////////////////////////////////////////////////

#endif
