#ifndef COOLFluiD_RadiativeTransfer_ParticleTracking2D_hh
#define COOLFluiD_RadiativeTransfer_ParticleTracking2D_hh

#include "ParticleTracking.hh"

namespace COOLFluiD {

namespace LagrangianSolver {

class ParticleTracking2D : public ParticleTracking
{
    public:

  void newDirection(RealVector &direction);

  void trackingStep();

  void getExitPoint(RealVector &exitPoint);

  CFreal getStepDistance();

  void setupAlgorithm();
  void newParticle(CommonData &particle);
  void getCommonData(CommonData &data);
  ParticleTracking2D(const std::string &name);
  ~ParticleTracking2D();

  void getNormals(CFuint faceID, RealVector &CartPosition, RealVector &faceNormal){
     getCartNormals(faceID, CartPosition, faceNormal);
  }

private:

  CFreal m_dx,m_dy,m_yy,m_A,m_xx,m_a,m_b,m_x0,m_y0,m_x1,m_x2,m_y1,m_y2;
  CFreal m_particle_t,m_particle_t_old, m_face_s,m_tt,m_ss, m_innerProd;
  RealVector faceOutNormal;
  RealVector particleTangent;

};


}

}
#endif
