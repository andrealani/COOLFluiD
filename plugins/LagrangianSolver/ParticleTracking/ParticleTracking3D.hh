#ifndef COOLFluiD_RadiativeTransfer_ParticleTracking3D_hh
#define COOLFluiD_RadiativeTransfer_ParticleTracking3D_hh

#include "ParticleTracking.hh"

namespace COOLFluiD {

namespace LagrangianSolver {

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

    void getExitPoint(RealVector &exitPoint){
       exitPoint = m_exitPoint;
    }

    void newParticle(CommonData &particle);

    void getNormals(CFuint faceID, RealVector &CartPosition, RealVector &faceNormal){
       getCartNormals(faceID, CartPosition, faceNormal);
    }

    CFreal getStepDistance(){
     return m_stepDist;
    }

private:
    std::vector<CFreal> m_centroids;
    CFuint m_maxNbFaces;
    RealVector m_exitPoint, m_entryPoint, m_direction;
    CFreal m_stepDist;
};



}

}
#endif
