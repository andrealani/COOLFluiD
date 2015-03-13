#ifndef COOLFluiD_RadiativeTransfer_ParticleTrackingAxi_hh
#define COOLFluiD_RadiativeTransfer_ParticleTrackingAxi_hh

#include <iostream>
#include <cmath>
#include <algorithm>

#include "ParticleTracking.hh"

namespace COOLFluiD {

namespace LagrangianSolver {

class ParticleTrackingAxi : public ParticleTracking
{
public:
    ParticleTrackingAxi(const std::string& name);

    ~ParticleTrackingAxi();

    void getCommonData(CommonData &data);

    void newDirection(RealVector &direction);

    void trackingStep();

    void getExitPoint(RealVector &exitPoint);

    void newParticle(CommonData &particle);

    void setupAlgorithm();

    void getNormals(CFuint faceID, RealVector &CartPosition, RealVector &faceNormal){
       getAxiNormals(faceID, CartPosition, faceNormal);
    }

    CFreal getStepDistance();

private:

    CFuint m_maxNbFaces;
    CFreal m_c,m_b,m_a,m_z0,m_y0,m_x0,m_invA,m_D,m_D1,m_D2,m_B1,m_B2,m_z1,m_z2,m_r1,m_r2,m_dr,m_dz,m_zz,m_dr_2,m_dz_2,m_t[2],m_s[2];
    bool isS1,isS2;
    CFreal m_D22, m_cx0, m_cy0,m_a1,m_a2,m_ca2,m_c_2;
    CFreal m_particle_t, m_particle_t_old;

    //std::vector<CFreal> m_tCantidates;
    //std::vector<CFreal> m_fCandidates;
    RealVector faceOutNormal;
    RealVector rayTangent;

};



}

}
#endif
