#ifndef COOLFluiD_RadiativeTransfer_ParticleTracking_hh
#define COOLFluiD_RadiativeTransfer_ParticleTracking_hh

//////////////////////////////////////////////////////////////////////////////

//#include "RadiativeTransfer/RadiativeTransfer.hh"

#include "Framework/SocketBundleSetter.hh"
#include "LagrangianSolver/ParticleData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace LagrangianSolver {

//////////////////////////////////////////////////////////////////////////////


class ParticleTracking : public Framework::SocketBundleSetter
{
public:
  
  /**
   * Constructor
   */
  ParticleTracking(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ParticleTracking();
  
  enum FaceType {WALL_FACE=1, INTERNAL_FACE = 2, COMP_DOMAIN_FACE=3, BOUNDARY_FACE=4};
  
  virtual void newParticle( CommonData particleCommonData )
  {
    m_particleCommonData = particleCommonData;
    m_particleCommonData.currentPoint = particleCommonData.currentPoint;
    m_particleCommonData.direction = particleCommonData.direction;
  }

  virtual void setupAlgorithm() = 0;

  virtual void getCommonData(CommonData &data)=0;

  virtual void newDirection(RealVector &direction)=0;

  virtual void trackingStep()=0;

  virtual void getExitPoint(RealVector &exitPoint) = 0;

  virtual CFreal getStepDistance() = 0;

  virtual void getNormals(CFuint faceID, RealVector &CartPosition, RealVector &faceNormal) = 0;

  CFuint getExitCellID() const {return m_exitCellID;}

  CFuint getExitFaceID() const {return m_exitFaceID;}

  void setFaceTypes(MathTools::CFMat<CFint> & wallTypes, std::vector<std::string>& wallNames,
                    std::vector<std::string>& boundaryNames);
  
protected: //data
  void getAxiNormals(CFuint faceID, RealVector &CartPosition, RealVector &faceNormal);
  
  void getCartNormals(CFuint faceID, RealVector &CartPosition, RealVector &faceNormal);

  MPI_Datatype m_particleDataType;

  CFuint m_exitFaceID, m_entryFaceID;

  CFuint m_exitCellID, m_entryCellID;

  CFuint m_cellIdx, m_faceIdx;
  
  CommonData m_particleCommonData;

};

} // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RadiativeTransfer_ParticleTracking_hh
