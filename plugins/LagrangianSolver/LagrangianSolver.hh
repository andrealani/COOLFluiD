#ifndef COOLFluiD_RadiativeTransfer_Lagrangian_Solver_hh
#define COOLFluiD_RadiativeTransfer_Lagrangian_Solver_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/MPI/MPIStructDef.hh"
#include "LagrangianSolverModule.hh"
#include "ParticleTracking/ParticleTracking.hh"
#include "Framework/SocketBundleSetter.hh"
#include "ParticleTracking/ParticleTrackingAxi.hh"
#include "ParticleTracking/ParticleTracking3D.hh"
#include "SendBuffer/SendBuffer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace LagrangianSolver {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class implements the algorithms for the computation of local parameters
   *
   * @author Pedro Santos
   * @author Alessandro Sanna
   *
   */

class ParticleTracking;

template<typename UserData, class PARTICLE_TRACKING>
class LagrangianSolver
{
public:

  /**
   * Constructor
   */
  LagrangianSolver(const std::string& name);

  /**
   * Default destructor
   */
   ~LagrangianSolver();

  inline void newParticle(Particle<UserData> particle)
  {
    m_particle = particle;
    m_particleTracking.newParticle(m_particle.commonData);
  }
  
  inline void setupSendBufferSize(CFuint sendBufferSize)
  {
    try {
      m_sendBuffer.reset(new SendBuffer<Particle<UserData> >());
      m_sendBufferSize = sendBufferSize;
      m_sendBuffer->reserve(m_sendBufferSize);
    } 
    catch (std::bad_alloc& ba) {
      std::cerr << "bad_alloc caught: " << ba.what() << '\n';
    }
  }
  
  inline void newDirection(RealVector direction){ m_particleTracking.newDirection(direction);}
  
  inline void trackingStep(){ m_particleTracking.trackingStep(); }

   inline void getExitPoint(RealVector &exitPoint){ m_particleTracking.getExitPoint(exitPoint); }

   inline CFreal getStepDistance(){ return m_particleTracking.getStepDistance(); }

   inline void getUserData(UserData &data){ data = m_particle.userData; }

   inline void setUserData(UserData &data){ m_particle.userData=data; }

   inline UserData& getUserDataPtr(){return (m_particle.userData); }

   inline void getParticle(Particle<UserData> &particle){ getUserData(particle.userData);
                                                          m_particleTracking.getCommonData(particle.commonData); }

   inline void getCommonData(CommonData &commonData){m_particleTracking.getCommonData(commonData);}

   inline CFuint getExitCellID(){ return m_particleTracking.getExitCellID(); }

   inline CFuint getExitFaceID(){ return m_particleTracking.getExitFaceID(); }

   void setupParallelization();

   inline bool sincronizeParticles(std::vector< Particle<UserData> >&particleBuffer, bool isLastPhoton);

   void setupParticleDatatype(MPI_Datatype ptrDatatype );

   inline MPI_Datatype getParticleDataType(){return m_particleDataType; }

  inline void setDataSockets( Framework::SocketBundle sockets )
  { 
    m_particleTracking.setDataSockets(sockets);
    m_particleTracking.setupAlgorithm();
  }

   inline void getNormals(CFuint faceID, RealVector& CartPosition, RealVector& faceNormal ){
       m_particleTracking.getNormals( faceID, CartPosition, faceNormal );
   }

  inline void setFaceTypes(std::vector<std::string>& wallNames, std::vector<std::string>& boundaryNames){
       m_particleTracking.setFaceTypes(m_wallTypes, wallNames, boundaryNames );
   }

   //inline CFuint getFaceStateID(CFuint faceID){return m_wallTypes(faceID,1);}

   inline CFuint getWallGhotsStateId(CFuint faceID){
     cf_assert(m_wallTypes(faceID,0) !=  ParticleTracking::INTERNAL_FACE );
     return m_wallTypes(faceID,2);
   }

   inline CFuint getWallStateId(CFuint faceID){
     cf_assert(m_wallTypes(faceID,0) !=  ParticleTracking::INTERNAL_FACE );
     return m_wallTypes(faceID,1);
   }

   inline CFuint getFaceType(CFuint faceID){
       cf_assert(m_wallTypes(faceID,0)!= -1 );
       return m_wallTypes(faceID,0);
   }

   CFuint getTargetProcess(CFuint faceID) const;

   void bufferCommitParticle(CFuint faceID);

private:

  void (ParticleTracking::*getNormalsPtr) (CFuint, RealVector, RealVector);

  MPI_Datatype m_particleDataType;

  Particle<UserData> m_particle;
  
  PARTICLE_TRACKING m_particleTracking;
  
  MathTools::CFMat<CFint> m_wallTypes;
  
  std::auto_ptr<SendBuffer<Particle<UserData> > > m_sendBuffer;
  
  CFuint m_sendBufferSize;

};

//////////////////////////////////////////////////////////////////////////////

template<typename UserData, class PARTICLE_TRACKING>
LagrangianSolver<UserData,PARTICLE_TRACKING>::LagrangianSolver(const std::string& name) :
  m_particleTracking(name)
{
}

//////////////////////////////////////////////////////////////////////////////

template<typename UserData, class PARTICLE_TRACKING>
LagrangianSolver<UserData,PARTICLE_TRACKING>::~LagrangianSolver()
{
}

//////////////////////////////////////////////////////////////////////////////

template<typename UserData, class PARTICLE_TRACKING>
CFuint LagrangianSolver<UserData, PARTICLE_TRACKING>::getTargetProcess(CFuint faceID) const
{
  return m_wallTypes(faceID,2);
}


//////////////////////////////////////////////////////////////////////////////
template<typename UserData, class PARTICLE_TRACKING>
void LagrangianSolver<UserData, PARTICLE_TRACKING>::bufferCommitParticle(CFuint faceID)
{
  const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
  
  static CFint rank = Common::PE::GetPE().GetRank(nsp);
  cf_assert(m_wallTypes(faceID,0) == ParticleTracking::COMP_DOMAIN_FACE );
  cf_assert(m_wallTypes(faceID,2) != rank );
  
  static Particle<UserData> sendParticle;
  getParticle(sendParticle);
  CFuint processRank = m_wallTypes(faceID,2);
  //CFLog(INFO, "processRank: "<<processRank<<"\n");
  sendParticle.commonData.cellID = m_wallTypes(faceID,3);
  m_sendBuffer->push_back(sendParticle, processRank );
}

//////////////////////////////////////////////////////////////////////////////

template<typename UserData, class PARTICLE_TRACKING>
bool LagrangianSolver<UserData,PARTICLE_TRACKING>::sincronizeParticles(std::vector< Particle<UserData> >&particleBuffer,
								       bool isLastPhoton)
{
  return m_sendBuffer->sincronize(particleBuffer, isLastPhoton);
}
    
//////////////////////////////////////////////////////////////////////////////

template<typename UserData, class PARTICLE_TRACKING>
void LagrangianSolver<UserData,PARTICLE_TRACKING>::setupParticleDatatype(MPI_Datatype ptrUserDatatype ){
  
    // *********************************************************
    //@TODO
    //extend MPIStructDef to be able to handle nested datatypes
    //
    // *********************************************************
  using namespace COOLFluiD::Common;
  
    // set the common data datatype
    CommonData commonData;
    MPIStruct commonDataType;

    int counts[3]={3,3,1};
    MPIStructDef::buildMPIStruct<CFreal,CFreal,CFuint>
      (&commonData.direction[0], &commonData.currentPoint[0], &commonData.cellID, &counts[0], commonDataType);
    
    //set the complete particle dataType
    
    Particle<UserData> particle;
    
    MPI_Aint displacements[2], start_address, address;
    
    // get block lenghts
    int block_lengths[2] = { 1,1 };
    MPI_Datatype typeList[2] = {commonDataType.type, ptrUserDatatype};
    
    displacements[0] =0;
    MPI_Get_address(&particle.commonData, &start_address); //Vatsalya: Changed from MPI_Address to make it compatible with MPI 4.0

    MPI_Get_address(&particle.userData, &address); //Vatsalya: Changed from MPI_Address to make it compatible with MPI 4.0
    displacements[1] = address - start_address;

    //create derived datatype
    MPI_Type_create_struct(2, block_lengths, displacements, typeList, &m_particleDataType); //Vatsalya: Changed from MPI_Type_struct to make it compatible with MPI 4.0

    // commit new datatype
    MPI_Type_commit( &m_particleDataType );
    
    m_sendBuffer->setMPIdatatype(m_particleDataType);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RadiativeTransfer_ComputeNonLocalParameters_hh
