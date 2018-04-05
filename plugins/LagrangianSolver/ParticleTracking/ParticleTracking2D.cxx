#include "ParticleTracking2D.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace LagrangianSolver {

//////////////////////////////////////////////////////////////////////////////

ParticleTracking2D::ParticleTracking2D(const std::string& name) :
    ParticleTracking(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ParticleTracking2D::~ParticleTracking2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParticleTracking2D::getCommonData(CommonData &data)
{
  static RealVector initialPoint(3);
  getExitPoint(initialPoint);
  
  data.currentPoint[0]=initialPoint[0];
  data.currentPoint[1]=initialPoint[1];
  
  data.direction[0] = m_particleCommonData.direction[0];
  data.direction[1] = m_particleCommonData.direction[1];
  
  data.cellID = m_particleCommonData.cellID;
}

//////////////////////////////////////////////////////////////////////////////

void ParticleTracking2D::newParticle(CommonData &particle)
{
  static RealVector buffer(2);
  ParticleTracking::newParticle(particle);
    
  m_particle_t_old=1e-8;
  m_particle_t=1e-8;
  
  m_entryCellID = m_particleCommonData.cellID;
  m_exitCellID = m_entryCellID;
  m_exitFaceID=-1;
  
  buffer[0] = m_particleCommonData.direction[0];
  buffer[1] = m_particleCommonData.direction[1];
  
  newDirection(buffer);
}
  
//////////////////////////////////////////////////////////////////////////////

void ParticleTracking2D::setupAlgorithm()
{
  ParticleTracking::setupAlgorithm();
}
  
//////////////////////////////////////////////////////////////////////////////

void ParticleTracking2D::getExitPoint(RealVector &exitPoint)
{
  cf_assert(exitPoint.size()<=3);
  exitPoint[0] = m_particleCommonData.currentPoint[0] + m_particle_t * m_particleCommonData.direction[0];
  exitPoint[1] = m_particleCommonData.currentPoint[1] + m_particle_t * m_particleCommonData.direction[1];
}
  
//////////////////////////////////////////////////////////////////////////////

CFreal ParticleTracking2D::getStepDistance()
{
  return (m_particle_t - m_particle_t_old);
}

//////////////////////////////////////////////////////////////////////////////

void ParticleTracking2D::newDirection(RealVector &direction)
{
  static RealVector initialPoint(2);
  cf_assert(direction.size() <= 3);
  getExitPoint(initialPoint);
  
  m_particle_t_old=1e-8;
  m_particle_t=1e-8;
  
  m_particleCommonData.currentPoint[0]=initialPoint[0];
  m_particleCommonData.currentPoint[1]=initialPoint[1];
  
  m_particleCommonData.direction[0]= direction[0];
  m_particleCommonData.direction[1]= direction[1];
  
  m_a =    direction[0]; m_b  =    direction[1];
  m_x0 = initialPoint[0]; m_y0 = initialPoint[1];
}
  
//////////////////////////////////////////////////////////////////////////////

void ParticleTracking2D::trackingStep()
{
  static DataHandle<CFint> faceIsOutwards= m_sockets.isOutward.getDataHandle();
  m_entryCellID = m_exitCellID;
  
  m_exitCellID=-1;
  m_exitFaceID=-1;
  
  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  
  cellData.idx = m_entryCellID;
  
  GeometricEntity *const cell = m_cellBuilder.buildGE();
  CFuint nFaces = cell->nbNeighborGeos();
  
  //        cout<<"****************************************************"<<endl;
  //        cout<<"New Step!"<<endl<<"on Cell with ID: "<<  m_entryCellID<<endl;
  
  //      cout<<"direction: "<<
  //          m_particleCommonData.direction[0]<<' '<<
  //          m_particleCommonData.direction[1]<<' '<<
  //          m_particleCommonData.direction[2]<< endl;
  
  //      cout<<"Particle t: "<<m_particle_t<<endl;
  
  //      cout<<"position: "<<
  //            m_particleCommonData.currentPoint[0] + m_particle_t * m_particleCommonData.direction[0]<<' '<<
  //            m_particleCommonData.currentPoint[1] + m_particle_t * m_particleCommonData.direction[1]<<' '<<
  //            m_particleCommonData.currentPoint[2] + m_particle_t * m_particleCommonData.direction[2]<< endl;
  
  static DataHandle<CFreal> normals= m_sockets.normals.getDataHandle();
  
  for(CFuint f=0; f<nFaces; ++f){
    GeometricEntity* const exitFace = cell->getNeighborGeo(f);
    vector<Node*>& nodes = *exitFace->getNodes();
    CFint faceID=exitFace->getID();
    m_x1 = (*(nodes[0]))[XX]; m_x2 = (*(nodes[1]))[XX];
    m_y1 = (*(nodes[0]))[YY]; m_y2 = (*(nodes[1]))[YY];
    
    m_dx=m_x2-m_x1; m_dy=m_y2-m_y1;
    m_xx=m_x0-m_x1; m_yy=m_y0-m_y1;
    
    m_A = 1./(m_a*m_dy-m_b*m_dx);
    m_ss =(m_a*m_yy-m_b*m_xx)*m_A;
    
    const CFreal circulation = (static_cast<CFuint>(faceIsOutwards[faceID])==m_entryCellID) ? -1.:1.;
    CFuint faceIdx = faceID*2;
    
    m_innerProd = normals[faceIdx+0]*circulation * m_a + normals[faceIdx+1]*circulation* m_b;
    
    // this way we assume that the parametrization runs clockwise over the cell's facets.
    // so that the face normal points outwards
    // CF doesn't impose this condition ( bug? ), so we have to correct it
    
    const bool isExitPoint = m_ss<1. && m_ss>0. && m_innerProd > 0;
    if(isExitPoint){
      m_particle_t_old=m_particle_t;
      m_particle_t  = (m_dx*m_yy-m_dy*m_xx)*m_A;
      
      m_exitFaceID = exitFace->getID();
      m_exitCellID = exitFace->getState(0)->getLocalID();
      
      m_exitCellID =  (m_exitCellID ==  m_entryCellID && !exitFace->getState(1)->isGhost() ) ?
	exitFace->getState(1)->getLocalID() : m_exitCellID;
      
      m_cellBuilder.releaseGE();
      //          cout<<"found Exit Point:"<<" faceID: "<<
      //                m_exitFaceID<<" newCell ID: "<<m_exitCellID<<" new t: "<<m_particle_t<<endl;
      return;
    }
  }
  
  CFLog(VERBOSE, "ParticleTracking2D::trackingStep() => Can't find an exit Point!!\n");
  m_exitCellID=-1;
  m_exitFaceID=-1;
  m_cellBuilder.releaseGE();
  //    //return out;
}
  
//////////////////////////////////////////////////////////////////////////////

}
}
