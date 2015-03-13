#include "ParticleTrackingAxi.hh"
#include "Framework/MeshData.hh"

namespace COOLFluiD {

namespace LagrangianSolver {

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

ParticleTrackingAxi::ParticleTrackingAxi(const std::string& name) :
    ParticleTracking(name),
    //m_tCantidates(),
    //m_fCandidates(),
    faceOutNormal(2),
    rayTangent(2)
{
    //m_tCantidates.reserve(5);
    //m_fCandidates.reserve(5);
}


ParticleTrackingAxi::~ParticleTrackingAxi()
{
}

void ParticleTrackingAxi::getCommonData(CommonData &data){
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


void ParticleTrackingAxi::newParticle(CommonData &particle)
{
    static RealVector buffer(3);
    ParticleTracking::newParticle(particle);

    m_particle_t_old=1e-8;
    m_particle_t=1e-8;

    m_entryCellID = m_particleCommonData.cellID;
    m_exitCellID = m_entryCellID;
    m_exitFaceID=-1;

    buffer[0] = m_particleCommonData.direction[0];
    buffer[1] = m_particleCommonData.direction[1];
    buffer[2] = m_particleCommonData.direction[2];

    newDirection(buffer);
}

void ParticleTrackingAxi::setupAlgorithm(){
  m_maxNbFaces = Framework::MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell();
}

//Maybe use the cell->getNeighborGeos() for caching.
//TODO: USE 2D algorithm with last ray tangent for first aproximation
//TODO: Don't need to calculate the intersections on an entry face(same as the last cell's exit face)


void ParticleTrackingAxi::getExitPoint(RealVector &exitPoint){
    cf_assert(exitPoint.size()==3);
    exitPoint[0] = m_particleCommonData.currentPoint[0] + m_particle_t * m_particleCommonData.direction[0];
    exitPoint[1] = m_particleCommonData.currentPoint[1] + m_particle_t * m_particleCommonData.direction[1];
    exitPoint[2] = m_particleCommonData.currentPoint[2] + m_particle_t * m_particleCommonData.direction[2];
}


CFreal ParticleTrackingAxi::getStepDistance()
{
  return (m_particle_t - m_particle_t_old);
}

void ParticleTrackingAxi::newDirection(RealVector &direction){
  static RealVector initialPoint(3);
  cf_assert(direction.size() == 3);
  getExitPoint(initialPoint);

  m_particle_t_old=1e-8;
  m_particle_t=1e-8;

  m_particleCommonData.currentPoint[0]=initialPoint[0];
  m_particleCommonData.currentPoint[1]=initialPoint[1];
  m_particleCommonData.currentPoint[2]=initialPoint[2];

  m_particleCommonData.direction[0]= direction[0];
  m_particleCommonData.direction[1]= direction[1];
  m_particleCommonData.direction[2]= direction[2];

   m_a =    direction[0]; m_b  =    direction[1]; m_c  =    direction[2];
  m_x0 = initialPoint[0]; m_y0 = initialPoint[1]; m_z0 = initialPoint[2];

  //Precalculate direction-dependent values
  m_D22=pow(m_c*m_y0-m_b*m_z0,2);
  m_cx0=m_a*m_z0;
  m_cy0=m_a*m_y0;
  m_a1=pow(m_c,2)+pow(m_b,2);
  m_a2=(m_z0*m_c + m_y0*m_b);
  m_ca2=m_a*m_a2;
  m_c_2=pow(m_a,2);
}

void ParticleTrackingAxi::trackingStep(){
  m_entryCellID = m_exitCellID;

  m_exitCellID=-1;
  m_exitFaceID=-1;

//  cout<<"****************************************************"<<endl;
//  cout<<"New Step!"<<endl<<"on Cell with ID: "<<  m_entryCellID<<endl;

//  cout<<"direction: "<<
//     m_particleCommonData.direction[0]<<' '<<
//     m_particleCommonData.direction[1]<<' '<<
//     m_particleCommonData.direction[2]<< endl;

//  cout<<"Particle t: "<<m_particle_t<<endl;

//  cout<<"position: "<<
//        m_particleCommonData.currentPoint[0] + m_particle_t * m_particleCommonData.direction[0]<<' '<<
//        m_particleCommonData.currentPoint[1] + m_particle_t * m_particleCommonData.direction[1]<<' '<<
//        m_particleCommonData.currentPoint[2] + m_particle_t * m_particleCommonData.direction[2]<< endl;

  static DataHandle<CFint> faceIsOutwards= m_sockets.isOutward.getDataHandle();

  static vector<CFreal> t_candidates(m_maxNbFaces*2);
  static vector<CFuint> f_candidates(m_maxNbFaces*2);

  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  //this->m_cellIdx = this->m_CellIDmap.find(this->m_entryCellID);
  m_cellIdx = m_entryCellID;
  cellData.idx = m_cellIdx;

  GeometricEntity *const cell = m_cellBuilder.buildGE();
//  vector<Node*>& myNodes = *cell->getNodes();
//  for(CFuint ii=0; ii<myNodes.size()-1; ++ii){
//    m_z1 = (*(myNodes[ii]))[XX]; m_z2 = (*(myNodes[ii+1]))[XX];
//    m_r1 = (*(myNodes[ii]))[YY]; m_r2 = (*(myNodes[ii+1]))[YY];
//    cout<<"line: ( "<< m_z1 <<" , "<< m_r1 <<" ) to ( "<< m_z2 <<" , "<< m_r2 <<" )"<<endl;
//  }

  CFuint candId = 0;
  CFuint nFaces = cell->nbNeighborGeos();
  for(CFuint f=0; f<nFaces; ++f){
    GeometricEntity* const face = cell->getNeighborGeo(f);
    vector<Node*>& nodes = *face->getNodes();
    CFint faceID=face->getID();
    m_z1 = (*(nodes[0]))[XX]; m_z2 = (*(nodes[1]))[XX];
    m_r1 = (*(nodes[0]))[YY]; m_r2 = (*(nodes[1]))[YY];
    //cout<<"line: ( "<< m_z1 <<" , "<< m_r1 <<" ) to ( "<< m_z2 <<" , "<< m_r2 <<" )"<<endl;
    //cout<<' ';
    m_dr=m_r2-m_r1; m_dz=m_z2-m_z1; m_zz=m_x0-m_z1;
    m_dr_2=pow(m_dr,2); m_dz_2=pow(m_dz,2);

    faceOutNormal[0] =  m_dr;
    faceOutNormal[1] = -m_dz;
    // this way we assume that the parametrization runs clockwise over the cell's facets.
    // so that the face normal points outwards
    // CF doesn't impose this condition ( bug? ), so we have to correct it
    faceOutNormal *= ( static_cast<CFuint>(faceIsOutwards[faceID])==m_cellIdx) ? 1.:-1.;
    const CFreal temp1 = m_dr*(m_c*m_zz-m_cx0)+m_c*m_dz*m_r1;
    const CFreal temp2 = m_b*m_dz*m_r1+m_dr*(m_b*m_zz-m_cy0);
    m_D=temp1*temp1 + temp2*temp2 - m_dz_2*m_D22;

    if (m_D>=0){
      m_invA=1/(m_dz_2*m_a1-m_c_2*m_dr_2);
      m_D=std::sqrt(m_D)*m_invA;
      m_D1=m_a*m_D;
      m_B1=(m_dz*(m_zz*m_a1-m_ca2)+m_c_2*m_dr*m_r1)*m_invA;

      m_s[0]=m_B1+m_D1; m_s[1]=m_B1-m_D1;

      isS1=m_s[0]<=1. && m_s[0]>=0.; // if the intersection point
      isS2=m_s[1]<=1. && m_s[1]>=0.; // lays between the vertices

      if(isS1|| isS2){
        //RealVector faceNormal(2);

        //cout<<"found something"<<endl;
        m_t[0]=-1.;m_t[1]=-1.;
        m_D2=m_dz*m_D;
        m_B2=(m_a*(m_zz*m_dr_2+m_r1*m_dz*m_dr)-m_dz_2*m_a2)*m_invA;

        CFuint kk=0;
        m_t[kk]=m_B2+m_D2;
        kk= isS1 ? kk+1 : kk;
        m_t[kk]=m_B2-m_D2;
        kk= isS2 ? kk+1 : kk;
        for (CFuint k=0;k<kk;++k){
          rayTangent[0]=m_a;
          rayTangent[1]=( m_c*( m_z0+m_c*m_t[k] ) + m_b*( m_y0+m_b*m_t[k]) )/
                         std::sqrt( pow(m_z0+m_c*m_t[k],2) + pow(m_y0+m_b*m_t[k],2) );
          //avoid use pow!! be cache friendly :)
          //cout<<"fount t= "<<m_t[k]<<endl;
          CFreal innerProd = MathFunctions::innerProd(rayTangent,faceOutNormal);
          //cout<<"InnerPRod: "<<innerProd;
          if (m_t[k]>m_particle_t && innerProd>=0){ //found exit point
             //cout<<"fount t:"<<m_t[k]<<endl;
             t_candidates[candId] = m_t[k];
             f_candidates[candId] = f;
             ++candId;
          }
        }
      }
    }
    //cout<<currentCellID<<" got the nodes"<<endl;
  }
  //cout<<"number of candidates: "<<tCantidates.size()<<endl;

  if(candId==0){
      CFLog(INFO,"Can't find an exit point! \n");
      m_exitCellID=-1;
      m_exitFaceID=-1;
      m_cellBuilder.releaseGE();
      return;
  }


  CFreal myT= t_candidates[0];
  CFuint index=0;
  bool condition;
  for(CFuint i=1;i< candId; ++i){
     // cout<<"looking tcandidate: " << m_tCantidates[i] <<endl;
      condition = ( t_candidates[i]<myT);
      //  cout<<"Choosen!"<<endl;
      myT = condition ? t_candidates[i] : myT;
      index = condition ? i : index;
  }

  //cout<<"Particle t: "<<myT<<endl;
  //cout<<"Particle t old : "<<m_particle_t_old<<endl;

  m_particle_t_old=m_particle_t;
  m_particle_t=myT;

  //cout<<"ray.t: "<<ray.tt<<endl;
  //cout<<"localGE face ID: "<<fCandidates[index]<<endl;
  //cout<<currentCellID<<" END2"<<endl;
  GeometricEntity* exitFace = cell->getNeighborGeo(f_candidates[index]);
  m_exitFaceID = exitFace->getID();

  m_exitFaceID = exitFace->getID();
  m_exitCellID = exitFace->getState(0)->getLocalID();

  m_exitCellID =  (m_exitCellID ==  m_entryCellID && !exitFace->getState(1)->isGhost() ) ?
    exitFace->getState(1)->getLocalID() : m_exitCellID;

  m_cellBuilder.releaseGE();


}

/////////////////////////////////////////////////////////////////////////////
}
}

