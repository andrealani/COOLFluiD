#include "ParticleTracking3D.hh"
#include <iostream>
#include <cmath>
#include <algorithm>
#include "Framework/MeshData.hh"

namespace COOLFluiD {

namespace LagrangianSolver {

ParticleTracking3D::ParticleTracking3D(const std::string& name):
    ParticleTracking(name),
    m_exitPoint(3),
    m_entryPoint(3),
    m_direction(3)
{

}

ParticleTracking3D::~ParticleTracking3D(){

}

void ParticleTracking3D::getCommonData(CommonData &data){
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



void ParticleTracking3D::setupAlgorithm(){
    m_maxNbFaces = Framework::MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell();

    CFuint nbFaces=Framework::MeshDataStack::getActive()->Statistics().getNbFaces();
    m_centroids.reserve(nbFaces*3);

    using namespace Framework;
    using namespace Common;

    CellTrsGeoBuilder::GeoData& cellsData = m_cellBuilder.getDataGE();
    SafePtr<TopologicalRegionSet> MediumCells = MeshDataStack::getActive()->getTrs("InnerCells");
    cellsData.trs = MediumCells;
    const CFuint nbCellsMedia = MediumCells->getLocalNbGeoEnts();
    for(CFuint i=0; i<nbCellsMedia; ++i){
      cellsData.idx = i;
      GeometricEntity *const cell = m_cellBuilder.buildGE();
      CFuint nFaces = cell->nbNeighborGeos();
      for(CFuint f=0; f<nFaces; ++f){
        GeometricEntity* const face = cell->getNeighborGeo(f);
        CFuint faceID = face->getID();
        const RealVector centroid = face->computeCentroid();
        m_centroids[faceID*3+0] = centroid[0];
        m_centroids[faceID*3+1] = centroid[1];
        m_centroids[faceID*3+2] = centroid[2];
      }
      m_cellBuilder.releaseGE();
    }
}


void ParticleTracking3D::trackingStep(){

  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;


 // cout<<"%\n%Start tracking; CellID: "<<m_exitCellID<<endl<<"%\n";
  m_exitFaceID=-1;
  //m_exitCellID=-1;

  const RealVector &initialPoint = m_exitPoint;

  static DataHandle<CFint> faceIsOutwards= m_sockets.isOutward.getDataHandle();
  static DataHandle<CFreal> normals= m_sockets.normals.getDataHandle();
  //static DataHandle<CFreal> faceCenters= m_sockets.faceCenters.getDataHandle();
  vector<CFreal> &faceCenters = m_centroids;

  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  m_entryCellID = m_exitCellID;
  cellData.idx = m_entryCellID;

  GeometricEntity *const cell = m_cellBuilder.buildGE();
  CFuint nFaces = cell->nbNeighborGeos();

//  cout<<"direction = ["<<m_direction<<"];"<<endl;
//  cout<<"initialPoint = ["<< initialPoint<<"];"<<endl;

  static CFreal temp2[3];
  static CFreal temp3[3];
  static CFreal faceOutNormal[3];
  static vector<CFreal> innerProds( m_maxNbFaces );

  for(CFuint f=0; f<nFaces; ++f){

    GeometricEntity* const face = cell->getNeighborGeo(f);
    CFint faceID=face->getID();
    vector<Node*>& myNodes = *face->getNodes();
    //RealVector centroid = face->computeCentroid();


    const CFuint nbNodes= myNodes.size();

    for(CFuint ii=0; ii<nbNodes; ++ii){
      CFuint ii_1 =( ii == nbNodes-1 ) ? 0 : ii +1;

      temp2[0] = (*myNodes[ii  ])[0] - initialPoint[0];
      temp2[1] = (*myNodes[ii  ])[1] - initialPoint[1];
      temp2[2] = (*myNodes[ii  ])[2] - initialPoint[2];

      temp3[0] = (*myNodes[ii_1])[0] - initialPoint[0];
      temp3[1] = (*myNodes[ii_1])[1] - initialPoint[1];
      temp3[2] = (*myNodes[ii_1])[2] - initialPoint[2];

      //MathFunctions::crossProd(temp2  , temp3, temp1);

      innerProds[ii] = ( temp2[1]*temp3[2] - temp2[2]*temp3[1] ) * m_direction[0] +
                       (-temp2[0]*temp3[2] + temp2[2]*temp3[0] ) * m_direction[1] +
                       ( temp2[0]*temp3[1] - temp2[1]*temp3[0] ) * m_direction[2];

      //innerProds[ii] = MathFunctions::innerProd(temp1, m_direction );


//     cout<<"pointsFace("<<ii+1<<",:,"<< f+1 <<") = [" << *myNodes[ii] <<"];"<<endl;
//     cout<<"pointsFace("<<ii+2<<",:,"<< f+1 <<") = [" << *myNodes[ii_1] <<"];"<<endl;
    }

    const CFreal circulation = (static_cast<CFuint>(faceIsOutwards[faceID])==m_entryCellID) ? 1.:-1.;
    CFuint faceIdx = faceID*3;
    faceOutNormal[0]=normals[faceIdx+0]*circulation;
    faceOutNormal[1]=normals[faceIdx+1]*circulation;
    faceOutNormal[2]=normals[faceIdx+2]*circulation;

    //centroid[0] = faceCenters[faceIdx+0];
    //centroid[1] = faceCenters[faceIdx+1];
    //centroid[2] = faceCenters[faceIdx+2];

//    cout<<"faceCenter("<< f+1<<",:) = ["<<centroid<<"];" <<endl;


//    RealVector faceOutNormal(3);
//    faceOutNormal.normalize();
//    faceOutNormal*=circulation;
    //bool isExitFaceCondition = MathFunctions::innerProd(faceOutNormal, m_direction) > 0;

    const CFreal dot_direction_normal = faceOutNormal[0] * m_direction[0] +
                                        faceOutNormal[1] * m_direction[1] +
                                        faceOutNormal[2] * m_direction[2] ;

    bool isExitFaceCondition = (dot_direction_normal> 0. );

//    cout<<"circulation("<< f+1<<",:) ="<<circulation<<";" <<endl;
//    cout<<"normal("<< f+1<<",:) = ["<< faceOutNormal[0]<<' '
//                                    << faceOutNormal[1]<<' '
//                                    << faceOutNormal[2] <<"];"<<endl;

//    cout<<"%innerProds: ";
//    for(CFuint ii=0; ii<innerProds.size(); ++ii){
//        cout<<innerProds[ii]*circulation<<' ';
//    }
//    cout<<endl;

    bool crossesFaceCondition = true;
    for(CFuint ii=0; ii<nbNodes-1; ++ii){
      crossesFaceCondition &= (innerProds[ii]*innerProds[ii+1]  > 0);
    }

    //cout<<"%Crosses face: "<<crossesFaceCondition<<" isExitFaceCondition: "<<isExitFaceCondition<<endl<<endl;
    if (isExitFaceCondition && crossesFaceCondition){
      //cout<<"%IS EXIT FACE"<<endl;


        //const RealVector temp4 =   centroid - initialPoint;

        //m_stepDist = MathFunctions::innerProd(temp4, faceOutNormal) /
        //    dot_direction_normal;

        //m_stepDist = dot(centroid - initialPoint, normal) / dot(direction,normal);


      m_stepDist =( (faceCenters[faceIdx+0] - initialPoint[0]) * faceOutNormal[0] +
                    (faceCenters[faceIdx+1] - initialPoint[1]) * faceOutNormal[1] +
                    (faceCenters[faceIdx+2] - initialPoint[2]) * faceOutNormal[2]
                  ) / ( dot_direction_normal );


      //cout<<"%stepDist = "<< m_stepDist <<endl;
      m_exitPoint = initialPoint + m_direction * m_stepDist;
      //cout<<"my_exitPoint = [ "<<m_exitPoint<<" ];"<<endl<<endl;


      m_exitFaceID = face->getID();
      m_exitCellID = face->getState(0)->getLocalID();

      m_exitCellID =  (m_exitCellID ==  m_entryCellID && !face->getState(1)->isGhost() ) ?
        face->getState(1)->getLocalID() : m_exitCellID;

      break;
    }
  }

  m_cellBuilder.releaseGE();
}

void ParticleTracking3D::newParticle(CommonData &particle){

    //std::cout<<"%*******************\n %NEW PARTICLE \n %************************************\n";

  static RealVector buffer(3);
  ParticleTracking::newParticle(particle);

  m_entryCellID = m_particleCommonData.cellID;
  m_exitCellID = m_entryCellID;
  m_exitFaceID=-1;

  buffer[0] = m_particleCommonData.direction[0];
  buffer[1] = m_particleCommonData.direction[1];
  buffer[2] = m_particleCommonData.direction[2];

  m_exitPoint[0] = m_particleCommonData.currentPoint[0];
  m_exitPoint[1] = m_particleCommonData.currentPoint[1];
  m_exitPoint[2] = m_particleCommonData.currentPoint[2];

  newDirection(buffer);
}


}
}
