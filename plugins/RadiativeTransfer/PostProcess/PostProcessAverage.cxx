#include <mpi.h>
#include <cmath>

#include "FiniteVolume/CellCenterFVM.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntity.hh"
#include "Environment/ObjectProvider.hh"

#include "PostProcessAverage.hh"
#include "RadiativeTransfer/RadiativeTransfer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

Environment::ObjectProvider<PostProcessAverage, PostProcess, RadiativeTransferModule, 1 >
PostProcessAverageProvider("PostProcessAverage");


//////////////////////////////////////////////////////////////////////////////
void PostProcessAverage::defineConfigOptions(Config::OptionList& options){

  options.addConfigOption< vector<CFreal> > ("RegectionDists","Distances from the boundaries prescribed above that the point is neglected");
  options.addConfigOption< vector<CFreal> >("AveragingDirection","Direction of whitch the averaging is done");
  options.addConfigOption< vector<string> >("RegectionTRS","TRS's where the regection distance is used");
  options.addConfigOption< CFuint >("nbIntervals","Number of averaging intervals");
  options.addConfigOption< bool >("mirrorOnOrigin","Mirror points on Origin");
  options.addConfigOption< bool >("outNormals","If the nornals point outwars the mesh");
  options.addConfigOption< bool >("isRevolutionVolume3d","Flag to tell if the 3D geomety has axial symmetry");
}

//////////////////////////////////////////////////////////////////////////////

void PostProcessAverage::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}

PostProcessAverage::PostProcessAverage(const string &name):
  PostProcess(name)
{
  addConfigOptionsTo(this);

  m_regectionDists = vector<CFreal>();
  setParameter("RegectionDists",&m_regectionDists);

  m_regectionTRS=std::vector<std::string>();
  setParameter("RegectionTRS",&m_regectionTRS);

  m_direction = vector<CFreal>();
  setParameter("AveragingDirection",&m_direction);

  m_nbIntervals = 10;
  setParameter("nbIntervals",&m_nbIntervals);

  m_isRevolutionVolume3d = false;
  setParameter("isRevolutionVolume3d",&m_isRevolutionVolume3d);

  m_mirrorOnOrigin = false;
  setParameter("mirrorOnOrigin",&m_mirrorOnOrigin);

  m_outNormals = true;
  setParameter("outNormals",&m_outNormals);
}

PostProcessAverage::~PostProcessAverage()
{
}

void PostProcessAverage::runPostProcess(DataHandle<CFreal> dataVector){

  FaceTrsGeoBuilder::GeoData& facesData = m_faceBuilder.getDataGE();
  DataHandle<CFreal> normals = m_sockets.normals.getDataHandle();
  CFuint dim = PhysicalModelStack::getActive()->getDim();

  cf_always_assert( m_regectionDists.size() == m_regectionTRS.size() );
  cf_always_assert( dim == m_direction.size() );

  //normalize vector
  CFreal invNorm= 0.;
  for (CFuint i=0; i< m_direction.size(); ++i){
    invNorm+=m_direction[i]*m_direction[i];
  }
  invNorm = 1/std::sqrt(invNorm);
  for (CFuint i=0; i< m_direction.size(); ++i){
    m_direction[i]*=invNorm;
  }

  vector<RealVector> faceNormal(0, dim);
  vector<RealVector> initialPoint(0, dim);
  vector<RealVector> maxPoint(0, dim);
  RealVector bufferFaceNormal(dim), bufferInitialPoint(dim), bufferMaxPoint(dim);

  m_faceBuilder.getDataGE().isBFace = true;
  for(CFuint j=0; j<m_regectionTRS.size(); ++j){

    SafePtr<TopologicalRegionSet> WallFaces = MeshDataStack::getActive()->getTrs( m_regectionTRS[j] );

    facesData.trs = WallFaces;
    const CFuint nbFacesWall = WallFaces->getLocalNbGeoEnts();

    if ( nbFacesWall > 0 ){

      facesData.trs = WallFaces;
      facesData.idx=0;

      //Get the normals (assumes all faces of the TRS have the same normal)
      GeometricEntity *const face = m_faceBuilder.buildGE();
      CFuint startID=face->getID()*dim;
      for(CFuint i=0;i< dim ;++i){
        bufferFaceNormal[i]= -normals[startID+i];
      }


      bufferFaceNormal.normalize();
      bufferInitialPoint= face->computeCentroid();

      bufferMaxPoint = bufferInitialPoint + bufferFaceNormal * m_regectionDists[j];

      faceNormal.push_back(bufferFaceNormal);
      initialPoint.push_back(bufferInitialPoint);
      maxPoint.push_back(bufferMaxPoint);

      m_faceBuilder.releaseGE();
    }

  }
  //Define the plane of regection for the boundaries

  vector<CFreal> cellPos;
  vector<CFreal> Q_out;

  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();


  cellPos.reserve(nCells);
  Q_out.reserve(nCells);

  for(CFuint c=0; c<nCells; ++c){
    cellData.idx = c;
    GeometricEntity* const cell = m_cellBuilder.buildGE();
    static Framework::DataHandle<Framework::State*, Framework::GLOBAL> states
        = m_sockets.states.getDataHandle();

    RealVector centroid = (*states[c]).getCoordinates();

    if (cell->getState(0)->isParUpdatable()){

      //Check is centroid is inside the regection area
      bool isRegected = false;

      for(CFuint i=0; i<maxPoint.size();++i){
        RealVector v1 = centroid - initialPoint[i];
        RealVector v2 =  faceNormal[i] * ( m_outNormals? -1. : 1. );

        isRegected |= MathFunctions::innerProd( v1,v2  ) < m_regectionDists[i] ;
      }

      if(!isRegected){
        CFreal distCell;
          if(m_isRevolutionVolume3d){
              distCell =  std::sqrt(centroid[0]*centroid[0] + centroid[1]*centroid[1]);
          }
          else{
             distCell = MathFunctions::innerProd(m_direction , centroid);
          }

        cellPos.push_back( distCell  );
        Q_out.push_back( dataVector[c] );

        if(m_mirrorOnOrigin){
          Q_out.push_back( dataVector[c] );
          cellPos.push_back( -distCell  );
        }

      }
    }
    m_cellBuilder.releaseGE();
  }

  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  
  CFuint nbProcesses = PE::GetPE().GetProcessorCount(nsp);
  MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
  CFuint myProcessRank = PE::GetPE().GetRank(nsp);
  
  //send all tuples to process 0 to be stored in a common file
  vector<int> nbsElems, disps;
  
  nbsElems.resize(nbProcesses,0);
  disps.resize(nbProcesses,0);
  
  int local_size = Q_out.size();
  MPI_Gather(&local_size,  1, Common::MPIStructDef::getMPIType(&local_size),
             &nbsElems[0], 1, Common::MPIStructDef::getMPIType(&nbsElems[0]),
             0, comm);
  
  CFuint globalSize=0;
  for (CFuint i=0; i< nbProcesses; ++i){
    disps[i] = globalSize;
    globalSize += nbsElems[i];
  }
  
  vector<CFreal> g_cellPos, g_Qout;
  
  g_cellPos.resize(globalSize,0);
  g_Qout.resize(globalSize,0);
  
  MPI_Gatherv(&Q_out[0],  Q_out.size(), MPI_DOUBLE,
              &g_Qout[0], &nbsElems[0], &disps[0],  MPI_DOUBLE, 0, comm);
  
  MPI_Gatherv(&cellPos[0],   cellPos.size(), MPI_DOUBLE,
              &g_cellPos[0], &nbsElems[0],   &disps[0],  MPI_DOUBLE, 0, comm);
  
  ofstream debug;
  if(myProcessRank == 0 ){
    debug.open("postprocessDebug.dat");

    for(CFuint i=0; i< g_cellPos.size();++i){
      debug<<g_cellPos[i]<<' '<<g_Qout[i]<<endl;
    }
    debug.close();
  }

  // calculate averages in parallel

  CFreal localPosMax= *max_element(cellPos.begin(),cellPos.end());
  CFreal localPosMin= *min_element(cellPos.begin(),cellPos.end());
  CFreal globalPosMax, globalPosMin;

  MPI_Allreduce(&localPosMax, &globalPosMax, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&localPosMin, &globalPosMin, 1, MPI_DOUBLE, MPI_MIN, comm);

  const CFreal deltaPos = (globalPosMax-globalPosMin);
  globalPosMax += deltaPos * 1e-4;
  globalPosMin -= deltaPos * 1e-4;

  CFreal nbSteps = m_nbIntervals +1 ;

  std::map<CFreal,CFuint> boundaries;
  vector<CFreal> qq(m_nbIntervals,0), global_qq(m_nbIntervals,0);
  vector<CFreal> nb_qq(m_nbIntervals,0), global_nb_qq(m_nbIntervals,0);
  vector<CFreal> rr(nbSteps,0);

  CFreal step=(globalPosMax-globalPosMin)/CFreal(m_nbIntervals);
  for (CFuint i=0; i<nbSteps; ++i){
    rr[i] = globalPosMin+step*i;
    boundaries.insert( boundaries.end(), pair<CFreal,CFuint>( rr[i], i) );
  }

  std::map<CFreal,CFuint>::iterator it;
  for(CFuint i=0; i<Q_out.size(); ++i){
    it = boundaries.upper_bound(cellPos[i]);
    it--;
    qq[ it->second ] += Q_out[i];
    nb_qq[it->second] += 1.;
  }


  MPI_Reduce(&qq[0]   , &global_qq[0]   , m_nbIntervals, MPI_DOUBLE, MPI_SUM, 0 ,comm);
  MPI_Reduce(&nb_qq[0], &global_nb_qq[0], m_nbIntervals, MPI_DOUBLE, MPI_SUM, 0 ,comm);

  if(myProcessRank == 0 ){

    debug.open("postprocess.plt");

    for(CFuint i=0; i< global_qq.size();++i){
        if( global_nb_qq[i] > 1e-1 ) {
          debug<< 0.5*(rr[i]+rr[i+1])<<' '<<global_qq[i]/global_nb_qq[i]<<endl;
        }
    }
    debug.close();
  }

}

}

}
