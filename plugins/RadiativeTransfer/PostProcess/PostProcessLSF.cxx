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

#include "PostProcessLSF.hh"
#include "RadiativeTransfer/RadiativeTransferModule.hh"

#include <gsl/gsl_multifit.h>


using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

Environment::ObjectProvider<PostProcessLSF, PostProcess, RadiativeTransferModule, 1 >
PostProcessLSFProvider("PostProcessLSF");

//////////////////////////////////////////////////////////////////////////////

void PostProcessLSF::defineConfigOptions(Config::OptionList& options){

    options.addConfigOption< vector<CFreal> > ("RegectionDists","Distances from the boundaries prescribed above that the point is neglected");
    options.addConfigOption< vector<CFreal> >("AveragingDirection","Direction of whitch the averaging is done");
    options.addConfigOption< vector<string> >("RegectionTRS","TRS's where the regection distance is used");
    options.addConfigOption< CFuint >("nbCoeffs","Number of coefficients in the interpolant polynomial");
    options.addConfigOption< CFuint >("nbPointsPlot","Number of Points in the Plot");
    options.addConfigOption< bool >("isRevolutionVolume3d","Flag to tell if the 3D geomety has axial symmetry");
    options.addConfigOption< bool >("mirrorOnOrigin","Mirror points on Origin");
    options.addConfigOption< bool >("outNormals","If the nornals point outwars the mesh");
    options.addConfigOption< CFreal >("ScalingFactor","Scaling Factor");
}

//////////////////////////////////////////////////////////////////////////////

void PostProcessLSF::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}

PostProcessLSF::PostProcessLSF(const string &name):
  PostProcess(name)
{
  addConfigOptionsTo(this);

  m_regectionDists = vector<CFreal>();
  setParameter("RegectionDists",&m_regectionDists);

  m_regectionTRS=std::vector<std::string>();
  setParameter("RegectionTRS",&m_regectionTRS);

  m_direction = vector<CFreal>();
  setParameter("AveragingDirection",&m_direction);

  m_nbCoeffs= 3;
  setParameter("nbCoeffs",&m_nbCoeffs);

  m_nbPointsPlot = 20;
  setParameter("nbPointsPlot",&m_nbPointsPlot);

  m_isRevolutionVolume3d = false;
  setParameter("isRevolutionVolume3d",&m_isRevolutionVolume3d);

  m_mirrorOnOrigin = false;
  setParameter("mirrorOnOrigin",&m_mirrorOnOrigin);

  m_outNormals = true;
  setParameter("outNormals",&m_outNormals);

  m_scalingFactor = 1.;
  setParameter("ScalingFactor",&m_scalingFactor);
}

PostProcessLSF::~PostProcessLSF()
{
}

void PostProcessLSF::runPostProcess(DataHandle<CFreal> dataVector){

  FaceTrsGeoBuilder::GeoData& facesData = m_faceBuilder.getDataGE();
  DataHandle<CFreal> normals = m_sockets.normals.getDataHandle();
  CFuint dim = PhysicalModelStack::getActive()->getDim();

  cf_always_assert( m_regectionDists.size() == m_regectionTRS.size() );
  cf_always_assert( dim == m_direction.size() );

  if (m_isRevolutionVolume3d){
    cf_always_assert( m_isRevolutionVolume3d && dim == 3 );
  }

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
        bufferFaceNormal[i]= normals[startID+i];
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
  //Define the planes of regection for the boundaries

  vector<CFreal> cellPos;
  vector<CFreal> Q_out;

  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();

  cellPos.reserve(nCells*2);
  Q_out.reserve(nCells*2);

  for(CFuint c=0; c<nCells; ++c){
    cellData.idx = c;
    GeometricEntity* const cell = this->m_cellBuilder.buildGE();
    static Framework::DataHandle<Framework::State*, Framework::GLOBAL> states
        = this->m_sockets.states.getDataHandle();

    RealVector centroid = (*states[c]).getCoordinates();

    if (cell->getState(0)->isParUpdatable()){

      //Check is centroid is inside the regection area
      bool isRegected = false;

      for(CFuint i=0; i<maxPoint.size();++i){
        RealVector v1 = centroid - initialPoint[i];
        RealVector v2 =  faceNormal[i] * ( m_outNormals? -1. : 1. );

        isRegected |= MathFunctions::innerProd( v1,v2  ) < m_regectionDists[i] ;
      }
      CFreal distCell;
      if(!isRegected){

          if(m_isRevolutionVolume3d){
              distCell =  std::sqrt(centroid[0]*centroid[0] + centroid[1]*centroid[1]);
          }
          else{
             distCell = MathFunctions::innerProd(m_direction , centroid);
          }

        cellPos.push_back( distCell  );
        Q_out.push_back( dataVector[c] / m_scalingFactor );

        if(m_mirrorOnOrigin){
          Q_out.push_back( dataVector[c] );
          cellPos.push_back( -distCell  );
        }

      }
    }
    this->m_cellBuilder.releaseGE();
  }

  CFuint nbProcesses = PE::GetPE().GetProcessorCount();
  MPI_Comm comm = PE::GetPE().GetCommunicator();
  CFuint myProcessRank = PE::GetPE().GetRank();


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

  CFreal localPosMax;
  CFreal localPosMin;

  CFreal tempAbsPos;
  for(CFuint i=0; i < cellPos.size(); ++i ){
      tempAbsPos = std::abs(cellPos[i]);
      localPosMin = (tempAbsPos < localPosMin) ? tempAbsPos : localPosMin;
      localPosMax = (tempAbsPos > localPosMax) ? tempAbsPos : localPosMax;
  }

  CFreal globalPosMax, globalPosMin;

  MPI_Reduce(&localPosMax, &globalPosMax, 1, MPI_DOUBLE, MPI_MAX, 0 , comm);
  MPI_Reduce(&localPosMin, &globalPosMin, 1, MPI_DOUBLE, MPI_MIN, 0 , comm);

  ofstream debug;
  if(myProcessRank == 0 ){
    debug.open("postprocessDebug.dat");

    for(CFuint i=0; i< g_cellPos.size();++i){
      debug<<g_cellPos[i]<<' '<<g_Qout[i]<<endl;
    }
    debug.close();

    // calculate the LSF matrices

    gsl_matrix *X, *cov;
    gsl_vector *y, *c;
    CFreal chisq;

    CFuint n = g_cellPos.size();

    X = gsl_matrix_alloc (n, m_nbCoeffs);
    y = gsl_vector_alloc (n);

    c = gsl_vector_alloc (m_nbCoeffs);
    cov = gsl_matrix_alloc (m_nbCoeffs, m_nbCoeffs);

    CFreal xi, yi;
    for (CFuint i = 0; i < n; i++)
    {
      gsl_matrix_set (X, i, 0, 1.0);
      xi = g_cellPos[i];
      for( CFuint j = 1 ; j< m_nbCoeffs ; ++j){
        gsl_matrix_set (X, i, j, xi);
        xi *= xi;
      }
      yi =  g_Qout[i];
      gsl_vector_set (y, i, yi);
    }

    {
      gsl_multifit_linear_workspace * work
          = gsl_multifit_linear_alloc (n, m_nbCoeffs);
      gsl_multifit_linear (X, y, c, cov, &chisq, work);
      gsl_multifit_linear_free (work);
    }

    debug.open("postprocessCoeffs.dat");
    vector<CFreal> coefficients(m_nbCoeffs);
    //output the regression coefficients to file
    CFreal buffer;
    for(CFuint i = 0; i< m_nbCoeffs; ++i){
      buffer=gsl_vector_get(c,i);
      coefficients[i] = buffer;
      debug<<buffer<<' ';
    }
    debug<<endl;
    debug.close();

    //output the plot to file

    //get the plot step
    vector<CFreal> positions(m_nbPointsPlot);
    CFreal step = (globalPosMax-globalPosMin)/m_nbPointsPlot;
    for(CFuint i=0; i<m_nbPointsPlot;++i){
        positions[i] = globalPosMin+step*i;
    }

    debug.open("postprocess.plt");

    CFreal temp,qPlot;
    for(CFuint i=0; i<m_nbPointsPlot;++i){
      temp = positions[i];
      qPlot = coefficients[0];
      for(CFuint j=1; j<m_nbCoeffs;++j){
        qPlot += temp * coefficients[j];
        temp *= temp;
      }
      debug<<positions[i]<<' '<<qPlot<<endl;
    }

    debug.close();

    gsl_matrix_free (X);
    gsl_vector_free (y);
    gsl_vector_free (c);
    gsl_matrix_free (cov);
  }

}

}

}
