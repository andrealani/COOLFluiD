#include "Common/PE.hh"
#include "Common/BadValueException.hh"

#include "MathTools/MathConsts.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PathAppender.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/PhysicalConsts.hh"

#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/GridConvergenceTwoFluid.hh"

#include "TecplotWriter/WriteTecplot.hh"


/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::TecplotWriter;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<GridConvergenceTwoFluid,
    DataProcessingData, FiniteVolumeMultiFluidMHDModule>
GridConvergenceTwoFluidProvider("GridConvergenceTwoFluid");

//////////////////////////////////////////////////////////////////////////////

void GridConvergenceTwoFluid::defineConfigOptions(Config::OptionList& options)
{
    options.addConfigOption< CFuint >("SaveRate","Save Output File every...iterations.");
    options.addConfigOption< bool >("AppendTime","Save each iteration to different file with suffix _time#.");
    options.addConfigOption< bool >("AppendIter","Save each iteration to different file with suffix _iter#.");
    options.addConfigOption< std::string >("OutputFileError","Name of Output File to write the electric and magnetic field divergence");
    options.addConfigOption< bool >("IsTransversal","The wave crosses the domain with an angle");
}

//////////////////////////////////////////////////////////////////////////////

GridConvergenceTwoFluid::GridConvergenceTwoFluid(const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_volumes("volumes"),
  socket_nodes("nodes"),
  socket_isOutward("isOutward"),
  socket_normals("normals"),
  socket_TheorySolution("TheorySolution"),
  m_ErrorL1(),
  m_ErrorL2(),
  m_ErrorTotalL1(),
  m_ErrorTotalL2()
{
  addConfigOptionsTo(this);

  m_saveRate = 1;
  setParameter("SaveRate",&m_saveRate);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);

  m_nameOutputFileError = "Error.plt";
  setParameter("OutputFileError",&m_nameOutputFileError);

  m_isTransversal = false;
  setParameter("IsTransversal",&m_isTransversal);

}

//////////////////////////////////////////////////////////////////////////////

GridConvergenceTwoFluid::~GridConvergenceTwoFluid()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
GridConvergenceTwoFluid::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  
  result.push_back(&socket_TheorySolution);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
GridConvergenceTwoFluid::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nstates);
  result.push_back(&socket_volumes);
  result.push_back(&socket_nodes);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_normals);  
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void GridConvergenceTwoFluid::setup()
{
  CFAUTOTRACE;
  
  // Get number of cells
  Common::SafePtr<Common::ConnectivityTable<CFuint> > cells =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");
  
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();

  const CFuint nbCells = cells->nbRows();

  DataHandle<CFreal> TheorySolution = socket_TheorySolution.getDataHandle();
  TheorySolution.resize(nbCells*nbEqs);
  TheorySolution = 0.0;

  m_ErrorL1.resize(nbEqs, 0.0);
  m_ErrorL2.resize(nbEqs, 0.0);
  m_ErrorTotalL1.resize(nbEqs, 0.0);
  m_ErrorTotalL2.resize(nbEqs, 0.0);


  m_geoBuilder.setup();

// AAL: To add in the future a chemical library to compute this
//  m_library = PhysicalModelStack::getActive()->getImplementor()->
//    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
//  cf_assert(m_library.isNotNull());
  
}

//////////////////////////////////////////////////////////////////////////////

void GridConvergenceTwoFluid::execute()
{
  CFout <<"GridConvergenceTwoFluid::computing error \n";
  CFAUTOTRACE;
  
  const std::string nsp = this->getMethodData().getNamespace();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  DataHandle<Framework::State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq(); 
   
  DataHandle<CFreal> TheorySolution = socket_TheorySolution.getDataHandle();
  
  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  std::vector<SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();

  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells"); 
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  // Computing Analytical Solution
  computeAnalyticalSolution();
  
  for(CFuint i = 0; i < nbEqs; i++){
    m_ErrorL1[i] = 0.;
    m_ErrorL2[i] = 0.;
    m_ErrorTotalL1[i] = 0.;
    m_ErrorTotalL2[i] = 0.; 
  }

  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();
    const CFuint elemID = currCell->getID();

    //Set the states
    State *currState = currCell->getState(0);

    if(currCell->getState(0)->isParUpdatable()) {
      for(CFuint i = 0; i < nbEqs; i++){
        m_ErrorL1[i] += abs((*currState)[i] - TheorySolution[nbEqs*iCell + i]);
        m_ErrorL2[i] += std::pow(abs((*currState)[i] - TheorySolution[nbEqs*iCell + i]),2.);   
      }
    }
    
    m_geoBuilder.releaseGE();
  }
  //Total sum of the processors

  const CFreal totalnbCells = MeshDataStack::getActive()->getTotalStateCount();  

  if (PE::GetPE().GetProcessorCount(nsp) > 1) {

#ifdef CF_HAVE_MPI
    for(CFuint i = 0; i<nbEqs ; i++){
      MPI_Allreduce(&m_ErrorL1[i], &m_ErrorTotalL1[i], 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
      MPI_Allreduce(&m_ErrorL2[i], &m_ErrorTotalL2[i], 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    }
#endif
  }
  else {
    for(CFuint i = 0; i<nbEqs ; i++){ 
      m_ErrorTotalL1[i] = m_ErrorL1[i];
      m_ErrorTotalL2[i] = m_ErrorL2[i];
    }
  }
  for(CFuint i = 0; i<nbEqs ; i++){
    m_ErrorTotalL1[i] /= totalnbCells;
    m_ErrorTotalL2[i] = sqrt(m_ErrorTotalL2[i]/totalnbCells); 
  }


  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  if(!(iter % m_saveRate)) {

    WriteTecplot::getInstance().setNodeExtrapolation(1);
      //we can use default 1
    writeOutputFile();
  }
}

//////////////////////////////////////////////////////////////////////////////

void GridConvergenceTwoFluid::computeAnalyticalSolution()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> TheorySolution = socket_TheorySolution.getDataHandle();
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();

  Common::SafePtr<TopologicalRegionSet> cells =
  MeshDataStack::getActive()->getTrs("InnerCells");
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;

  const CFuint nbCells = cells->getLocalNbGeoEnts();
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
  CFreal t = SubSystemStatusStack::getActive()->getCurrentTimeDim();; 

  if(!m_isTransversal) {
    for (CFuint iCell = 0; iCell < nbCells; ++iCell){
      //Set the cartesian and radial coordinates
      geoData.idx = iCell;
      GeometricEntity* currCell = m_geoBuilder.buildGE();
      const CFuint elemID = currCell->getID();
      Node& coordinate = currCell->getState(0)->getCoordinates();
      const CFreal x = coordinate[XX];
      const CFreal y = coordinate[YY];
    
      TheorySolution[iCell*nbEqs + 0] = 6e-09; //Bx
      TheorySolution[iCell*nbEqs + 1] = 6e-10*cos(3.10528461412e-05*x - 3.121781382662538*t); //By
      TheorySolution[iCell*nbEqs + 2] = 6e-10*sin(3.10528461412e-05*x - 3.121781382662538*t); //Bz
      TheorySolution[iCell*nbEqs + 3] = 0.0; //Ex
      TheorySolution[iCell*nbEqs + 4] = -6.03187489185e-05*sin(3.10528461412e-05*x - 3.121781382662538*t); //Ey
      TheorySolution[iCell*nbEqs + 5] = 6.03187489185e-05*cos(3.10528461412e-05*x - 3.121781382662538*t); //Ez
      TheorySolution[iCell*nbEqs + 6] = 0.0; //Psi
      TheorySolution[iCell*nbEqs + 7] = 0.0; //Phi
      TheorySolution[iCell*nbEqs + 8] = 6.5336788193e-23; //rhoe
      TheorySolution[iCell*nbEqs + 9] = 1.67262177774e-20; //rhoi
      TheorySolution[iCell*nbEqs + 10] = 0.0; //Ue
      TheorySolution[iCell*nbEqs + 11] = 10271.0531054*cos(3.10528461412e-05*x - 3.121781382662538*t); //Ve
      TheorySolution[iCell*nbEqs + 12] = 10271.0531054*sin(3.10528461412e-05*x - 3.121781382662538*t); //We
      TheorySolution[iCell*nbEqs + 13] = 0.0; //Ui
      TheorySolution[iCell*nbEqs + 14] = 1563.05016333*cos(3.10528461412e-05*x - 3.121781382662538*t); //Vi
      TheorySolution[iCell*nbEqs + 15] = 1563.05016333*sin(3.10528461412e-05*x - 3.121781382662538*t); //Wi
      TheorySolution[iCell*nbEqs + 16] = 5187.39732495; //Te
      TheorySolution[iCell*nbEqs + 17] = 5187.39732495; //Ti
      m_geoBuilder.releaseGE();
    }
  }
  else {
    for (CFuint iCell = 0; iCell < nbCells; ++iCell){
      geoData.idx = iCell;
      GeometricEntity* currCell = m_geoBuilder.buildGE();
      const CFuint elemID = currCell->getID();
      Node& coordinate = currCell->getState(0)->getCoordinates();
      const CFreal x = coordinate[XX];
      const CFreal y = coordinate[YY];
      const CFreal sinalpha = sqrt(4./5.);
      const CFreal cosalpha = sqrt(1./5.);
      const CFreal xStar = (x*cosalpha+y*sinalpha);
      const CFreal Bpar  = 6e-09;
      const CFreal Bperp = 6e-10*cos(3.10528461412e-05*xStar - 3.121781382662538*t);
      const CFreal Epar  = 0.;
      const CFreal Eperp = -6.03187489185e-05*sin(3.10528461412e-05*xStar - 3.121781382662538*t);
      const CFreal Uepar = 0.;
      const CFreal Ueperp= 10271.0531054*cos(3.10528461412e-05*xStar - 3.121781382662538*t);
      const CFreal Uipar = 0.;
      const CFreal Uiperp= 1563.05016333*cos(3.10528461412e-05*xStar - 3.121781382662538*t);
    
      TheorySolution[iCell*nbEqs + 0] = Bpar*cosalpha - Bperp*sinalpha; //Bx
      TheorySolution[iCell*nbEqs + 1] = Bpar*sinalpha + Bperp*cosalpha; //By
      TheorySolution[iCell*nbEqs + 2] = 6e-10*sin(3.10528461412e-05*xStar - 3.121781382662538*t); //Bz
      TheorySolution[iCell*nbEqs + 3] = Epar*cosalpha - Eperp*sinalpha; //Ex
      TheorySolution[iCell*nbEqs + 4] = Epar*sinalpha + Eperp*cosalpha; //Ey
      TheorySolution[iCell*nbEqs + 5] = 6.03187489185e-05*cos(3.10528461412e-05*xStar - 3.121781382662538*t); //Ez
      TheorySolution[iCell*nbEqs + 6] = 0.0; //Psi
      TheorySolution[iCell*nbEqs + 7] = 0.0; //Phi
      TheorySolution[iCell*nbEqs + 8] = 6.5336788193e-23; //rhoe
      TheorySolution[iCell*nbEqs + 9] = 1.67262177774e-20; //rhoi
      TheorySolution[iCell*nbEqs + 10] = Uepar*cosalpha - Ueperp*sinalpha; //Ue
      TheorySolution[iCell*nbEqs + 11] = Uepar*sinalpha + Ueperp*cosalpha; //Ve
      TheorySolution[iCell*nbEqs + 12] = 10271.0531054*sin(3.10528461412e-05*xStar - 3.121781382662538*t); //We
      TheorySolution[iCell*nbEqs + 13] = Uipar*cosalpha - Uiperp*sinalpha; //Ui
      TheorySolution[iCell*nbEqs + 14] = Uipar*sinalpha + Uiperp*cosalpha; //Vi
      TheorySolution[iCell*nbEqs + 15] = 1563.05016333*sin(3.10528461412e-05*xStar - 3.121781382662538*t); //Wi
      TheorySolution[iCell*nbEqs + 16] = 5187.39732495; //Te
      TheorySolution[iCell*nbEqs + 17] = 5187.39732495; //Ti
      m_geoBuilder.releaseGE();
    }   
  }
}

/////////////////////////////////////////////////////////////////////////////

boost::filesystem::path GridConvergenceTwoFluid::constructFilename()
{
  const bool isParallel = PE::GetPE().IsParallel ();
  const std::string nsp = this->getMethodData().getNamespace();

  if (isParallel) {
    std::ostringstream fname;
#ifdef CF_HAVE_BOOST_1_85
    fname << boost::filesystem::path(m_nameOutputFileError).stem().string()
          << "-" << PE::GetPE().GetRank("Default")
          << boost::filesystem::path(m_nameOutputFileError).extension().string();
#else
    fname << boost::filesystem::basename(boost::filesystem::path(m_nameOutputFileError))
          << "-" << PE::GetPE().GetRank("Default")
          << boost::filesystem::extension(boost::filesystem::path(m_nameOutputFileError));
#endif
  }

  CFout << "Writing Error file to: " << m_nameOutputFileError << "\n";
  return boost::filesystem::path(m_nameOutputFileError);
}

//////////////////////////////////////////////////////////////////////////////
void GridConvergenceTwoFluid::prepareOutputFile(std::ofstream& outputFile)
{
  outputFile << "TITLE  =  Error Monitoring\n";
  outputFile << "VARIABLES = Iter Bx-L1 By-L1 Bz-L1 Ex-L1 Ey-L1 Ez-L1 Psi-L1 Phi-L1 Rhoe-L1 Rhoi-L1 Uex-L1 Uey-L1 Uez-L1 Uix-L1 Uiy-L1 Uiz-L1 Te-L1 Ti-L1 Bx-L2 By-L2 Bz-L2 Ex-L2 Ey-L2 Ez-L2 Psi-L2 Phi-L2 Rhoe-L2 Rhoi-L2 Uex-L2 Uey-L2 Uez-L2 Uix-L2 Uiy-L2 Uiz-L2 Te-L2 Ti-L2 Time nbCells\n";
}

//////////////////////////////////////////////////////////////////////////////

void GridConvergenceTwoFluid::writeOutputFile()
{
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;

  const std::string nsp = this->getMethodData().getNamespace();
  MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
  int nbProcessors, myRank;
  MPI_Comm_size(comm, &nbProcessors);
  MPI_Comm_rank(comm, &myRank);

  if (myRank==0) {
//  const CFuint nbCells = cells->getLocalNbGeoEnts();
  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  const CFreal CurrentTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
  const CFreal totalnbCells = MeshDataStack::getActive()->getTotalStateCount();

  //std::cout << "IN writeOutputFile() \n";

  SelfRegistPtr<Environment::FileHandlerOutput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();
  
  const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();

    if (iter == 1) {
      ofstream& outputFile = (*fhandle)->open(constructFilename());
      prepareOutputFile(outputFile);
      outputFile << iter <<" ";
      for(CFuint i = 0; i < nbEqs; i++){
        outputFile << m_ErrorTotalL1[i] <<" ";
      }
      for(CFuint i = 0; i < nbEqs; i++){
        outputFile << m_ErrorTotalL2[i] <<" ";
      }
      outputFile << CurrentTime << " " << totalnbCells << "\n";
      outputFile.close();
    }
    else {
      ofstream& outputFile = (*fhandle)->open(constructFilename(), ios::app);
      outputFile << iter <<" ";
      for(CFuint i = 0; i < nbEqs; i++){
        outputFile << m_ErrorTotalL1[i] <<" ";
      }
      for(CFuint i = 0; i < nbEqs; i++){
        outputFile << m_ErrorTotalL2[i] <<" ";
      }
      outputFile << CurrentTime << " " << totalnbCells << "\n";
      outputFile.close();
    }
    delete fhandle;
  }
  CFLog(VERBOSE, "Writing Error file finished." << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void GridConvergenceTwoFluid::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

