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
#include "FiniteVolumeMultiFluidMHD/GridConvergence.hh"

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

MethodCommandProvider<GridConvergence,
    DataProcessingData, FiniteVolumeMultiFluidMHDModule>
GridConvergenceProvider("GridConvergence");

//////////////////////////////////////////////////////////////////////////////

void GridConvergence::defineConfigOptions(Config::OptionList& options)
{
    options.addConfigOption< CFuint >("SaveRate","Save Output File every...iterations.");
    options.addConfigOption< bool >("AppendTime","Save each iteration to different file with suffix _time#.");
    options.addConfigOption< bool >("AppendIter","Save each iteration to different file with suffix _iter#.");
    options.addConfigOption< std::string >("OutputFileError","Name of Output File to write the electric and magnetic field divergence");
}

//////////////////////////////////////////////////////////////////////////////

GridConvergence::GridConvergence(const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_volumes("volumes"),
  socket_nodes("nodes"),
  socket_isOutward("isOutward"),
  socket_normals("normals"),
  socket_BxTheory("BxTheory"),
  socket_ByTheory("ByTheory"),
  socket_EzTheory("EzTheory"),
  socket_RhoTheory("RhoTheory"),
  socket_UxTheory("UxTheory"),
  socket_UyTheory("UyTheory"),
  socket_TTheory("TTheory")
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

}

//////////////////////////////////////////////////////////////////////////////

GridConvergence::~GridConvergence()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
GridConvergence::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  
  result.push_back(&socket_BxTheory);
  result.push_back(&socket_ByTheory);
  result.push_back(&socket_EzTheory);
  result.push_back(&socket_RhoTheory);
  result.push_back(&socket_UxTheory);
  result.push_back(&socket_UyTheory);
  result.push_back(&socket_TTheory);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
GridConvergence::needsSockets()
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

void GridConvergence::setup()
{
  CFAUTOTRACE;

  DataProcessingCom::setup();
  
  // Get number of cells
  Common::SafePtr<Common::ConnectivityTable<CFuint> > cells =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  const CFuint nbCells = cells->nbRows();

  DataHandle<CFreal> BxTheory = socket_BxTheory.getDataHandle();
  BxTheory.resize(nbCells);
  BxTheory = 0.0;

  DataHandle<CFreal> ByTheory = socket_ByTheory.getDataHandle();
  ByTheory.resize(nbCells);
  ByTheory = 0.0;

  DataHandle<CFreal> EzTheory = socket_EzTheory.getDataHandle();
  EzTheory.resize(nbCells);
  EzTheory = 0.0;

  DataHandle<CFreal> RhoTheory = socket_RhoTheory.getDataHandle();
  RhoTheory.resize(nbCells);
  RhoTheory = 0.0;

  DataHandle<CFreal> UxTheory = socket_UxTheory.getDataHandle();
  UxTheory.resize(nbCells);
  UxTheory = 0.0;
  
  DataHandle<CFreal> UyTheory = socket_UyTheory.getDataHandle();
  UyTheory.resize(nbCells);
  UyTheory = 0.0;

  DataHandle<CFreal> TTheory = socket_TTheory.getDataHandle();
  TTheory.resize(nbCells);
  TTheory = 0.0;

  m_geoBuilder.setup();
  
  // AAL: To add in the future a chemical library to compute this
  //  m_library = PhysicalModelStack::getActive()->getImplementor()->
  //    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  //  cf_assert(m_library.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void GridConvergence::execute()
{
  CFout <<"GridConvergence::computing properties \n";
  CFAUTOTRACE;
  
  const std::string nsp = this->getMethodData().getNamespace();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  DataHandle<Framework::State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  DataHandle<CFreal> BxTheory = socket_BxTheory.getDataHandle();
  DataHandle<CFreal> ByTheory = socket_ByTheory.getDataHandle();
  DataHandle<CFreal> EzTheory = socket_EzTheory.getDataHandle();
  DataHandle<CFreal> RhoTheory = socket_RhoTheory.getDataHandle();
  DataHandle<CFreal> UxTheory = socket_UxTheory.getDataHandle();
  DataHandle<CFreal> UyTheory = socket_UyTheory.getDataHandle();
  DataHandle<CFreal> TTheory = socket_TTheory.getDataHandle();

  
  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  std::vector<SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();

  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells"); 
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  // Computing Analytical Solution
  computeAnalyticalSolution();
  m_ErrorBxL1 = 0.;
  m_ErrorByL1 = 0.;
  m_ErrorEzL1 = 0.;
  m_ErrorRhoL1 = 0.;
  m_ErrorUxL1 = 0.;
  m_ErrorUyL1 = 0.;
  m_ErrorTL1 = 0.;
  m_ErrorTotalL1 = 0.;

  m_ErrorBxL2 = 0.;
  m_ErrorByL2 = 0.;
  m_ErrorEzL2 = 0.;
  m_ErrorRhoL2 = 0.;
  m_ErrorUxL2 = 0.;
  m_ErrorUyL2 = 0.;
  m_ErrorTL2 = 0.;
  m_ErrorTotalL2 =0.;


  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();
    const CFuint elemID = currCell->getID();

    CFreal Bx = 0.;
    CFreal By = 0.;
    CFreal Ez = 0.;
    CFreal Rho = 0.;
    CFreal Ux = 0.;
    CFreal Uy = 0.;
    CFreal T = 0.;

    //Set the states
    State *currState = currCell->getState(0);
    Bx = (*currState)[0];
    By = (*currState)[1];
    Ez = (*currState)[5];
    Rho = (*currState)[8];
    Ux = (*currState)[9];
    Uy = (*currState)[10];
    T = (*currState)[11];


    if(currCell->getState(0)->isParUpdatable()) {
      m_ErrorBxL1 += abs(Bx - BxTheory[iCell]);
      m_ErrorByL1 += abs(By - ByTheory[iCell]);
      m_ErrorEzL1 += abs(Ez - EzTheory[iCell]);
      m_ErrorRhoL1 += abs(Rho - RhoTheory[iCell]);
      m_ErrorUxL1 += abs(Ux - UxTheory[iCell]);
      m_ErrorUyL1 += abs(Uy - UyTheory[iCell]);
      m_ErrorTL1 += abs(T - TTheory[iCell]);

      m_ErrorBxL2 += abs(Bx*Bx - BxTheory[iCell]*BxTheory[iCell]);
      m_ErrorByL2 += abs(By*By - ByTheory[iCell]*ByTheory[iCell]);
      m_ErrorEzL2 += abs(Ez*Ez - EzTheory[iCell]*EzTheory[iCell]);
      m_ErrorRhoL2 += abs(Rho*Rho - RhoTheory[iCell]*RhoTheory[iCell]);
      m_ErrorUxL2 += abs(Ux*Ux - UxTheory[iCell]*UxTheory[iCell]);
      m_ErrorUyL2 += abs(Uy*Uy - UyTheory[iCell]*UyTheory[iCell]);
      m_ErrorTL2 += abs(T*T - TTheory[iCell]*TTheory[iCell]);
    }


    //m_ErrorBxL1 *= volumes[elemID];
    //m_ErrorByL1 *= volumes[elemID];
    //m_ErrorUxL1 *= volumes[elemID];
    //m_ErrorUyL1 *= volumes[elemID];

    //m_ErrorBxL2 *= volumes[elemID];
    //m_ErrorByL2 *= volumes[elemID];
    //m_ErrorUxL2 *= volumes[elemID];
    //m_ErrorUyL2 *= volumes[elemID];
    
    m_geoBuilder.releaseGE();
  }
  //Total sum of the processors
  CFdouble totalErrorBxL1 = 0.;
  CFdouble totalErrorByL1 = 0.;
  CFdouble totalErrorEzL1 = 0.;
  CFdouble totalErrorRhoL1 = 0.;
  CFdouble totalErrorUxL1 = 0.;
  CFdouble totalErrorUyL1 = 0.;
  CFdouble totalErrorTL1 = 0.;

  CFdouble totalErrorBxL2 = 0.;
  CFdouble totalErrorByL2 = 0.;
  CFdouble totalErrorEzL2 = 0.;
  CFdouble totalErrorRhoL2 = 0.;
  CFdouble totalErrorUxL2 = 0.;
  CFdouble totalErrorUyL2 = 0.;
  CFdouble totalErrorTL2 = 0.;

  if (PE::GetPE().GetProcessorCount(nsp) > 1) {

#ifdef CF_HAVE_MPI
    MPI_Allreduce(&m_ErrorBxL1, &totalErrorBxL1, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorByL1, &totalErrorByL1, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorEzL1, &totalErrorEzL1, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorRhoL1, &totalErrorRhoL1, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorUxL1, &totalErrorUxL1, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorUyL1, &totalErrorUyL1, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorTL1, &totalErrorTL1, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorBxL2, &totalErrorBxL2, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorByL2, &totalErrorByL2, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorEzL2, &totalErrorEzL2, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorRhoL2, &totalErrorRhoL2, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorUxL2, &totalErrorUxL2, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorUyL2, &totalErrorUyL2, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&m_ErrorTL2, &totalErrorTL2, 1, MPI_DOUBLE, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));

#endif
  }
  else {
    totalErrorBxL1 = m_ErrorBxL1;
    totalErrorByL1 = m_ErrorByL1;
    totalErrorEzL1 = m_ErrorEzL1;
    totalErrorRhoL1 = m_ErrorRhoL1;
    totalErrorUxL1 = m_ErrorUxL1;
    totalErrorUyL1 = m_ErrorUyL1;
    totalErrorTL1 = m_ErrorTL1;
    totalErrorBxL2 = m_ErrorBxL2;
    totalErrorByL2 = m_ErrorByL2;
    totalErrorEzL2 = m_ErrorEzL2;
    totalErrorRhoL2 = m_ErrorRhoL2;
    totalErrorUxL2 = m_ErrorUxL2;
    totalErrorUyL2 = m_ErrorUyL2;
    totalErrorTL2 = m_ErrorTL2;
  }

  const CFreal totalnbCells = MeshDataStack::getActive()->getTotalStateCount();
  const CFreal B_ref = 0.112099824326;
  const CFreal E_ref = 11.2099824326;
  const CFreal Rho_ref = 1.;
  const CFreal U_ref = 100.;
  const CFreal T_ref = 1211.47544325;

  totalErrorBxL1 /= totalnbCells;
  totalErrorByL1 /= totalnbCells;
  totalErrorEzL1 /= totalnbCells;
  totalErrorRhoL1 /= totalnbCells;
  totalErrorUxL1 /= totalnbCells;
  totalErrorUyL1 /= totalnbCells;
  totalErrorTL1 /= totalnbCells;
  //Total error needs to be everything in non-dimensional form
  m_ErrorTotalL1 = (totalErrorBxL1 + totalErrorByL1)/B_ref
          + totalErrorEzL1/E_ref + totalErrorRhoL1/Rho_ref
          + (totalErrorUxL1 + totalErrorUyL1)/U_ref + totalErrorTL1/T_ref;

  //Second norm is done different
  m_ErrorTotalL2 = (totalErrorBxL2 + totalErrorByL2)/(B_ref*B_ref)
          + totalErrorEzL2/(E_ref*E_ref) + totalErrorRhoL2/(Rho_ref*Rho_ref)
          + (totalErrorUxL2 + totalErrorUyL2)/(U_ref*U_ref) + totalErrorTL2/(T_ref*T_ref);
  m_ErrorTotalL2 = sqrt(m_ErrorTotalL2)/totalnbCells;

  totalErrorBxL2 = sqrt(totalErrorBxL2)/totalnbCells;
  totalErrorByL2 = sqrt(totalErrorByL2)/totalnbCells;
  totalErrorEzL2 = sqrt(totalErrorEzL2)/totalnbCells;
  totalErrorRhoL2 = sqrt(totalErrorRhoL2)/totalnbCells;
  totalErrorUxL2 = sqrt(totalErrorUxL2)/totalnbCells;
  totalErrorUyL2 = sqrt(totalErrorUyL2)/totalnbCells;
  totalErrorTL2 = sqrt(totalErrorTL2)/totalnbCells;

  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  if(!(iter % m_saveRate)) {

    m_ErrorBxL1 = totalErrorBxL1;
    m_ErrorByL1 = totalErrorByL1;
    m_ErrorEzL1 = totalErrorEzL1;
    m_ErrorRhoL1 = totalErrorRhoL1;
    m_ErrorUxL1 = totalErrorUxL1;
    m_ErrorUyL1 = totalErrorUyL1;
    m_ErrorTL1 = totalErrorTL1;

    m_ErrorBxL2 = totalErrorBxL2;
    m_ErrorByL2 = totalErrorByL2;
    m_ErrorEzL2 = totalErrorEzL2;
    m_ErrorRhoL2 = totalErrorRhoL2;
    m_ErrorUxL2 = totalErrorUxL2;
    m_ErrorUyL2 = totalErrorUyL2;
    m_ErrorTL2 = totalErrorTL2;

    WriteTecplot::getInstance().setNodeExtrapolation(1);
      //we can use default 1
    writeOutputFile();
  }
}

//////////////////////////////////////////////////////////////////////////////

void GridConvergence::computeAnalyticalSolution()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> BxTheory = socket_BxTheory.getDataHandle();
  DataHandle<CFreal> ByTheory = socket_ByTheory.getDataHandle();
  DataHandle<CFreal> EzTheory = socket_EzTheory.getDataHandle();
  DataHandle<CFreal> RhoTheory = socket_RhoTheory.getDataHandle();
  DataHandle<CFreal> UxTheory = socket_UxTheory.getDataHandle();
  DataHandle<CFreal> UyTheory = socket_UyTheory.getDataHandle();
  DataHandle<CFreal> TTheory = socket_TTheory.getDataHandle();
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();

  Common::SafePtr<TopologicalRegionSet> cells =
  MeshDataStack::getActive()->getTrs("InnerCells");
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;

  const CFuint nbCells = cells->getLocalNbGeoEnts();

  for (CFuint iCell = 0; iCell < nbCells; ++iCell){
    //Set the cartesian and radial coordinates
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();
    const CFuint elemID = currCell->getID();
    Node& coordinate = currCell->getState(0)->getCoordinates();
    const CFreal x = coordinate[XX];
    const CFreal y = coordinate[YY];
    const CFreal r2 = (x - 500.)*(x - 500.) + (y - 500.)*(y - 500.);
    const CFreal r = std::sqrt(r2)/100.;

    BxTheory[iCell] = -0.000178412411613*(y-500.)*std::exp(0.5*(1-r*r));
    ByTheory[iCell] = 0.000178412411613*(x-500.)*std::exp(0.5*(1-r*r));
    EzTheory[iCell] = (100.+0.159154943092*(x-500.)*exp(0.5*(1-r*r)))*(-0.000178412411613*(y-500.)*exp(0.5*(1-r*r)))
            -(100.-0.159154943092*(y-500.)*exp(0.5*(1-r*r)))*(0.000178412411613*(x-500.)*exp(0.5*(1-r*r)));
    RhoTheory[iCell] = 1.0;
    UxTheory[iCell] = 100.-0.159154943092*(y-500.)*std::exp(0.5*(1-r*r));
    UyTheory[iCell] = 100.+0.159154943092*(x-500.)*std::exp(0.5*(1-r*r));
    TTheory[iCell] = 1211.47544325+15.343515733*(1-r*r)*exp((1-r*r))-15.343515733*exp((1-r*r));
    m_geoBuilder.releaseGE();
  }
}

/////////////////////////////////////////////////////////////////////////////

boost::filesystem::path GridConvergence::constructFilename()
{
  const bool isParallel = PE::GetPE().IsParallel ();
  const std::string nsp = this->getMethodData().getNamespace();

  if (isParallel) {
    std::ostringstream fname;
    fname << boost::filesystem::basename(boost::filesystem::path(m_nameOutputFileError))
          << "-" << PE::GetPE().GetRank("Default")
          << boost::filesystem::extension(boost::filesystem::path(m_nameOutputFileError));
  }

  CFout << "Writing Error file to: " << m_nameOutputFileError << "\n";
  return boost::filesystem::path(m_nameOutputFileError);
}

//////////////////////////////////////////////////////////////////////////////
void GridConvergence::prepareOutputFile(std::ofstream& outputFile)
{
  outputFile << "TITLE  =  Error Monitoring\n";
  outputFile << "VARIABLES = Iter Bx-L1 By-L1 Ez-L1 Rho-L1 Ux-L1 Uy-L1 T-L1 Total-L1 Bx-L2 By-L2 Ez-L2 Rho-L2 Ux-L2 Uy-L2 T-L2 Total-L2 Time nbCells\n";
}

//////////////////////////////////////////////////////////////////////////////

void GridConvergence::writeOutputFile()
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

    if (iter == 1) {
      ofstream& outputFile = (*fhandle)->open(constructFilename());
      prepareOutputFile(outputFile);
      outputFile << iter
               << " "
               << m_ErrorBxL1
               << " "
               << m_ErrorByL1
               << " "
               << m_ErrorEzL1
               << " "
               << m_ErrorRhoL1
               << " "
               << m_ErrorUxL1
               << " "
               << m_ErrorUyL1
               << " "
               << m_ErrorTL1
               << " "
               << m_ErrorTotalL1
               << " "
               << m_ErrorBxL2
               << " "
               << m_ErrorByL2
               << " "
               << m_ErrorEzL2
               << " "
               << m_ErrorRhoL2
               << " "
               << m_ErrorUxL2
               << " "
               << m_ErrorUyL2
               << " "
               << m_ErrorTL2
               << " "
               << m_ErrorTotalL2
               << " "
               << CurrentTime
               << " "
               << totalnbCells
               << "\n";
      outputFile.close();
    }
    else {
      ofstream& outputFile = (*fhandle)->open(constructFilename(), ios::app);
      outputFile << iter
                 << " "
                 << m_ErrorBxL1
                 << " "
                 << m_ErrorByL1
                 << " "
                 << m_ErrorEzL1
                 << " "
                 << m_ErrorRhoL1
                 << " "
                 << m_ErrorUxL1
                 << " "
                 << m_ErrorUyL1
                 << " "
                 << m_ErrorTL1
                 << " "
                 << m_ErrorTotalL1
                 << " "
                 << m_ErrorBxL2
                 << " "
                 << m_ErrorByL2
                 << " "
                 << m_ErrorEzL2
                 << " "
                 << m_ErrorRhoL2
                 << " "
                 << m_ErrorUxL2
                 << " "
                 << m_ErrorUyL2
                 << " "
                 << m_ErrorTL2
                 << " "
                 << m_ErrorTotalL2
                 << " "
                 << CurrentTime
                 << " "
                 << totalnbCells
                 << "\n";
      outputFile.close();
    }
    delete fhandle;
  }
  CFout << "Writing of Electric and Magnetic Field divergence finished." << "\n";
//  outputFile.close();
}

//////////////////////////////////////////////////////////////////////////////

void GridConvergence::unsetup()
{
  CFAUTOTRACE;

  DataProcessingCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

