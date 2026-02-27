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
#include "FiniteVolumeMultiFluidMHD/EnergyIntegration.hh"

#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"

#include "TecplotWriter/WriteTecplot.hh"


/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::TecplotWriter;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace COOLFluiD::Physics::Maxwell;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<EnergyIntegration,
    DataProcessingData, FiniteVolumeMultiFluidMHDModule>
EnergyIntegrationProvider("EnergyIntegration");

//////////////////////////////////////////////////////////////////////////////

void EnergyIntegration::defineConfigOptions(Config::OptionList& options)
{
    options.addConfigOption< CFuint >("SaveRate","Save Output File every...iterations.");
    options.addConfigOption< bool >("AppendTime","Save each iteration to different file with suffix _time#.");
    options.addConfigOption< bool >("AppendIter","Save each iteration to different file with suffix _iter#.");
    options.addConfigOption< std::string >("OutputFileError","Name of Output File to write the electric and magnetic field divergence");
}

//////////////////////////////////////////////////////////////////////////////

EnergyIntegration::EnergyIntegration(const std::string& name) :
  //m_updateVarSet(CFNULL),
  DataProcessingCom(name),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_volumes("volumes"),
  socket_nodes("nodes"),
  socket_isOutward("isOutward"),
  socket_normals("normals")
{
  addConfigOptionsTo(this);

  m_saveRate = 1;
  setParameter("SaveRate",&m_saveRate);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);

  m_nameOutputFileError = "EnergyIntegration.plt";
  setParameter("OutputFileError",&m_nameOutputFileError);

}

//////////////////////////////////////////////////////////////////////////////

EnergyIntegration::~EnergyIntegration()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
EnergyIntegration::needsSockets()
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

void EnergyIntegration::setup()
{
  CFAUTOTRACE;
  
  //m_updateVarSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();


  m_geoBuilder.setup();

// AAL: To add in the future a chemical library to compute this
//  m_library = PhysicalModelStack::getActive()->getImplementor()->
//    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
//  cf_assert(m_library.isNotNull());
  
}

//////////////////////////////////////////////////////////////////////////////

void EnergyIntegration::execute()
{
  CFout <<"EnergyIntegration::computing properties \n";
  CFAUTOTRACE;
  
  const std::string nsp = this->getMethodData().getNamespace();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  DataHandle<Framework::State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  std::vector<SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();

  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells"); 
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  CFreal MagneticEner = 0.;
  CFreal ElectricEner = 0.;
  CFreal PsiEner = 0.;
  CFreal PhiEner = 0.;
  CFreal KinElecEner = 0.;
  CFreal KinIonsEner = 0.;
  CFreal IntElecEner = 0.;
  CFreal IntIonsEner = 0.;
  CFreal totalEMEner = 0.;
  CFreal totalElecEner = 0.;
  CFreal totalIonsEner = 0.;
  CFreal TotalEner = 0.;
  CFuint totalNbCells = 0.;

  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();
    const CFuint elemID = currCell->getID();
    
    //Set the states
    State *currState = currCell->getState(0);
    CFreal Bx   = (*currState)[0];
    CFreal By   = (*currState)[1];
    CFreal Bz   = (*currState)[2];
    CFreal Ex   = (*currState)[3];
    CFreal Ey   = (*currState)[4];
    CFreal Ez   = (*currState)[5];
    CFreal Psi  = (*currState)[6];
    CFreal Phi  = (*currState)[7];
    CFreal Rho0 = (*currState)[8];
    CFreal Rho1 = (*currState)[9];
    CFreal Ux0  = (*currState)[10];
    CFreal Uy0  = (*currState)[11];
    CFreal Uz0  = (*currState)[12];
    CFreal Ux1  = (*currState)[13];
    CFreal Uy1  = (*currState)[14];
    CFreal Uz1  = (*currState)[15];
    CFreal T0   = (*currState)[16];
    CFreal T1   = (*currState)[17];
    //mu_0 = m_updateVarSet->getModel()->getPermeability();
    const CFreal mu_0      = Framework::PhysicalConsts::VacuumPermeability();
    const CFreal c2        = 413854.274193*413854.274193;
    const CFreal epsilon_0 = 1./(mu_0*c2);
    const CFreal chi       = 1.;
    const CFreal gamma     = 1.; 

    if(currCell->getState(0)->isParUpdatable()) {
      MagneticEner += 1./(2.*mu_0)*(Bx*Bx + By*By + Bz*Bz);
      ElectricEner += 1./2.*epsilon_0*(Ex*Ex + Ey*Ey + Ez*Ez);
      PsiEner += 1./2.*gamma*epsilon_0*Psi*Psi;
      PhiEner += chi/(2.*mu_0)*Phi*Phi;
      KinElecEner += 1./2.*Rho0*(Ux0*Ux0 + Uy0*Uy0 + Uz0*Uz0);
      KinIonsEner += 1./2.*Rho1*(Ux1*Ux1 + Uy1*Uy1 + Uz1*Uz1);
      IntElecEner += Rho0*Framework::PhysicalConsts::Boltzmann()/6.69048711096e-29*T0/(5./3. - 1.);
      IntIonsEner += Rho1*Framework::PhysicalConsts::Boltzmann()/1.67262177774e-27*T1/(5./3. - 1.);
      totalNbCells += 1;

    }
    totalEMEner   = MagneticEner + ElectricEner + PsiEner + PhiEner;
    totalElecEner = KinElecEner + IntElecEner;
    totalIonsEner = KinIonsEner + IntIonsEner;
    TotalEner = MagneticEner + ElectricEner + PsiEner + PhiEner + KinElecEner + KinIonsEner + IntElecEner + IntIonsEner;

    m_geoBuilder.releaseGE();
  }
  //Total sum of the processors
  _totalMagneticEner = 0.;
  _totalElectricEner = 0.;
  _totalPsiEner = 0.;
  _totalPhiEner = 0.;
  _totalKinElecEner = 0.;
  _totalKinIonsEner = 0.;
  _totalIntElecEner = 0.;
  _totalIntIonsEner = 0.;
  _totalTotalEMEner = 0.;
  _totalTotalElecEner = 0.;
  _totalTotalIonsEner = 0.;
  _totalTotalEner = 0.;
  _totalTotalNbCells = 0.;

  if (PE::GetPE().GetProcessorCount(nsp) > 1) {

#ifdef CF_HAVE_MPI
    MPI_Barrier(PE::GetPE().GetCommunicator(nsp));

    MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&MagneticEner);
    MPI_Datatype MPI_CFUINT = Common::MPIStructDef::getMPIType(&totalNbCells);

    MPI_Allreduce(&MagneticEner, &_totalMagneticEner, 1, MPI_CFREAL, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&ElectricEner, &_totalElectricEner, 1, MPI_CFREAL, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&PsiEner, &_totalPsiEner, 1, MPI_CFREAL, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&PhiEner, &_totalPhiEner, 1, MPI_CFREAL, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&KinElecEner, &_totalKinElecEner, 1, MPI_CFREAL, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&KinIonsEner, &_totalKinIonsEner, 1, MPI_CFREAL, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&IntElecEner, &_totalIntElecEner, 1, MPI_CFREAL, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&IntIonsEner, &_totalIntIonsEner, 1, MPI_CFREAL, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&totalEMEner, &_totalTotalEMEner, 1, MPI_CFREAL, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&totalElecEner, &_totalTotalElecEner, 1, MPI_CFREAL, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&totalIonsEner, &_totalTotalIonsEner, 1, MPI_CFREAL, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&TotalEner, &_totalTotalEner, 1, MPI_CFREAL, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));
    MPI_Allreduce(&totalNbCells, &_totalTotalNbCells, 1, MPI_CFUINT, MPI_SUM,
              PE::GetPE().GetCommunicator(nsp));

    MPI_Barrier(PE::GetPE().GetCommunicator(nsp));
#endif
  }
  else {
    _totalMagneticEner = MagneticEner;
    _totalElectricEner = ElectricEner;
    _totalPsiEner = PsiEner;
    _totalPhiEner = PhiEner;
    _totalKinElecEner = KinElecEner;
    _totalKinIonsEner = KinIonsEner;
    _totalIntElecEner = IntElecEner;
    _totalIntIonsEner = IntIonsEner;
    _totalTotalEMEner = totalEMEner;
    _totalTotalElecEner = totalElecEner;
    _totalTotalIonsEner = totalIonsEner;
    _totalTotalEner = TotalEner;
    _totalTotalNbCells = totalNbCells;
  }

  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  if(!(iter % m_saveRate)) {
    
    WriteTecplot::getInstance().setNodeExtrapolation(1);
    //we can use default 1
    writeOutputFile();
  }
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path EnergyIntegration::constructFilename()
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

  CFout << "Writing Energy monitoring file to: " << m_nameOutputFileError << "\n";
  return boost::filesystem::path(m_nameOutputFileError);
}

//////////////////////////////////////////////////////////////////////////////
void EnergyIntegration::prepareOutputFile(std::ofstream& outputFile)
{
  outputFile << "TITLE  =  Energy Monitoring\n";
  outputFile << "VARIABLES = Iter MagneticEnergy ElectricEnergy PsiEnergy PhiEnergy ElecKinEnergy IonsKinEnergy ElecIntEnergy IonIntEnergy TotalEMEnergy TotalElecEnergy TotalIonsEnergy TotalEnergy Time nbCells\n";
}

//////////////////////////////////////////////////////////////////////////////

void EnergyIntegration::writeOutputFile()
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
    CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
    const CFreal CurrentTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();

    SelfRegistPtr<Environment::FileHandlerOutput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();

    if (iter == 1) {
      ofstream& outputFile = (*fhandle)->open(constructFilename());
      prepareOutputFile(outputFile);
      outputFile << iter
               << " "
               << _totalMagneticEner
               << " "
               << _totalElectricEner
               << " "
               << _totalPsiEner
               << " "
               << _totalPhiEner
               << " "
               << _totalKinElecEner
               << " "
               << _totalKinIonsEner
               << " "
               << _totalIntElecEner
               << " "
               << _totalIntIonsEner
	       << " "
               << _totalTotalEMEner
               << " "
               << _totalTotalElecEner
               << " "
               << _totalTotalIonsEner
               << " "
               << _totalTotalEner
               << " "
               << CurrentTime
               << " "
               << _totalTotalNbCells
               << "\n";
      outputFile.close();
    }
    else {
      ofstream& outputFile = (*fhandle)->open(constructFilename(), ios::app);
      outputFile << iter
                 << " "
                 << _totalMagneticEner
                 << " "
                 << _totalElectricEner
                 << " "
                 << _totalPsiEner
                 << " "
                 << _totalPhiEner
                 << " "
                 << _totalKinElecEner
                 << " "
                 << _totalKinIonsEner
                 << " "
                 << _totalIntElecEner
                 << " "
                 << _totalIntIonsEner
                 << " "
                 << _totalTotalEMEner
                 << " "
                 << _totalTotalElecEner
                 << " "
                 << _totalTotalIonsEner
                 << " "
                 << _totalTotalEner
                 << " "
                 << CurrentTime
                 << " "
                 << _totalTotalNbCells
                 << "\n";
      outputFile.close();
    }
    delete fhandle;
  }
  CFout << "Writing of Energy integration" << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void EnergyIntegration::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

