#include "Common/PE.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "MathTools/MathConsts.hh"
#include "FiniteVolumeMaxwell/DivMonitoring.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Common/BadValueException.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"
#include "FiniteVolumeMaxwell/FiniteVolumeMaxwell.hh"
#include "Maxwell/MaxwellVarSet.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "TecplotWriter/WriteTecplot.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Maxwell;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::TecplotWriter;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeMaxwell {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<DivMonitoring, DataProcessingData, FiniteVolumeMaxwellModule>
divMonitoringProvider("DivMonitoring");

//////////////////////////////////////////////////////////////////////////////

void DivMonitoring::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("SaveRate","Save Output File every...iterations.");
   options.addConfigOption< bool >("AppendTime","Save each iteration to different file with suffix _time#.");
   options.addConfigOption< bool >("AppendIter","Save each iteration to different file with suffix _iter#.");
   options.addConfigOption< std::string >("OutputFileDivMonitoring","Name of Output File to write the electric and magnetic field divergence");
   options.addConfigOption< std::string >("Command","something to run...");
}

//////////////////////////////////////////////////////////////////////////////
      
DivMonitoring::DivMonitoring(const std::string& name) :
  DataProcessingCom(name),
  socket_avgBxFace("avgBxFace"),
  socket_avgByFace("avgByFace"),
  socket_avgBzFace("avgBzFace"),
  socket_avgExFace("avgExFace"),
  socket_avgEyFace("avgEyFace"),
  socket_avgEzFace("avgEzFace"),  
  socket_divB("divB"),
  socket_curlB("curlB"),
  socket_divE("divE"),
  socket_theta("theta"),
  socket_Bradial("Bradial"),
  socket_Btheta("Btheta"),
  socket_BthetaTheory("BthetaTheory"),
  socket_ExPWTheory("ExPWTheory"),
  socket_EyPWTheory("EyPWTheory"),
  socket_EzPWTheory("EzPWTheory"),
  socket_BxPWTheory("BxPWTheory"),
  socket_ByPWTheory("ByPWTheory"),
  socket_BzPWTheory("BzPWTheory"),
  socket_Error("Error"), 
  socket_volumes("volumes"),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_isOutward("isOutward"),
  socket_normals("normals"),
  socket_ErrorExPW("ErrorExPW"),
  socket_ErrorEyPW("ErrorEyPW"),
  socket_ErrorEzPW("ErrorEzPW"),
  socket_ErrorBxPW("ErrorBxPW"),
  socket_ErrorByPW("ErrorByPW"),
  socket_ErrorBzPW("ErrorBzPW"),
  m_geoBuilder(),
  m_faceBuilder(),  
  m_updateVarSet(CFNULL),
  m_divBL2Norm(CFNULL),
  m_divEL2Norm(CFNULL),
  m_TEPWErrorL2Norm(CFNULL),
  m_TMPWErrorL2Norm(CFNULL),
  m_fullPWErrorL2Norm(CFNULL)

{
  addConfigOptionsTo(this);
 
   m_saveRate = 1;
   setParameter("SaveRate",&m_saveRate);

   m_appendIter = false;
   setParameter("AppendIter",&m_appendIter);

   m_appendTime = false;
   setParameter("AppendTime",&m_appendTime);

   m_nameOutputFileDivMonitoring = "DivMonitoring.plt";
   setParameter("OutputFileDivMonitoring",&m_nameOutputFileDivMonitoring);

   m_toRun = "";
   setParameter("Command",&m_toRun);
}
      
//////////////////////////////////////////////////////////////////////////////

DivMonitoring::~DivMonitoring()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> > DivMonitoring::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  
  result.push_back(&socket_avgBxFace);
  result.push_back(&socket_avgByFace);
  result.push_back(&socket_avgBzFace);
  result.push_back(&socket_avgExFace);
  result.push_back(&socket_avgEyFace);
  result.push_back(&socket_avgEzFace);  
  result.push_back(&socket_divB);
  result.push_back(&socket_curlB);
  result.push_back(&socket_divE);
  result.push_back(&socket_theta);
  result.push_back(&socket_Bradial);
  result.push_back(&socket_Btheta);
  result.push_back(&socket_BthetaTheory);  
  result.push_back(&socket_ExPWTheory);  
  result.push_back(&socket_EyPWTheory);  
  result.push_back(&socket_EzPWTheory);  
  result.push_back(&socket_BxPWTheory);  
  result.push_back(&socket_ByPWTheory);  
  result.push_back(&socket_BzPWTheory);
  result.push_back(&socket_Error); 
  result.push_back(&socket_ErrorExPW);
  result.push_back(&socket_ErrorEyPW);
  result.push_back(&socket_ErrorEzPW);
  result.push_back(&socket_ErrorBxPW);
  result.push_back(&socket_ErrorByPW);
  result.push_back(&socket_ErrorBzPW);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
DivMonitoring::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_volumes);
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_nstates);
  result.push_back(&socket_gstates);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_normals);
  
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////

void DivMonitoring::setup()
{
  CFAUTOTRACE;
  
  // run WriteTecplot setup
  WriteTecplot::getInstance().setup();
  WriteTecplot::getInstance().setDataSockets(socket_volumes, socket_nodes);
  
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
   
  // geometry builder setup
  m_geoBuilder.setup();
  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  m_geoBuilder.getDataGE().trs = cells;
  
  // geometry builder setup
  m_faceBuilder.setup();
  m_faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
//   m_faceBuilder.getDataGE().trs = faces;  
  
  m_updateVarSet = getMethodData().getUpdateVarSet().d_castTo<MaxwellVarSet>();
  
  if (m_toRun.size() != 0)
  {
     m_toRun = m_toRun + " & ";
     system(m_toRun.c_str());
  }  
    
  DataHandle<CFreal> avgBxFace  = socket_avgBxFace.getDataHandle();
  avgBxFace.resize(nbFaces);
    
  DataHandle<CFreal> avgByFace  = socket_avgByFace.getDataHandle();
  avgByFace.resize(nbFaces);
  
  DataHandle<CFreal> avgExFace  = socket_avgExFace.getDataHandle();
  avgExFace.resize(nbFaces);
    
  DataHandle<CFreal> avgEyFace  = socket_avgEyFace.getDataHandle();
  avgEyFace.resize(nbFaces);  
  
  DataHandle<CFreal> avgBzFace  = socket_avgBzFace.getDataHandle();
  avgBzFace.resize(nbFaces);
    
  DataHandle<CFreal> avgEzFace  = socket_avgEzFace.getDataHandle();
  avgEzFace.resize(nbFaces);    
    
  // preparing divB socket
  DataHandle<CFreal> divB = socket_divB.getDataHandle();
  divB.resize(nbCells);
  
  // preparing curlB socket
  DataHandle<CFreal> curlB = socket_curlB.getDataHandle();
  curlB.resize(nbCells);  
  
  // preparing divE socket
  DataHandle<CFreal> divE = socket_divE.getDataHandle();
  divE.resize(nbCells);
  
  // preparing theta socket
  DataHandle<CFreal> theta = socket_theta.getDataHandle();
  theta.resize(nbCells);  
  
  DataHandle<CFreal> Bradial = socket_Bradial.getDataHandle();
  Bradial.resize(nbCells); 
  
  DataHandle<CFreal> Btheta = socket_Btheta.getDataHandle();
  Btheta.resize(nbCells);
  
  DataHandle<CFreal> BthetaTheory = socket_BthetaTheory.getDataHandle();
  BthetaTheory.resize(nbCells);

  DataHandle<CFreal> ExPWTheory = socket_ExPWTheory.getDataHandle();
  ExPWTheory.resize(nbCells);

  DataHandle<CFreal> EyPWTheory = socket_EyPWTheory.getDataHandle();
  EyPWTheory.resize(nbCells);

  DataHandle<CFreal> EzPWTheory = socket_EzPWTheory.getDataHandle();
  EzPWTheory.resize(nbCells);

  DataHandle<CFreal> BxPWTheory = socket_BxPWTheory.getDataHandle();
  BxPWTheory.resize(nbCells);

  DataHandle<CFreal> ByPWTheory = socket_ByPWTheory.getDataHandle();
  ByPWTheory.resize(nbCells);

  DataHandle<CFreal> BzPWTheory = socket_BzPWTheory.getDataHandle();
  BzPWTheory.resize(nbCells);

  DataHandle<CFreal> Error = socket_Error.getDataHandle();
  Error.resize(nbCells);  
  
  // preparing Error in Planar wave testcase
  DataHandle<CFreal> ErrorExPW = socket_ErrorExPW.getDataHandle();
  ErrorExPW.resize(nbCells);  

  DataHandle<CFreal> ErrorEyPW = socket_ErrorEyPW.getDataHandle();
  ErrorEyPW.resize(nbCells);  

  DataHandle<CFreal> ErrorEzPW = socket_ErrorEzPW.getDataHandle();
  ErrorEzPW.resize(nbCells);  
  
  DataHandle<CFreal> ErrorBxPW = socket_ErrorBxPW.getDataHandle();
  ErrorBxPW.resize(nbCells);  

  DataHandle<CFreal> ErrorByPW = socket_ErrorByPW.getDataHandle();
  ErrorByPW.resize(nbCells);  

  DataHandle<CFreal> ErrorBzPW = socket_ErrorBzPW.getDataHandle();
  ErrorBzPW.resize(nbCells);  
  
  for (CFuint iCell = 0; iCell < nbCells; iCell++){
    divB[iCell] = 0.0;
    curlB[iCell] = 0.0;
    divE[iCell] = 0.0;
    theta[iCell] = 0.0;
    Bradial[iCell] = 0.0;
    Btheta[iCell] = 0.0;
    BthetaTheory[iCell] = 0.0;    
    ExPWTheory[iCell] = 0.0; 
    EyPWTheory[iCell] = 0.0; 
    EzPWTheory[iCell] = 0.0; 
    BxPWTheory[iCell] = 0.0; 
    ByPWTheory[iCell] = 0.0; 
    BzPWTheory[iCell] = 0.0;
    Error[iCell] = 0.0;    
    ErrorExPW[iCell] = 0.0;
    ErrorEyPW[iCell] = 0.0;
    ErrorEzPW[iCell] = 0.0;
    ErrorBxPW[iCell] = 0.0;
    ErrorByPW[iCell] = 0.0;
    ErrorBzPW[iCell] = 0.0;
//    divBL2Norm[iCell] = 0.0;
//    divEL2Norm[iCell] = 0.0;
//    TEPWErrorL2Norm[iCell] = 0.0;    
  }
}
//////////////////////////////////////////////////////////////////////////////

void DivMonitoring::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DivMonitoring::execute()
{
  CFAUTOTRACE;
  
  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> divB = socket_divB.getDataHandle();
  DataHandle<CFreal> curlB = socket_curlB.getDataHandle();
  DataHandle<CFreal> divE = socket_divE.getDataHandle();  
  DataHandle<CFreal> theta = socket_theta.getDataHandle();
  DataHandle<CFreal> Bradial = socket_Bradial.getDataHandle();  
  DataHandle<CFreal> Btheta = socket_Btheta.getDataHandle();
  DataHandle<CFreal> BthetaTheory = socket_BthetaTheory.getDataHandle();  
  DataHandle<CFreal> ExPWTheory = socket_ExPWTheory.getDataHandle();  
  DataHandle<CFreal> EyPWTheory = socket_EyPWTheory.getDataHandle();  
  DataHandle<CFreal> EzPWTheory = socket_EzPWTheory.getDataHandle();  
  DataHandle<CFreal> BxPWTheory = socket_BxPWTheory.getDataHandle();  
  DataHandle<CFreal> ByPWTheory = socket_ByPWTheory.getDataHandle();  
  DataHandle<CFreal> BzPWTheory = socket_BzPWTheory.getDataHandle(); 
  DataHandle<CFreal> Error = socket_Error.getDataHandle();  
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();  
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  DataHandle<CFreal> ErrorExPW = socket_ErrorExPW.getDataHandle();  
  DataHandle<CFreal> ErrorEyPW = socket_ErrorEyPW.getDataHandle();  
  DataHandle<CFreal> ErrorEzPW = socket_ErrorEzPW.getDataHandle();
  DataHandle<CFreal> ErrorBxPW = socket_ErrorBxPW.getDataHandle();  
  DataHandle<CFreal> ErrorByPW = socket_ErrorByPW.getDataHandle();  
  DataHandle<CFreal> ErrorBzPW = socket_ErrorBzPW.getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
//   std::cout << "DivMonitoring::execute => before computing \n";
  
  const CFreal pi = 3.14159265358979323846; 
  const CFreal c_e = 299792458;

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");  
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;  
    
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  //CFLog( VERBOSE, "nbCells" << nbCells << "\n");  

  m_divBL2Norm.resize(1);
  m_divBL2Norm = 0.;
  m_divEL2Norm.resize(1);
  m_divEL2Norm = 0.;
  m_TEPWErrorL2Norm.resize(1);
  m_TEPWErrorL2Norm = 0.;
  m_TMPWErrorL2Norm.resize(1);
  m_TMPWErrorL2Norm = 0.;
  m_fullPWErrorL2Norm.resize(1);
  m_fullPWErrorL2Norm = 0.;  
 
  computeDivergence ();
  computeAnaliticalSolution();
  computeCylindricalComponents();
  
//   std::cout << "DivMonitoring::execute => after computing \n";  
  
  ///Computation of Errors
  CFreal TESumTheoryNorm2 = 0;		//Summation of the L2 Norm of theoretical solution of TE Planar wave
  CFreal TMSumTheoryNorm2 = 0;		//Summation of the L2 Norm of theoretical solution of TM Planar wave
  CFreal fullSumTheoryNorm2 = 0;	//Summation of the L2 Norm of theoretical solution of Full Planar wave 
  CFreal TESumDifferenceNorm2 = 0;	//Summation of the L2 Norm of difference between theoretical and numerical solution of TE Planar wave  
  CFreal TMSumDifferenceNorm2 = 0;	//Summation of the L2 Norm of difference between theoretical and numerical solution of TM Planar wave  
  CFreal fullSumDifferenceNorm2 = 0;	//Summation of the L2 Norm of difference between theoretical and numerical solution of Full Planar wave  
  CFreal TotalVolume = 0;		//Total Volume of the domain
  CFreal SumDivE2 = 0;
  CFreal SumDivB2 = 0;

  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    
    // set the builder data and build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();
    const CFuint elemID = currCell->getID(); 
    
    //Error Wire Testcase
    CFreal ErrorNonAbs = 0;
    //Error Coaxial or Planar wave testcase
    CFreal ErrorExPWNonAbs = 0;
    CFreal ErrorEyPWNonAbs = 0;    
    CFreal ErrorEzPWNonAbs = 0; 
    CFreal ErrorBxPWNonAbs = 0;
    CFreal ErrorByPWNonAbs = 0;    
    CFreal ErrorBzPWNonAbs = 0;     
   
    CFreal Bx = 0;
    CFreal By = 0;    
    CFreal Bz = 0; 
    CFreal Ex = 0;   
    CFreal Ey = 0;
    CFreal Ez = 0;
 
    //Set the states
    State *currState = currCell->getState(0);
    Bx = (*currState)[0];
    By = (*currState)[1];
    Bz = (*currState)[2];
    Ex = (*currState)[3];    
    Ey = (*currState)[4];     
    Ez = (*currState)[5]; 

    ///Error Computation Every cell in Wire testcase
    ErrorNonAbs = ((BthetaTheory[iCell] - Btheta[iCell])/BthetaTheory[iCell])*100;
    Error[iCell] = std::abs(ErrorNonAbs);
    
    ///Error Computation every cell in Coaxial or Plane wave testcases
    ErrorExPWNonAbs = (ExPWTheory[iCell] - Ex);
    ErrorExPW[iCell] = std::abs(100*ErrorExPWNonAbs/ExPWTheory[iCell]);    

    ErrorEyPWNonAbs = (EyPWTheory[iCell] - Ey);
    ErrorEyPW[iCell] = std::abs(100*ErrorEyPWNonAbs/EyPWTheory[iCell]);  

    ErrorEzPWNonAbs = (EzPWTheory[iCell] - Ez);
    ErrorEzPW[iCell] = std::abs(100*ErrorEzPWNonAbs/EzPWTheory[iCell]);  

    ErrorBxPWNonAbs = (BxPWTheory[iCell] - Bx);
    ErrorBxPW[iCell] = std::abs(100*ErrorBxPWNonAbs/BxPWTheory[iCell]);    

    ErrorByPWNonAbs = (ByPWTheory[iCell] - By);
    ErrorByPW[iCell] = std::abs(100*ErrorByPWNonAbs/ByPWTheory[iCell]);  

    ErrorBzPWNonAbs = (BzPWTheory[iCell] - Bz);
    ErrorBzPW[iCell] = std::abs(100*ErrorBzPWNonAbs/BzPWTheory[iCell]);      


///ADimensional case
// 
//     TESumTheoryNorm2 += (ExPWTheory[iCell]*ExPWTheory[iCell] + EyPWTheory[iCell]*EyPWTheory[iCell] + BzPWTheory[iCell]*BzPWTheory[iCell])*volumes[elemID];
//     TESumDifferenceNorm2 += (ErrorExPWNonAbs*ErrorExPWNonAbs + ErrorEyPWNonAbs*ErrorEyPWNonAbs + ErrorBzPWNonAbs*ErrorBzPWNonAbs)*volumes[elemID];


//     TMSumTheoryNorm2 += (BxPWTheory[iCell]*BxPWTheory[iCell] + ByPWTheory[iCell]*ByPWTheory[iCell] + EzPWTheory[iCell]*EzPWTheory[iCell])*volumes[elemID];
//     TMSumDifferenceNorm2 += (ErrorBxPWNonAbs*ErrorBxPWNonAbs + ErrorByPWNonAbs*ErrorByPWNonAbs + ErrorEzPWNonAbs*ErrorEzPWNonAbs)*volumes[elemID];
    
//     fullSumTheoryNorm2 += (BxPWTheory[iCell]*BxPWTheory[iCell] + ByPWTheory[iCell]*ByPWTheory[iCell] + BzPWTheory[iCell]*BzPWTheory[iCell] + ExPWTheory[iCell]*ExPWTheory[iCell] + EyPWTheory[iCell]*EyPWTheory[iCell] + EzPWTheory[iCell]*EzPWTheory[iCell])*volumes[elemID];
//     fullSumDifferenceNorm2 += (ErrorBxPWNonAbs*ErrorBxPWNonAbs + ErrorByPWNonAbs*ErrorByPWNonAbs + ErrorBzPWNonAbs*ErrorBzPWNonAbs + ErrorExPWNonAbs*ErrorExPWNonAbs + ErrorEyPWNonAbs*ErrorEyPWNonAbs + ErrorEzPWNonAbs*ErrorEzPWNonAbs)*volumes[elemID];    
 
///Dimensional case

    TESumTheoryNorm2 += (ExPWTheory[iCell]*ExPWTheory[iCell] + EyPWTheory[iCell]*EyPWTheory[iCell] + BzPWTheory[iCell]*BzPWTheory[iCell]*c_e*c_e)*volumes[elemID];
    TESumDifferenceNorm2 += (ErrorExPWNonAbs*ErrorExPWNonAbs + ErrorEyPWNonAbs*ErrorEyPWNonAbs + ErrorBzPWNonAbs*ErrorBzPWNonAbs*c_e*c_e)*volumes[elemID];


    TMSumTheoryNorm2 += (BxPWTheory[iCell]*BxPWTheory[iCell]*c_e*c_e + ByPWTheory[iCell]*ByPWTheory[iCell]*c_e*c_e + EzPWTheory[iCell]*EzPWTheory[iCell])*volumes[elemID];
    TMSumDifferenceNorm2 += (ErrorBxPWNonAbs*ErrorBxPWNonAbs*c_e*c_e + ErrorByPWNonAbs*ErrorByPWNonAbs*c_e*c_e + ErrorEzPWNonAbs*ErrorEzPWNonAbs)*volumes[elemID];
    
    fullSumTheoryNorm2 += (BxPWTheory[iCell]*BxPWTheory[iCell]*c_e*c_e + ByPWTheory[iCell]*ByPWTheory[iCell]*c_e*c_e + BzPWTheory[iCell]*BzPWTheory[iCell]*c_e*c_e + ExPWTheory[iCell]*ExPWTheory[iCell] + EyPWTheory[iCell]*EyPWTheory[iCell] + EzPWTheory[iCell]*EzPWTheory[iCell])*volumes[elemID];
    fullSumDifferenceNorm2 += (ErrorBxPWNonAbs*ErrorBxPWNonAbs*c_e*c_e + ErrorByPWNonAbs*ErrorByPWNonAbs*c_e*c_e + ErrorBzPWNonAbs*ErrorBzPWNonAbs*c_e*c_e + ErrorExPWNonAbs*ErrorExPWNonAbs + ErrorEyPWNonAbs*ErrorEyPWNonAbs + ErrorEzPWNonAbs*ErrorEzPWNonAbs)*volumes[elemID];    
    
    
    SumDivB2 += divB[iCell]*divB[iCell]*volumes[elemID];
    SumDivE2 += divE[iCell]*divE[iCell]*volumes[elemID];
    TotalVolume += volumes[elemID];
         
    m_geoBuilder.releaseGE();
    
  }
  //std::cout <<"after the cycle \n";
    //CFLog( VERBOSE, "after the cycle" << "\n");
    
    
    
    
  //Total sum of the processors  
  CFdouble totalTESumTheoryNorm2 = 0.0;
  CFdouble totalTESumDifferenceNorm2 = 0.0;
  CFdouble totalTMSumTheoryNorm2 = 0.0;
  CFdouble totalTMSumDifferenceNorm2 = 0.0; 
  CFdouble totalfullSumTheoryNorm = 0.0; 
  CFdouble totalfullSumDifferenceNorm2 = 0.0;  
  CFdouble totalSumDivB2 = 0.0;  
  CFdouble totalSumDivE2 = 0.0;  
  CFdouble totalTotalVolume = 0.0;  
  
  const std::string nsp = this->getMethodData().getNamespace();
  
  if (PE::GetPE().GetProcessorCount(nsp) > 1) {

#ifdef CF_HAVE_MPI
	MPI_Allreduce(&TESumTheoryNorm2, &totalTESumTheoryNorm2, 1, MPI_DOUBLE, MPI_SUM,
		      PE::GetPE().GetCommunicator(nsp));
	MPI_Allreduce(&TESumDifferenceNorm2, &totalTESumDifferenceNorm2, 1, MPI_DOUBLE, MPI_SUM,
		      PE::GetPE().GetCommunicator(nsp));
	MPI_Allreduce(&TMSumTheoryNorm2, &totalTMSumTheoryNorm2, 1, MPI_DOUBLE, MPI_SUM,
		      PE::GetPE().GetCommunicator(nsp));
	MPI_Allreduce(&TMSumDifferenceNorm2, &totalTMSumDifferenceNorm2, 1, MPI_DOUBLE, MPI_SUM,
		      PE::GetPE().GetCommunicator(nsp));
	MPI_Allreduce(&fullSumTheoryNorm2, &totalfullSumTheoryNorm, 1, MPI_DOUBLE, MPI_SUM,
		      PE::GetPE().GetCommunicator(nsp));
	MPI_Allreduce(&fullSumDifferenceNorm2, &totalfullSumDifferenceNorm2, 1, MPI_DOUBLE, MPI_SUM,
		      PE::GetPE().GetCommunicator(nsp));
	MPI_Allreduce(&SumDivB2, &totalSumDivB2, 1, MPI_DOUBLE, MPI_SUM,
		      PE::GetPE().GetCommunicator(nsp));
	MPI_Allreduce(&SumDivE2, &totalSumDivE2, 1, MPI_DOUBLE, MPI_SUM,
		      PE::GetPE().GetCommunicator(nsp));
	MPI_Allreduce(&TotalVolume, &totalTotalVolume, 1, MPI_DOUBLE, MPI_SUM,
		      PE::GetPE().GetCommunicator(nsp));
	
	TESumTheoryNorm2 = totalTESumTheoryNorm2;
	TESumDifferenceNorm2 = totalTESumDifferenceNorm2;
	TMSumTheoryNorm2 = totalTMSumTheoryNorm2;
	TMSumDifferenceNorm2 = totalTMSumDifferenceNorm2;
	fullSumTheoryNorm2 = totalfullSumTheoryNorm;
	fullSumDifferenceNorm2 = totalfullSumDifferenceNorm2;
	SumDivB2 = totalSumDivB2;
	SumDivE2 = totalSumDivE2; 
	TotalVolume = totalTotalVolume;
    
#endif
  }
     
  CFreal TEPWErrorL2Norm =  std::sqrt(TESumDifferenceNorm2/TESumTheoryNorm2);
  CFreal TMPWErrorL2Norm =  std::sqrt(TMSumDifferenceNorm2/TMSumTheoryNorm2);
  CFreal fullPWErrorL2Norm = std::sqrt(fullSumDifferenceNorm2/fullSumTheoryNorm2);
  CFreal divBL2Norm = std::sqrt(SumDivB2/TotalVolume);
  CFreal divEL2Norm = std::sqrt(SumDivE2/TotalVolume);      
    
    
  if(!(iter % m_saveRate)) {
    m_divBL2Norm = divBL2Norm;
    m_divEL2Norm = divEL2Norm;
    m_TEPWErrorL2Norm = TEPWErrorL2Norm*100;
    m_TMPWErrorL2Norm = TMPWErrorL2Norm*100;
    m_fullPWErrorL2Norm = fullPWErrorL2Norm*100;

    WriteTecplot::getInstance().setNodeExtrapolation(1); 
      //we can use default 1
    writeOutputFile();
  }  
  //std::cout << "OUTexectue \n";
}
//////////////////////////////////////////////////////////////////////////////
void DivMonitoring::computeDivergence()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> avgBxFace  = socket_avgBxFace.getDataHandle();
  DataHandle<CFreal> avgByFace  = socket_avgByFace.getDataHandle();
  DataHandle<CFreal> avgBzFace  = socket_avgBzFace.getDataHandle();
  DataHandle<CFreal> avgExFace  = socket_avgExFace.getDataHandle();
  DataHandle<CFreal> avgEyFace  = socket_avgEyFace.getDataHandle();
  DataHandle<CFreal> avgEzFace  = socket_avgEzFace.getDataHandle();  
  DataHandle<CFreal> divB = socket_divB.getDataHandle();
  DataHandle<CFreal> curlB = socket_curlB.getDataHandle();
  DataHandle<CFreal> divE = socket_divE.getDataHandle();  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();  
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  
       
  CFuint faceIdx = 0;
  
  CFreal faceAvBx = 0;
  CFreal faceAvBy = 0;
  CFreal faceAvBz = 0;
  CFreal faceAvEx = 0;
  CFreal faceAvEy = 0;
  CFreal faceAvEz = 0;

  vector<SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trs.size();
  FaceTrsGeoBuilder::GeoData& faceData = m_faceBuilder.getDataGE();
  
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    
    SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];
   
    if (currTrs->getName() != "PartitionFaces" && currTrs->getName() != "InnerCells") {    
       if (currTrs->hasTag("writable")) {
       }
       else {
	faceData.trs = currTrs; 
	const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
	for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace, ++faceIdx) {
	  faceData.idx = iFace;
	  GeometricEntity *const currFace = m_faceBuilder.buildGE();  

	  if (currFace->getState(0)->isParUpdatable() || 
	    (!currFace->getState(1)->isGhost() && currFace->getState(1)->isParUpdatable())) {
	    const CFuint currFaceID = currFace->getID();	    
	    avgBxFace[currFaceID] = 0;
	    avgByFace[currFaceID] = 0;
	    avgBzFace[currFaceID] = 0;
	    avgExFace[currFaceID] = 0;
	    avgEyFace[currFaceID] = 0;
	    avgEzFace[currFaceID] = 0; 
	    State *lState = currFace->getState(0);
	    State *rState = currFace->getState(1);       
	    faceAvBx = 0.5*((*rState)[0] + (*lState)[0]);
	    faceAvBy = 0.5*((*rState)[1] + (*lState)[1]);
	    faceAvBz = 0.5*((*rState)[2] + (*lState)[2]);
	    faceAvEx = 0.5*((*rState)[3] + (*lState)[3]);
	    faceAvEy = 0.5*((*rState)[4] + (*lState)[4]);
	    faceAvEz = 0.5*((*rState)[5] + (*lState)[5]);

	    avgBxFace[currFaceID] = faceAvBx;
	    avgByFace[currFaceID] = faceAvBy;
	    avgBzFace[currFaceID] = faceAvBz;
	    avgExFace[currFaceID] = faceAvEx;
	    avgEyFace[currFaceID] = faceAvEy;
	    avgEzFace[currFaceID] = faceAvEz; 
	  }
	  m_faceBuilder.releaseGE();      
	}
      }
    }
        
    if (currTrs->getName() == "PartitionFaces") {  
      faceData.trs = currTrs; 
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace, ++faceIdx) {   
        faceData.idx = iFace;       
        GeometricEntity* currFace = m_faceBuilder.buildGE();       
        const CFuint currFaceID = currFace->getID();      
        avgBxFace[currFaceID] = 0;
        avgByFace[currFaceID] = 0;
        avgBzFace[currFaceID] = 0;
        avgExFace[currFaceID] = 0;
        avgEyFace[currFaceID] = 0;
        avgEzFace[currFaceID] = 0; 
        m_faceBuilder.releaseGE();      
      }
    }   
  }
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
    
  Common::SafePtr<TopologicalRegionSet> cells =
  MeshDataStack::getActive()->getTrs("InnerCells");
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;    
    
  const CFuint nbCells = cells->getLocalNbGeoEnts(); 
     
  for (CFuint iCell = 0; iCell < nbCells; ++iCell){
    divB[iCell] = 0;
    curlB[iCell] = 0;
    divE[iCell] = 0;
    CFuint atPartitionCell = 0;
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();
    const CFuint elemID = currCell->getID(); 
    const vector<GeometricEntity*>& faces = *currCell->getNeighborGeos();
    const CFuint nbFaces = faces.size();
//       std::vector<Node*> localNodes = *currCell->getNodes();
//       const CFuint nbNodesInElem = localNodes.size();
    
    for (CFuint i = 0; i < nbFaces; ++i){
	
      const CFuint faceID = faces[i]->getID();
      const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
      CFreal nx = normals[startID];
      CFreal ny = normals[startID + 1];
      //CFLog( VERBOSE, "nx = " << nx << " ny = " << ny <<" norm = " << sqrt(nx*nx + ny*ny) << "\n");     

      faceAvBx = avgBxFace[faceID];
      faceAvBy = avgByFace[faceID];
      faceAvBz = avgBzFace[faceID];
      faceAvEx = avgExFace[faceID];
      faceAvEy = avgEyFace[faceID];
      faceAvEz = avgEzFace[faceID];      
	
      if (faceAvBx == 0 && faceAvBy == 0 && faceAvBz == 0 && 
	faceAvEx == 0 && faceAvEy == 0 && faceAvEz == 0 ){
	  atPartitionCell = 1;
      }
	
      if (static_cast<CFuint>(isOutward[faceID]) != elemID) {
	nx *= -1.;
	ny *= -1.;
      }
      
      CFreal faceArea = std::sqrt(nx*nx + ny*ny);
      curlB[iCell] += -ny*faceAvBx + nx*faceAvBy;
      divB[iCell] += nx*faceAvBx + ny*faceAvBy;
      divE[iCell] += nx*faceAvEx + ny*faceAvEy;
      
      if(dim == 3){
	
	CFreal nz = normals[startID + 2];
	
	if (static_cast<CFuint>(isOutward[faceID]) != elemID) {
	  nz *= -1.;
	}
	
	divB[iCell] += nz*faceAvBz;
	divE[iCell] += nz*faceAvEz;
      }
    }
    
    curlB[iCell] *= 1/volumes[elemID];
    divB[iCell] *= 1/volumes[elemID];
    divE[iCell] *= 1/volumes[elemID]; 
  
    if (atPartitionCell == 1) {
      curlB[iCell] = 0;
      divB[iCell] = 0;
      divE[iCell] = 0;
    } 
    
    m_geoBuilder.releaseGE();
  }  
  
  
/*  else {
    //std::cout <<"Dim = 3 chosen \n";

    CFreal faceAvBx = 0;
    CFreal faceAvBy = 0;
    CFreal faceAvBz = 0;
    CFreal faceAvEx = 0;
    CFreal faceAvEy = 0;
    CFreal faceAvEz = 0;
    
    // set the list of faces
    
    vector<SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
    const CFuint nbTRSs = trs.size();    
    
    Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");
    CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
    geoData.trs = cells;    
    
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    
    for (CFuint iCell = 0; iCell < nbCells; ++iCell){
      
      divB[iCell] = 0;
      divE[iCell] = 0;
      
//       //std::cout <<"In Element to compute the divergence 3D \n";

      
      // set the builder data and build the GeometricEntity
      geoData.idx = iCell;
      GeometricEntity* currCell = m_geoBuilder.buildGE();
      const CFuint elemID = currCell->getID(); 
      const vector<GeometricEntity*>& faces = *currCell->getNeighborGeos();
      const CFuint nbFaces = faces.size();
    
      std::vector<Node*> localNodes = *currCell->getNodes();
      const CFuint nbNodesInElem = localNodes.size();
      
//       //std::cout <<"In Element to compute the divergence \n";
//       //std::cout <<"Nb of Nodes in Element= \t "<< nbNodesInElem <<"\n";
//       //std::cout <<"Nb of Faces in Element= \t "<< nbFaces <<"\n";
      
      for (CFuint i = 0; i < nbFaces; ++i){
	
// 	//std::cout <<"In face to compute the divergence 3D \n";
	
	State *lState = faces[i]->getState(0);
	State *rState = faces[i]->getState(1);
//  	//std::cout <<"At the face to compute the face average \n";
	
	const CFuint faceID = faces[i]->getID();
	const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
	CFreal nx = normals[startID];
	CFreal ny = normals[startID + 1];
	CFreal nz = normals[startID + 2];
// 	CFLog( VERBOSE, "nx = " << nx << " ny = " << ny <<" norm = " << sqrt(nx*nx + ny*ny) << "\n");     
      
	if (static_cast<CFuint>(isOutward[faceID]) != elemID) {
	  nx *= -1.;
	  ny *= -1.;
	  nz *= -1.;
	}
//  	//std::cout <<"After taking the values of the cell and the neighbour\n";
	
	faceAvBx = 0.5*((*rState)[0] + (*lState)[0]);
	faceAvBy = 0.5*((*rState)[1] + (*lState)[1]);
	faceAvBz = 0.5*((*rState)[2] + (*lState)[2]);
	faceAvEx = 0.5*((*rState)[3] + (*lState)[3]);
	faceAvEy = 0.5*((*rState)[4] + (*lState)[4]);
	faceAvEz = 0.5*((*rState)[5] + (*lState)[5]);
	
	divB[iCell] += nx*faceAvBx + ny*faceAvBy + nz*faceAvBz;
	divE[iCell] += nx*faceAvEx + ny*faceAvEy + nz*faceAvEz;
      }
      divB[iCell] *= 1/volumes[elemID];
      divE[iCell] *= 1/volumes[elemID]; 
//        //std::cout <<"After calculating div*[iCell]\n";
      m_geoBuilder.releaseGE();
    }
  }   */
      
//   const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
//   if (dim == 2){
//     
//     CFreal faceAvBx = 0;
//     CFreal faceAvBy = 0;
//     CFreal faceAvBz = 0;
//     CFreal faceAvEx = 0;
//     CFreal faceAvEy = 0;
//     CFreal faceAvEz = 0;
//     
//     Common::SafePtr<TopologicalRegionSet> cells =
//     MeshDataStack::getActive()->getTrs("InnerCells");
//     CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
//     geoData.trs = cells;    
//     
//     const CFuint nbCells = cells->getLocalNbGeoEnts();   
//     
//     for (CFuint iCell = 0; iCell < nbCells; ++iCell){
//       
//       divB[iCell] = 0;
//       divE[iCell] = 0;
//       
//       // set the builder data and build the GeometricEntity
//       geoData.idx = iCell;
//       GeometricEntity* currCell = m_geoBuilder.buildGE();
//       const CFuint elemID = currCell->getID(); 
//       const vector<GeometricEntity*>& faces = *currCell->getNeighborGeos();
//       std::vector<Node*> localNodes = *currCell->getNodes();
//       const CFuint nbNodesInElem = localNodes.size();
//       
// //       //std::cout <<"In Element to compute the divergence \n";
//       
//       for (CFuint i = 0; i < nbNodesInElem; ++i){
// 	
// 	State *lState = faces[i]->getState(0);
// 	State *rState = faces[i]->getState(1);
// // 	//std::cout <<"At the face to compute the face average \n";
// 	
// 	const CFuint faceID = faces[i]->getID();
// 	const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
// 	CFreal nx = normals[startID];
// 	CFreal ny = normals[startID + 1];
// 	//CFLog( VERBOSE, "nx = " << nx << " ny = " << ny <<" norm = " << sqrt(nx*nx + ny*ny) << "\n");     
//       
// 	if (static_cast<CFuint>(isOutward[faceID]) != elemID) {
// 	  nx *= -1.;
// 	  ny *= -1.;
// 	}
// // 	//std::cout <<"After taking the values of the cell and the neighbour\n";
// 	
// 	faceAvBx = 0.5*((*rState)[0] + (*lState)[0]);
// 	faceAvBy = 0.5*((*rState)[1] + (*lState)[1]);
// 	faceAvBz = 0.5*((*rState)[2] + (*lState)[2]);
// 	faceAvEx = 0.5*((*rState)[3] + (*lState)[3]);
// 	faceAvEy = 0.5*((*rState)[4] + (*lState)[4]);
// 	faceAvEz = 0.5*((*rState)[5] + (*lState)[5]);
// 	
// 	divB[iCell] += nx*faceAvBx + ny*faceAvBy;
// 	divE[iCell] += nx*faceAvEx + ny*faceAvEy;
//       }
//       divB[iCell] *= 1/volumes[elemID];
//       divE[iCell] *= 1/volumes[elemID]; 
// //       //std::cout <<"After calculating div*[iCell]\n";
//       m_geoBuilder.releaseGE();
//     }
//   }  
//   
//   
//   else {
//     //std::cout <<"Dim = 3 chosen \n";
// 
//     CFreal faceAvBx = 0;
//     CFreal faceAvBy = 0;
//     CFreal faceAvBz = 0;
//     CFreal faceAvEx = 0;
//     CFreal faceAvEy = 0;
//     CFreal faceAvEz = 0;
//     
//     // set the list of faces
//     
//     vector<SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
//     const CFuint nbTRSs = trs.size();    
//     
//     Common::SafePtr<TopologicalRegionSet> cells =
//     MeshDataStack::getActive()->getTrs("InnerCells");
//     CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
//     geoData.trs = cells;    
//     
//     const CFuint nbCells = cells->getLocalNbGeoEnts();
//     
//     for (CFuint iCell = 0; iCell < nbCells; ++iCell){
//       
//       divB[iCell] = 0;
//       divE[iCell] = 0;
//       
// //       //std::cout <<"In Element to compute the divergence 3D \n";
// 
//       
//       // set the builder data and build the GeometricEntity
//       geoData.idx = iCell;
//       GeometricEntity* currCell = m_geoBuilder.buildGE();
//       const CFuint elemID = currCell->getID(); 
//       const vector<GeometricEntity*>& faces = *currCell->getNeighborGeos();
//       const CFuint nbFaces = faces.size();
//     
//       std::vector<Node*> localNodes = *currCell->getNodes();
//       const CFuint nbNodesInElem = localNodes.size();
//       
// //       //std::cout <<"In Element to compute the divergence \n";
// //       //std::cout <<"Nb of Nodes in Element= \t "<< nbNodesInElem <<"\n";
// //       //std::cout <<"Nb of Faces in Element= \t "<< nbFaces <<"\n";
//       
//       for (CFuint i = 0; i < nbFaces; ++i){
// 	
// // 	//std::cout <<"In face to compute the divergence 3D \n";
// 	
// 	State *lState = faces[i]->getState(0);
// 	State *rState = faces[i]->getState(1);
// //  	//std::cout <<"At the face to compute the face average \n";
// 	
// 	const CFuint faceID = faces[i]->getID();
// 	const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
// 	CFreal nx = normals[startID];
// 	CFreal ny = normals[startID + 1];
// 	CFreal nz = normals[startID + 2];
// // 	CFLog( VERBOSE, "nx = " << nx << " ny = " << ny <<" norm = " << sqrt(nx*nx + ny*ny) << "\n");     
//       
// 	if (static_cast<CFuint>(isOutward[faceID]) != elemID) {
// 	  nx *= -1.;
// 	  ny *= -1.;
// 	  nz *= -1.;
// 	}
// //  	//std::cout <<"After taking the values of the cell and the neighbour\n";
// 	
// 	faceAvBx = 0.5*((*rState)[0] + (*lState)[0]);
// 	faceAvBy = 0.5*((*rState)[1] + (*lState)[1]);
// 	faceAvBz = 0.5*((*rState)[2] + (*lState)[2]);
// 	faceAvEx = 0.5*((*rState)[3] + (*lState)[3]);
// 	faceAvEy = 0.5*((*rState)[4] + (*lState)[4]);
// 	faceAvEz = 0.5*((*rState)[5] + (*lState)[5]);
// 	
// 	divB[iCell] += nx*faceAvBx + ny*faceAvBy + nz*faceAvBz;
// 	divE[iCell] += nx*faceAvEx + ny*faceAvEy + nz*faceAvEz;
//       }
//       divB[iCell] *= 1/volumes[elemID];
//       divE[iCell] *= 1/volumes[elemID]; 
// //        //std::cout <<"After calculating div*[iCell]\n";
//       m_geoBuilder.releaseGE();
//     }
//   }
}
//////////////////////////////////////////////////////////////////////////////
void DivMonitoring::computeAnaliticalSolution()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> BthetaTheory = socket_BthetaTheory.getDataHandle();  
  DataHandle<CFreal> ExPWTheory = socket_ExPWTheory.getDataHandle();  
  DataHandle<CFreal> EyPWTheory = socket_EyPWTheory.getDataHandle();  
  DataHandle<CFreal> EzPWTheory = socket_EzPWTheory.getDataHandle();  
  DataHandle<CFreal> BxPWTheory = socket_BxPWTheory.getDataHandle();  
  DataHandle<CFreal> ByPWTheory = socket_ByPWTheory.getDataHandle();  
  DataHandle<CFreal> BzPWTheory = socket_BzPWTheory.getDataHandle(); 
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();  
  //std::cout << "INcomputeAnaliticalSolution()\n";


  const CFreal pi = 3.14159265358979323846;  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal c_e = 299792458;
  
  if (dim == 2) {
    
  ///Compute Solution for Plane wave and Wire Testcase
    
    Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");
    CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
    geoData.trs = cells;    
    
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    const CFreal CurrentTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
    
    for (CFuint iCell = 0; iCell < nbCells; ++iCell){
      
      CFreal x = 0;
      CFreal y = 0;
      CFreal r = 0;
      CFreal r2 = 0;     
      
      CFreal Bx = 0;
      CFreal By = 0;    
      CFreal Bz = 0; 
      CFreal Ex = 0;   
      CFreal Ey = 0;
      CFreal Ez = 0;      
      
      // set the builder data and build the GeometricEntity
      geoData.idx = iCell;
      GeometricEntity* currCell = m_geoBuilder.buildGE();
      const CFuint elemID = currCell->getID(); 
      
      
      //Set the cartesian and radial coordinates 
      Node& coordinate = currCell->getState(0)->getCoordinates();
      x = coordinate[XX];
      y = coordinate[YY];
      r2 = x*x + y*y;
      r = std::sqrt(r2);         
      
      //Set the states
      State *currState = currCell->getState(0);
      Bx = (*currState)[0];
      By = (*currState)[1];
      Bz = (*currState)[2];
      Ex = (*currState)[3];    
      Ey = (*currState)[4];     
      Ez = (*currState)[5];       
          
      ///Adimensional Solution of the Wire Problem
      BthetaTheory[iCell] = 1/r ; 
      
      ///Tranverse Electric Plane Wave Dimensional(TEPW)
      ExPWTheory[iCell] = -sin(pi*y)*cos(pi*x - std::sqrt(2)*pi*CurrentTime*c_e) ;			//Only for cases    
      EyPWTheory[iCell] = cos(pi*y)*sin(pi*x - std::sqrt(2)*pi*CurrentTime*c_e) ; 			//Only for cases
      BzPWTheory[iCell] = (std::sqrt(2)/c_e)*cos(pi*y)*sin(pi*x - std::sqrt(2)*pi*CurrentTime*c_e) ;	//Only for cases
      ///Tranverse Electric Plane Wave Adimensional(TEPW)
//     	ExPWTheory[iCell] = -sin(pi*y)*cos(pi*x - std::sqrt(2)*pi*CurrentTime) ;			//Only for adimensional cases    
//     	EyPWTheory[iCell] = cos(pi*y)*sin(pi*x - std::sqrt(2)*pi*CurrentTime) ; 			//Only for adimensional cases
//     	BzPWTheory[iCell] = std::sqrt(2)*cos(pi*y)*sin(pi*x - std::sqrt(2)*pi*CurrentTime) ;	//Only for adimensional cases

      ///Tranverse Magnetic Plane Wave Dimensional(TEPW)
      BxPWTheory[iCell] = -(1/c_e)*cos(pi*y)*sin(pi*x - std::sqrt(2)*pi*CurrentTime*c_e) ; 			//Only for dimensional cases    
      ByPWTheory[iCell] = (1/c_e)*sin(pi*y)*cos(pi*x - std::sqrt(2)*pi*CurrentTime*c_e) ; 			//Only for dimensional cases
      EzPWTheory[iCell] = -std::sqrt(2)*sin(pi*y)*cos(pi*x - std::sqrt(2)*pi*CurrentTime*c_e) ;	//Only for dimensional cases        
      ///Transverse Magnetic Plane Wave Adimensional(TEPW)
//     	BxPWTheory[iCell] = -cos(pi*y)*sin(pi*x - std::sqrt(2)*pi*CurrentTime) ; 		//Only for adimensional cases    
//     	ByPWTheory[iCell] = sin(pi*y)*cos(pi*x - std::sqrt(2)*pi*CurrentTime) ; 			//Only for adimensional cases
//     	EzPWTheory[iCell] = -std::sqrt(2)*sin(pi*y)*cos(pi*x - std::sqrt(2)*pi*CurrentTime) ;	//Only for adimensional cases
      
      m_geoBuilder.releaseGE();
    }
  } 
  else{
  ///Compute Solution for Coaxial Waveguide Testcase
    
    Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");
    CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
    geoData.trs = cells;    
    
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    const CFreal CurrentTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
    
    for (CFuint iCell = 0; iCell < nbCells; ++iCell){
      
      CFreal x = 0;
      CFreal y = 0; 
      CFreal z = 0;
      
      CFreal Bx = 0;
      CFreal By = 0;    
      CFreal Ex = 0;   
      CFreal Ey = 0;
      
      // set the builder data and build the GeometricEntity
      geoData.idx = iCell;
      GeometricEntity* currCell = m_geoBuilder.buildGE();
      const CFuint elemID = currCell->getID(); 
      
      
      //Set the cartesian and polar coordinates 
      Node& coordinate = currCell->getState(0)->getCoordinates();
      x = coordinate[XX];
      y = coordinate[YY];
      z = coordinate[ZZ];

      
      //Set the states
      State *currState = currCell->getState(0);
      Bx = (*currState)[0];
      By = (*currState)[1];
      Ex = (*currState)[3];    
      Ey = (*currState)[4];     
           
      
      ///Coaxial Waveguide Dimensional (TEPW)
      ExPWTheory[iCell] = (x/(x*x + y*y))*sin(2*pi*z - 2*pi*CurrentTime*c_e) ;			//Only for cases    
      EyPWTheory[iCell] = (y/(x*x + y*y))*sin(2*pi*z - 2*pi*CurrentTime*c_e) ; 			//Only for cases
      BzPWTheory[iCell] = 0 ;	//Only for cases      ///Tranverse Electric Plane Wave Adimensional(TEPW)

      BxPWTheory[iCell] = -(1/c_e)*(y/(x*x + y*y))*sin(2*pi*z - 2*pi*CurrentTime*c_e) ; 			//Only for dimensional cases    
      ByPWTheory[iCell] = (1/c_e)*(x/(x*x + y*y))*sin(2*pi*z - 2*pi*CurrentTime*c_e) ; 			//Only for dimensional cases
      EzPWTheory[iCell] = 0 ;	//Only for dimensional cases        
      
      m_geoBuilder.releaseGE();
    }
  }
  //std::cout << "OUTcomputeAnaliticalSolution()\n";
}
/////////////////////////////////////////////////////////////////////////////
void DivMonitoring::computeCylindricalComponents()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> theta = socket_theta.getDataHandle();
  DataHandle<CFreal> Bradial = socket_Bradial.getDataHandle();  
  DataHandle<CFreal> Btheta = socket_Btheta.getDataHandle();
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();  
  
  //std::cout << "INcomputeCylindricalComponents()\n";

     
  Common::SafePtr<TopologicalRegionSet> cells =
  MeshDataStack::getActive()->getTrs("InnerCells");
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;  
  
  
  const CFreal pi = 3.14159265358979323846;  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    
    // set the builder data and build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();
    
    CFreal x = 0;
    CFreal y = 0;
    
    CFreal Bx = 0;
    CFreal By = 0;    
    CFreal Ex = 0;   
    CFreal Ey = 0;
    
    //Set the cartesian and polar coordinates
    Node& coordinate = currCell->getState(0)->getCoordinates();
    x = coordinate[XX];
    y = coordinate[YY];
    
    //Set the thetas
    if ((y > 0) && (x > 0.0)) {
	theta[iCell] = atan(y/x);
    }
    else if ((y > 0) && (x < 0.0)) {
	theta[iCell] = atan(y/x) + pi;
    }
    else if ((y < 0) && (x < 0.0)) {
	theta[iCell] = atan(y/x) + pi;
    }    
    else if ((y < 0) && (x > 0.0)) {
	theta[iCell] = atan(y/x) + 2*pi;
    }
    
    //Set the states
    State *currState = currCell->getState(0);
    Bx = (*currState)[0];
    By = (*currState)[1];
    Ex = (*currState)[3];    
    Ey = (*currState)[4]; 
    
    Bradial[iCell] = Bx*cos(theta[iCell]) + By*sin(theta[iCell])  ;
    Btheta[iCell] = -Bx*sin(theta[iCell]) + By*cos(theta[iCell])  ;
    
    m_geoBuilder.releaseGE();
  }
  //std::cout << "OUTcomputeCylindricalComponents()\n"; 
}
/////////////////////////////////////////////////////////////////////////////

boost::filesystem::path DivMonitoring::constructFilename()
{
  const bool isParallel = PE::GetPE().IsParallel ();
  const std::string nsp = this->getMethodData().getNamespace();
  
  if (isParallel) {
    std::ostringstream fname;
#ifdef CF_HAVE_BOOST_1_85
    fname << boost::filesystem::path(m_nameOutputFileDivMonitoring).stem().string()
          << "-" << PE::GetPE().GetRank("Default")
          << boost::filesystem::path(m_nameOutputFileDivMonitoring).extension().string();
#else
    fname << boost::filesystem::basename(boost::filesystem::path(m_nameOutputFileDivMonitoring))
          << "-" << PE::GetPE().GetRank("Default")
          << boost::filesystem::extension(boost::filesystem::path(m_nameOutputFileDivMonitoring));
#endif  
  }
  
  CFout << "Writing Electric and Magnetic Field Divergence to : " << m_nameOutputFileDivMonitoring << "\n";
  return boost::filesystem::path(m_nameOutputFileDivMonitoring);
}

//////////////////////////////////////////////////////////////////////////////

void DivMonitoring::prepareOutputFile(std::ofstream& outputFile)
{
  outputFile << "TITLE  = Error Norms\n";
  outputFile << "VARIABLES = Iter  TEL2Error TML2Error FullL2Error L2DivB L2DivE Time NbCells\n";
}

//////////////////////////////////////////////////////////////////////////////

void DivMonitoring::writeOutputFile()
{
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");
    
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;    
  
//  const CFuint nbCells = cells->getLocalNbGeoEnts();  
  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  const CFreal CurrentTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
  const CFuint nbCells = cells->getLocalNbGeoEnts();  
  //std::cout << "IN writeOutputFile() \n";
  
  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  if (iter == 1) {
    ofstream& outputFile = fhandle->open(constructFilename());
    prepareOutputFile(outputFile); 
    outputFile << iter
               << " "
               << m_TEPWErrorL2Norm
               << " "
               << m_TMPWErrorL2Norm
               << " "
               << m_fullPWErrorL2Norm
               << " "		  
               << m_divBL2Norm
               << " "
               << m_divEL2Norm
               << " "
               << CurrentTime
               << " "
	       << nbCells
               << "\n"; 
    outputFile.close();     
  }
  else {

    
    ofstream& outputFile = fhandle->open(constructFilename(), ios::app);
    outputFile << iter
               << " "
               << m_TEPWErrorL2Norm
               << " "
               << m_TMPWErrorL2Norm
               << " "
               << m_fullPWErrorL2Norm
               << " "		  
               << m_divBL2Norm
               << " "
               << m_divEL2Norm
               << " "
               << CurrentTime  
               << " "
	       << nbCells               
               << "\n"; 
    outputFile.close();     
      
  }
 
  CFout << "Writing of Electric and Magnetic Field divergence finished." << "\n";
//  outputFile.close();
}

//////////////////////////////////////////////////////////////////////////////

// void DivMonitoring::printDivMonitoringToFile()
// {
//   CFAUTOTRACE;
// 
//   using namespace std;
// 
//   boost::filesystem::path fpath (m_nameOutputFileDivMonitoring);
//   fpath = PathAppender::getInstance().appendAllInfo( fpath, m_appendIter, m_appendTime);
// 
//   SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
//   ofstream& outputFile = fhandle->open(fpath,ios::app);
// 
//   const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
//   const CFreal actualDivB = divB;
// 
//   if (iter == 1) {
//     outputFile << "TITLE  =  Divergence Monitoring\n";
//     outputFile << "VARIABLES = Iter I\n";
//   }
//   outputFile << iter << " " << actualdivB << "\n";
// 
//   fhandle->close();
// }

//////////////////////////////////////////////////////////////////////////////

void DivMonitoring::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeMaxwell

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
