#include "FiniteVolumeMHD/iPic3DCoupler.hh"
#include "Common/PE.hh"
#include "Common/OSystem.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "MathTools/MathConsts.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Common/BadValueException.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "FiniteVolume/CellCenterFVM.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeMHD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<iPic3DCoupler, DataProcessingData, FiniteVolumeMHDModule>
divMonitoringProvider("iPic3DCoupler");

//////////////////////////////////////////////////////////////////////////////

void iPic3DCoupler::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("SaveRate","Save Output File every...iterations.");
   options.addConfigOption< bool >("AppendTime","Save each iteration to different file with suffix _time#.");
   options.addConfigOption< bool >("AppendIter","Save each iteration to different file with suffix _iter#.");
   options.addConfigOption< std::string >("PathToInputFile","Name of input file for iPic3D");
}

//////////////////////////////////////////////////////////////////////////////
      
iPic3DCoupler::iPic3DCoupler(const std::string& name) :
  DataProcessingCom(name),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_uX("uX"),
  socket_uY("uY"), 
  socket_uZ("uZ"),
  socket_curlB("curlB"),
  m_geoBuilder(),
  m_updateVarSet(CFNULL)
{
  addConfigOptionsTo(this);

  m_saveRate = 1;
  setParameter("SaveRate",&m_saveRate);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);
  
  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);
  
  m_pathToInputFile = "";
  setParameter("PathToInputFile",&m_pathToInputFile);
}
      
//////////////////////////////////////////////////////////////////////////////

iPic3DCoupler::~iPic3DCoupler()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
iPic3DCoupler::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_uX);
  result.push_back(&socket_uY);
  result.push_back(&socket_uZ);
  
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> > iPic3DCoupler::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_curlB);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void iPic3DCoupler::setup()
{
  CFAUTOTRACE;
  
  Common::SafePtr<TopologicalRegionSet> cells = 
    MeshDataStack::getActive()->getTrs("InnerCells");  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  // geometry builder setup
  m_geoBuilder.setup();
  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  m_geoBuilder.getDataGE().trs = cells;
  
  m_updateVarSet = getMethodData().getUpdateVarSet().d_castTo<MHD3DProjectionVarSet>();
  
  socket_curlB.getDataHandle().resize(nbCells*3);
}
      
//////////////////////////////////////////////////////////////////////////////

void iPic3DCoupler::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void iPic3DCoupler::execute()
{
  CFAUTOTRACE;
  
  // CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  // DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < CFreal > curlB = socket_curlB.getDataHandle();
  
   Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");  
  // CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  // geoData.trs = cells;  
  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  //  const CFuint nbProc = PE::GetPE().GetProcessorCount();
  //   if (PE::GetPE().GetRank() == 0) {
  //     // std::string command1 = "mpirun -np " + StringOps::to_str(nbProc) + " ./iPIC3D " + m_pathToInputFile;
  //     // cout << "Executing: " << command1 << endl;
  //     std::string command1 = "mpirun -np 2 ./iPIC3D " + m_pathToInputFile;
  //     Common::OSystem::getInstance().executeCommand(command1);
  //    } 
  
  //   PE::GetPE().setBarrier(); 
  
  // fill in the curlB array
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> uZ = socket_uZ.getDataHandle();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    if (states[iCell]->isParUpdatable()) {
      const CFuint start  = iCell*3;
      const CFuint startB = iCell*nbEqs;
      const CFuint BxID  = startB + 4;
      const CFuint ByID  = startB + 5;
      const CFuint BzID  = startB + 6;
      curlB[start + XX] =    uY[BzID] - uZ[ByID];
      curlB[start + YY] = - (uX[BzID] - uZ[BxID]);
      curlB[start + ZZ] =    uX[ByID] - uY[BxID];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void iPic3DCoupler::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeMHD

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
