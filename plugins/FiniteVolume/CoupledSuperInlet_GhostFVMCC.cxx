#include "FiniteVolume/FiniteVolume.hh"
#include "CoupledSuperInlet_GhostFVMCC.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CoupledSuperInlet_GhostFVMCC,
                                 CellCenterFVMData,
                                 FiniteVolumeModule>
CoupledSuperInlet_GhostFVMCCFVMCCProvider("CoupledSuperInlet_GhostFVMCC");

//////////////////////////////////////////////////////////////////////////////

void CoupledSuperInlet_GhostFVMCC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Interface","Name of the Interface.");
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

CoupledSuperInlet_GhostFVMCC::CoupledSuperInlet_GhostFVMCC(const std::string& name) :
  FVMCC_BC(name),
  _sockets()
{
  addConfigOptionsTo(this);

  _interfaceName = "";
  setParameter("Interface",&_interfaceName);

  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);
}

//////////////////////////////////////////////////////////////////////////////

CoupledSuperInlet_GhostFVMCC::~CoupledSuperInlet_GhostFVMCC()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CoupledSuperInlet_GhostFVMCC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = FVMCC_BC::needsSockets();

//  copy(_sockets.getAllSinkSockets().begin(), _sockets.getAllSinkSockets().end(), back_inserter(result));

  std::vector<Common::SafePtr<BaseDataSocketSink> > result2 = _sockets.getAllSinkSockets();

  for(CFuint i=0;i<result2.size();++i)
  {
    result.push_back(result2[i]);
  }

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CoupledSuperInlet_GhostFVMCC::setup()
{
  CFAUTOTRACE;

  FVMCC_BC::setup();
}

//////////////////////////////////////////////////////////////////////////////

void CoupledSuperInlet_GhostFVMCC::setGhostState(GeometricEntity *const face)
 {
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  const CFuint faceID = face->getID();
  const CFuint faceIdx = MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs")->getIdxInTrs(faceID);
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  const std::string trsName = trs->getName();
  const std::string currentSubSystem = SubSystemStatusStack::getActive()->getSubSystemName();
  const std::string nameSpace = getMethodData().getNamespace();

  const std::string baseSocketName =
    "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Ghost_";

  DataHandle<RealVector> interfaceData =
    _sockets.getSocketSink<RealVector>(baseSocketName + "DATA")->getDataHandle();
  DataHandle<CFreal> isAccepted =
    _sockets.getSocketSink<CFreal>(baseSocketName + "ISACCEPTED")->getDataHandle();
  DataHandle<RealVector> coordinates =
    _sockets.getSocketSink<RealVector>(baseSocketName + "COORD")->getDataHandle();

  ///Check that the datahandle has same size as the TRS
  cf_assert(isAccepted.size() == trs->getLocalNbGeoEnts());

  //Then get the values at the face nodes from the DataHandle
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State faceValue;
  RealVector adimFaceValue(nbEqs);
  RealVector faceCoord = coordinates[faceIdx];

  if((isAccepted[faceIdx]>=0.) && (interfaceData[faceIdx].size() > 0))
  {
    faceValue = interfaceData[faceIdx];
  }
  else
  {
    _vFunction.evaluate(faceCoord, faceValue);
  }

  //Adimensionalize
  getMethodData().getUpdateVar()->setAdimensionalValues(faceValue, adimFaceValue);

  // set the ghost state
  for(CFuint iEq=0; iEq < nbEqs; ++iEq)
  {
    (*ghostState)[iEq] = 2.0*adimFaceValue[iEq] - (*innerState)[iEq];
  }

}

//////////////////////////////////////////////////////////////////////////////

void CoupledSuperInlet_GhostFVMCC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  FVMCC_BC::configure(args);

  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);

  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

  const std::string nameSpace = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(nameSpace);
  Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

  const std::string currentSubSystem = subsystemStatus->getSubSystemName();
  const std::vector<std::string>& trsNames = getTrsNames();

  for(CFuint iTRS = 0; iTRS < trsNames.size(); iTRS++)
  {
    const std::string trsName = trsNames[iTRS];
    const std::string baseSocketName =
      "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Ghost_";

    _sockets.createSocketSink<CFreal>(baseSocketName + "ISACCEPTED");
    _sockets.createSocketSink<RealVector>(baseSocketName + "DATA");
    _sockets.createSocketSink<RealVector>(baseSocketName + "COORD");
  }

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
