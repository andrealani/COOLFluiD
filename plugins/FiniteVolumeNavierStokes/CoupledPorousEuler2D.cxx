#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/CoupledPorousEuler2D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DCons.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CoupledPorousEuler2D, CellCenterFVMData, FiniteVolumeNavierStokesModule> CoupledPorousEuler2DFVMCCProvider("CoupledPorousEuler2DFVMCC");

//////////////////////////////////////////////////////////////////////////////

void CoupledPorousEuler2D::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Interface","Name of the Interface.");
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

CoupledPorousEuler2D::CoupledPorousEuler2D(const std::string& name) :
  FVMCC_BC(name),
  _sockets(),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
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

CoupledPorousEuler2D::~CoupledPorousEuler2D()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CoupledPorousEuler2D::needsSockets()
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

void CoupledPorousEuler2D::setup()
{
  CFAUTOTRACE;

  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DCons>();

  cf_assert(_varSet.isNotNull());

  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

}

//////////////////////////////////////////////////////////////////////////////

void CoupledPorousEuler2D::setGhostState(GeometricEntity *const face)
 {
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  const CFuint faceID = face->getID();
  const CFuint faceIdx = MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs")->getIdxInTrs(faceID);
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  const std::string trsName = trs->getName();
  const std::string currentSubSystem = SubSystemStatusStack::getActive()->getSubSystemName();
  const std::string nameSpace = getMethodData().getNamespace();

  std::string socketName = "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_DATA";
  DataHandle<RealVector> interfaceData = _sockets.getSocketSink<RealVector>(socketName)->getDataHandle();

  socketName = "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_ISACCEPTED";
  DataHandle<CFreal> isAccepted = _sockets.getSocketSink<CFreal>(socketName)->getDataHandle();

  socketName = "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_COORD";
  DataHandle<RealVector> coordinates = _sockets.getSocketSink<RealVector>(socketName)->getDataHandle();

  ///Check that the datahandle has same size as the TRS
  ///@todo HERE we should check if the values passed are the Nodal values!!
  cf_assert(isAccepted.size() == trs->getLocalNbGeoEnts());

  //Then get the values at the face nodes from the DataHandle
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State faceValue;
  RealVector adimFaceValue(nbEqs);
  RealVector faceCoord = coordinates[faceIdx];

  if((isAccepted[faceIdx] >= 0.) && (interfaceData[faceIdx].size() > 0))
  {
    faceValue = interfaceData[faceIdx];
  }
  else
  {
    _vFunction.evaluate(faceCoord, faceValue);
  }

  _varSet->setAdimensionalValues(faceValue, adimFaceValue);

  //Then compute the value on the face from the nodal values (average of the nodal values)
  //This is for Euler2DCons
  const CFreal rho  = adimFaceValue[0];
  const CFreal rhoU = adimFaceValue[1];
  const CFreal rhoV = adimFaceValue[2];
  const CFreal rhoE = adimFaceValue[3];

  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
  // unused // const CFreal R = _varSet->getModel()->getR();
  const CFreal rhoV2 = (rhoU*rhoU + rhoV*rhoV)/rho;

  const CFreal pressure = (gamma - 1.)*(rhoE - (0.5*rhoV2));
  // unused // const CFreal temperature = pressure/(R*rho);
  // unused // const CFreal u = rhoU/rho;
  // unused // const CFreal v = rhoV/rho;

  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  const CFreal uInner = _dataInnerState[EulerTerm::VX];
  const CFreal vInner = _dataInnerState[EulerTerm::VY];
  const CFreal rhoInner = _dataInnerState[EulerTerm::RHO];
  // unused // const CFreal vnInner = uInner*nx + vInner*ny;
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];
  // unused // const CFreal aInnerState = _dataInnerState[EulerTerm::A];
  // unused // const CFreal machInner = vnInner / aInnerState;

  // unused // const CFreal pOuterState = _dataInnerState[EulerTerm::P];
  ///Method 1: Pressure Loss Model
  //Compute the normal velocity from the resistance
//   CFreal resistance = 0.07;
//   CFreal porosity = 1./(1.+resistance);
//or
  CFreal porosity = 0.07;
  CFreal resistance = (1./porosity) - 1.;

  //Compute the normal velocity
  // unused // CFreal newVn = vnInner + (pressure - pInnerState)/resistance;
  CFreal Dvn = (pressure - pInnerState)/resistance;

  ///We need to compute u and v from Vn and Vt!!!!
  // pressure is kept the same...
  //We assume tangential velocity is unchanged
  // or newVt = oldVt;

  //if face is horizontal  v_new = v_old + Dvn
  //if face is vertical  u_new = u_old + Dvn

  ///in general: u_new = u_old + Dvn*nx
  ///          : v_new = v_old + Dvn*ny
  CFreal uInner_new = uInner + Dvn*nx;
  CFreal vInner_new = vInner + Dvn*ny;

  // set all the physical data corresponding to the ghost state
  _dataGhostState[EulerTerm::RHO] = rhoInner;
  _dataGhostState[EulerTerm::VX] = 2.0*uInner_new - uInner;
  _dataGhostState[EulerTerm::VY] = 2.0*vInner_new - vInner;
  _dataGhostState[EulerTerm::V] = _dataInnerState[EulerTerm::V];
  _dataGhostState[EulerTerm::P] = _dataInnerState[EulerTerm::P];
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
       + 0.5*_dataGhostState[EulerTerm::RHO]*_dataGhostState[EulerTerm::V]*
       _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
				       _dataGhostState[EulerTerm::RHO]);
  
  _dataGhostState[EulerTerm::T] = _dataGhostState[EulerTerm::P]/(_dataGhostState[EulerTerm::RHO]*_varSet->getModel()->getR());
  
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

void CoupledPorousEuler2D::configure ( Config::ConfigArgs& args )
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
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(nameSpace);
  Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

  const std::string currentSubSystem = subsystemStatus->getSubSystemName();
  const std::vector<std::string>& trsNames = getTrsNames();

  for(CFuint iTRS = 0; iTRS < trsNames.size(); iTRS++)
  {
    const std::string trsName = trsNames[iTRS];
    std::string socketName = "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_ISACCEPTED";
    _sockets.createSocketSink<CFreal>(socketName);

    socketName = "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_DATA";
    _sockets.createSocketSink<RealVector>(socketName);

    socketName = "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Nodes_COORD";
    _sockets.createSocketSink<RealVector>(socketName);
  }

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
