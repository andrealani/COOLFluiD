#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "ConcurrentCoupler/ConcurrentCoupler.hh"
#include "ConcurrentCoupler/ConcurrentCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ConcurrentCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<ConcurrentCouplerData>, 
		      ConcurrentCouplerData, ConcurrentCouplerModule> 
nullCouplerComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<bool> >("NonMatchingGeometry","True if the interfaces geometry are not matching.");
  options.addConfigOption< std::vector<CFreal> >("NonMatchingGeometryRotation","Rotation Angle around 00 for geometrically non-matching interfaces.");
  options.addConfigOption< std::vector<CFreal> >("NonMatchingGeometryVector","Translation Vectors for geometrically non-matching interfaces.");
  options.addConfigOption< std::vector<CFreal> >("NonMatchingGeometryThreshold","Threshold to use to exclude part of the interface when the interfaces geometry are not fully matching.");
  
  options.addConfigOption< std::vector<std::string> >("CoordType","Type of coordinates: nodes/states/gauss/ghost/nodalgauss");
  options.addConfigOption< bool >("FileTransfer","Transfer data using files (otherwise, though datahandles");
}
      
//////////////////////////////////////////////////////////////////////////////

ConcurrentCouplerData::ConcurrentCouplerData(Common::SafePtr<Framework::Method> owner)
  : CouplerData(owner),
    _stdTrsGeoBuilder(),
    _faceTrsGeoBuilder(),
    _coupledInterfaces(),
    _interfaces()
{
  addConfigOptionsTo(this);

  _isGeometryNonMatching = std::vector<bool>();
  setParameter("NonMatchingGeometry",&_isGeometryNonMatching);

  _defaultNonMatchingGeometryThreshold  = 0.0000000001;
  _nonMatchingGeometryThreshold = std::vector<CFreal>();
  setParameter("NonMatchingGeometryThreshold",&_nonMatchingGeometryThreshold);

  _nonMatchingGeometryRotation = std::vector<CFreal>();
  setParameter("NonMatchingGeometryRotation",&_nonMatchingGeometryRotation);

  _nonMatchingGeometryVector = std::vector<CFreal>();
  setParameter("NonMatchingGeometryVector",&_nonMatchingGeometryVector);

  _coordTypeStr = std::vector<std::string>();
  setParameter("CoordType",&_coordTypeStr);

  _isTransferFiles = true;
  setParameter("FileTransfer",&_isTransferFiles);
}


//////////////////////////////////////////////////////////////////////////////

ConcurrentCouplerData::~ConcurrentCouplerData()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CouplerData::configure(args);

  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  cf_assert(_isGeometryNonMatching.size() == _nonMatchingGeometryRotation.size());
  cf_assert(physModel->getDim() * _nonMatchingGeometryRotation.size() == _nonMatchingGeometryVector.size());

  ///Transform angle into radians
  for(CFuint j=0;j<_nonMatchingGeometryRotation.size();j++)
  {
  _nonMatchingGeometryRotation[j] *= (MathTools::MathConsts::CFrealPi()/180.);
  }
  
  ///Configure the coordinates type
  /* cf_assert(_coordTypeStr.size()!=0.);
  _coordType.resize(_coordTypeStr.size());
  for(CFuint i = 0; i< _coordTypeStr.size(); i++)
  {

    if((_coordTypeStr[i] != "NodesGauss")
        && (_coordTypeStr[i] != "StatesGauss")
        && (_coordTypeStr[i] != "NodesStates")
        && (_coordTypeStr[i] != "StatesNodes"))
    {
      cf_assert((_coordTypeStr[i] == "Nodes")||(_coordTypeStr[i] == "States")||(_coordTypeStr[i] == "Gauss")||(_coordTypeStr[i] == "Ghost"));
      _coordType[i].resize(1);
      (_coordType[i])[0] = _coordTypeStr[i];
    }
    else
    {
      if(_coordTypeStr[i] == "NodesGauss")
      {
        _coordType[i].resize(2);
        (_coordType[i])[0] = "Nodes";
        (_coordType[i])[1] = "Gauss";
      }
      if(_coordTypeStr[i] == "StatesGauss")
      {
        _coordType[i].resize(2);
        (_coordType[i])[0] = "States";
        (_coordType[i])[1] = "Gauss";
      }
      if((_coordTypeStr[i] == "NodesStates") || (_coordTypeStr[i] == "StatesNodes"))
      {
        _coordType[i].resize(2);
        (_coordType[i])[0] = "Nodes";
        (_coordType[i])[1] = "States";
      }
    }
    }*/
}

//////////////////////////////////////////////////////////////////////////////

void ConcurrentCouplerData::setup()
{
  // set up the GeometricEntity builders
  _stdTrsGeoBuilder.setup();
  _faceTrsGeoBuilder.setup();
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> ConcurrentCouplerData::getThisCoupledCoordName
(const std::string& interfaceName, const std::string& trs)
{
  const std::string dataType = "COORD";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize(_coordType[idx].size());
  for(CFuint iType = 0; iType < _coordType[idx].size(); iType++)
  {
    _tempStringVector[iType] = thisNameBuilder(interfaceName, trs, (_coordType[idx])[iType],dataType);
  }

  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConcurrentCouplerData::getThisCoupledCoordName
(const std::string& interfaceName, const std::string& trs, const std::string& coordType)
{
  const std::string dataType = "COORD";
  return thisNameBuilder(interfaceName, trs, coordType,dataType);
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> ConcurrentCouplerData::getOtherCoupledCoordName
(const std::string& interfaceName, const std::string& trs, const CFuint& iProcessor)
{
  const std::string dataType = "COORD";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize((_coupledSubSystemsCoordType[idx]).size());
  for(CFuint iType = 0; iType < (_coupledSubSystemsCoordType[idx]).size(); iType++) {
    _tempStringVector[iType] =
      otherNameBuilder(interfaceName, trs, (_coupledSubSystemsCoordType[idx])[iType],dataType);
    
    if(Common::PE::GetPE().IsParallel()) _tempStringVector[iType] += ".P" + StringOps::to_str(iProcessor);
  }
  
  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConcurrentCouplerData::getOtherCoupledCoordName
(const std::string& interfaceName,
 const std::string& trs,
 const std::string& coordType,
 const CFuint& iProcessor)
{
  const std::string dataType = "COORD";
  std::string name = otherNameBuilder(interfaceName, trs, coordType,dataType);
  if(Common::PE::GetPE().IsParallel()) name += ".P" + StringOps::to_str(iProcessor);
  return name;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> ConcurrentCouplerData::getThisCoupledAcceptedName
(const std::string& interfaceName,
 const std::string& trs,
 const CFuint& iProcessor)
{
  const std::string dataType = "ISACCEPTED";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize(_coordType[idx].size());
  for(CFuint iType = 0; iType < _coordType[idx].size(); iType++) {
    _tempStringVector[iType] = thisNameBuilder(interfaceName, trs, (_coordType[idx])[iType],dataType);
    
    if(Common::PE::GetPE().IsParallel()) {
      const CFuint iProcLocal = Common::PE::GetPE().GetRank(getNamespace());
      //accepted by processor
      _tempStringVector[iType] += ".P" + StringOps::to_str(iProcessor);
      //data from processor
      _tempStringVector[iType] += "P" + StringOps::to_str(iProcLocal);
    }
  }
  
  return _tempStringVector;
}
      
//////////////////////////////////////////////////////////////////////////////

vector<std::string> ConcurrentCouplerData::getThisCoupledAcceptedName
(const std::string& interfaceName, const std::string& trs)
{
  const std::string dataType = "ISACCEPTED";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize(_coordType[idx].size());
  for(CFuint iType = 0; iType < _coordType[idx].size(); iType++) {
    _tempStringVector[iType] = thisNameBuilder(interfaceName, trs, (_coordType[idx])[iType],dataType);
  }
  
  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConcurrentCouplerData::getThisCoupledAcceptedName
(const std::string& interfaceName,
 const std::string& trs,
 const std::string& coordType,
 const CFuint& iProcessor)
{
  const std::string dataType = "ISACCEPTED";
  std::string name = thisNameBuilder(interfaceName, trs, coordType,dataType);
  if(Common::PE::GetPE().IsParallel()) {  
    const std::string nsp = getNamespace();
    const CFuint iProcLocal = Common::PE::GetPE().GetRank(nsp);
    //accepted by processor
    name += ".P" + StringOps::to_str(iProcessor);
    //data from processor
    name += "P" + StringOps::to_str(iProcLocal);
  }
  
  return name;
}

//////////////////////////////////////////////////////////////////////////////
      
std::string ConcurrentCouplerData::getThisCoupledAcceptedName
(const std::string& interfaceName,
 const std::string& trs,
 const std::string& coordType)
{
  const std::string dataType = "ISACCEPTED";
  std::string name = thisNameBuilder(interfaceName, trs, coordType,dataType);
  return name;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> ConcurrentCouplerData::getOtherCoupledAcceptedName
(const std::string& interfaceName,
 const std::string& trs,
 const CFuint& iProcessor)
{
  const std::string dataType = "ISACCEPTED";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize((_coupledSubSystemsCoordType[idx]).size());
  for(CFuint iType = 0; iType < (_coupledSubSystemsCoordType[idx]).size(); iType++) {
    _tempStringVector[iType] =
      otherNameBuilder(interfaceName, trs, (_coupledSubSystemsCoordType[idx])[iType],dataType);
    
    if(Common::PE::GetPE().IsParallel()) {
      const CFuint iProcLocal = Common::PE::GetPE().GetRank(getNamespace());
      //accepted by processor
      _tempStringVector[iType] += ".P" + StringOps::to_str(iProcLocal);
      //data from processor
      _tempStringVector[iType] += "P" + StringOps::to_str(iProcessor);
    }
  }
  
  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConcurrentCouplerData::getOtherCoupledAcceptedName
(const std::string& interfaceName,
 const std::string& trs,
 const std::string& coordType,
 const CFuint& iProcessor)
{
  const std::string dataType = "ISACCEPTED";
  std::string name = otherNameBuilder(interfaceName, trs, coordType,dataType);
  if(Common::PE::GetPE().IsParallel()) {
    const CFuint iProcLocal = Common::PE::GetPE().GetRank(getNamespace());
    //accepted by processor
    name += ".P" + StringOps::to_str(iProcLocal);
    //data from processor
    name += "P" + StringOps::to_str(iProcessor);
  }
  
  return name;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> ConcurrentCouplerData::getThisCoupledDataName
(const std::string& interfaceName, const std::string& trs)
{
  const std::string dataType = "DATA";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize(_coordType[idx].size());
  for(CFuint iType = 0; iType < _coordType[idx].size(); iType++) {
    _tempStringVector[iType] = thisNameBuilder(interfaceName, trs, (_coordType[idx])[iType],dataType);
  }
  
  return _tempStringVector;
}


//////////////////////////////////////////////////////////////////////////////

vector<std::string> ConcurrentCouplerData::getThisCoupledDataName
(const std::string& interfaceName,
 const std::string& trs,
 const CFuint& iProcessor)
{
  const std::string dataType = "DATA";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize(_coordType[idx].size());
  for(CFuint iType = 0; iType < _coordType[idx].size(); iType++) {
    _tempStringVector[iType] = thisNameBuilder(interfaceName, trs, (_coordType[idx])[iType],dataType);

    if(Common::PE::GetPE().IsParallel()) {
      const CFuint iProcLocal = Common::PE::GetPE().GetRank(getNamespace());
      //accepted by processor
      _tempStringVector[iType] += ".P" + StringOps::to_str(iProcessor);
      //data from processor
      _tempStringVector[iType] += "P" + StringOps::to_str(iProcLocal);
    }
  }
  
  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////
      
std::string ConcurrentCouplerData::getThisCoupledDataName(const std::string& interfaceName,
							  const std::string& trs,
							  const std::string& coordType,
							  const CFuint& iProcessor)
{
  const std::string dataType = "DATA";
  std::string name = thisNameBuilder(interfaceName, trs, coordType,dataType);
  if(Common::PE::GetPE().IsParallel()) {
    const CFuint iProcLocal = Common::PE::GetPE().GetRank(getNamespace());
    //accepted by processor
    name += ".P" + StringOps::to_str(iProcessor);
    //data from processor
    name += "P" + StringOps::to_str(iProcLocal);
  }
  
  return name;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConcurrentCouplerData::getThisCoupledDataName
(const std::string& interfaceName,
 const std::string& trs,
 const std::string& coordType)
{
  const std::string dataType = "DATA";
  std::string name = thisNameBuilder(interfaceName, trs, coordType,dataType);
  return name;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> ConcurrentCouplerData::getOtherCoupledDataName
(const std::string& interfaceName,
 const std::string& trs,
 const CFuint& iProcessor)
{
  const std::string dataType = "DATA";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize((_coupledSubSystemsCoordType[idx]).size());
  for(CFuint iType = 0; iType < (_coupledSubSystemsCoordType[idx]).size(); iType++) {
    _tempStringVector[iType] = otherNameBuilder(interfaceName, trs, (_coupledSubSystemsCoordType[idx])[iType],dataType);
    if(Common::PE::GetPE().IsParallel()) {
      const CFuint iProcLocal = Common::PE::GetPE().GetRank(getNamespace());
      //accepted by processor
      _tempStringVector[iType] += ".P" + StringOps::to_str(iProcLocal);
      //data from processor
      _tempStringVector[iType] += "P" + StringOps::to_str(iProcessor);
    }
  }
  
  return _tempStringVector;
}
      
//////////////////////////////////////////////////////////////////////////////

std::string ConcurrentCouplerData::getOtherCoupledDataName
(const std::string& interfaceName,
 const std::string& trs,
 const std::string& coordType,
 const CFuint& iProcessor)
{
  const std::string dataType = "DATA";
  std::string name = otherNameBuilder(interfaceName, trs, coordType,dataType);
  if(Common::PE::GetPE().IsParallel()) {
    const CFuint iProcLocal = Common::PE::GetPE().GetRank(getNamespace());
    //accepted by processor
    name += ".P" + StringOps::to_str(iProcLocal);
    //data from processor
    name += "P" + StringOps::to_str(iProcessor);
  }
  
  return name;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> ConcurrentCouplerData::getThisCoupledConnectedName
(const std::string& interfaceName, const std::string& trs)
{
  const std::string dataType = "CONNECTED";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize(_coordType[idx].size());
  for(CFuint iType = 0; iType < _coordType[idx].size(); iType++) {
    _tempStringVector[iType] = thisNameBuilder(interfaceName, trs, (_coordType[idx])[iType],dataType);
  }
  
  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConcurrentCouplerData::getThisCoupledConnectedName
(const std::string& interfaceName,
 const std::string& trs,
 const std::string& coordType)
{
  const std::string dataType = "CONNECTED";
  return thisNameBuilder(interfaceName, trs, coordType,dataType);
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> ConcurrentCouplerData::getOtherCoupledConnectedName
(const std::string& interfaceName,
 const std::string& trs)
{
  throw Common::ShouldNotBeHereException (FromHere(),"getOtherCoupledConnectedName");
  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConcurrentCouplerData::getOtherCoupledConnectedName
(const std::string& interfaceName,
 const std::string& trs,
 const std::string& coordType)
{
  throw Common::ShouldNotBeHereException (FromHere(),"getOtherCoupledConnectedName");
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConcurrentCouplerData::thisNameBuilder(const std::string& interface,
						   const std::string& trs,
						   const std::string& coordType,
						   const std::string& type)
{
  const std::string nameSpace = getNamespace();
  
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(nameSpace);
  Common::SafePtr<SubSystemStatus> subSys = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);
  const std::string subsystem = subSys->getSubSystemName();
  
  return "COUPLING_" + interface + "_" + trs + "_" + nameSpace + "_" + subsystem + "_" + coordType + "_" + type;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConcurrentCouplerData::otherNameBuilder(const std::string& interface,
						    const std::string& trs,
						    const std::string& coordType,
						    const std::string& type)
{
  const std::string nameSpace = getCoupledNameSpaceName(interface);
  const std::string subsystem = getCoupledSubSystemName(interface);
  
  return "COUPLING_" + interface + "_" + trs + "_" + nameSpace + "_" + subsystem + "_" + coordType + "_" + type;
}
      
//////////////////////////////////////////////////////////////////////////////


    } // namespace SubSystem

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
