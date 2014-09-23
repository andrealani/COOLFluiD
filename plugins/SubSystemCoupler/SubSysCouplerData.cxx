#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "SubSystemCoupler/SubSystemCoupler.hh"
#include "SubSystemCoupler/SubSysCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<SubSysCouplerData>, SubSysCouplerData, SubSystemCouplerModule> nullCouplerComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void SubSysCouplerData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<bool> >("NonMatchingGeometry","True if the interfaces geometry are not matching.");
   options.addConfigOption< std::vector<CFreal> >("NonMatchingGeometryRotation","Rotation Angle around 00 for geometrically non-matching interfaces.");
   options.addConfigOption< std::vector<CFreal> >("NonMatchingGeometryVector","Translation Vectors for geometrically non-matching interfaces.");
   options.addConfigOption< std::vector<CFreal> >("NonMatchingGeometryThreshold","Threshold to use to exclude part of the interface when the interfaces geometry are not fully matching.");

   options.addConfigOption< std::vector<std::string> >
      ("PreVariableTransformers","Variable Transformers for each interface. Transformation before sending data");
   options.addConfigOption< std::vector<std::string> >
      ("PostVariableTransformers","Variable Transformers for each interface. Transformation after receiving data");

   options.addConfigOption< std::vector<std::string> >("CoordType","Type of coordinates: nodes/states/gauss/ghost/nodalgauss");
   options.addConfigOption< bool >("FileTransfer","Transfer data using files (otherwise, though datahandles");
}

//////////////////////////////////////////////////////////////////////////////

SubSysCouplerData::SubSysCouplerData(Common::SafePtr<Framework::Method> owner)
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

  _preTransformersStr = std::vector<std::string>();
  setParameter("PreVariableTransformers",&_preTransformersStr);

  _postTransformersStr = std::vector<std::string>();
  setParameter("PostVariableTransformers",&_postTransformersStr);

  _coordTypeStr = std::vector<std::string>();
  setParameter("CoordType",&_coordTypeStr);

  _isTransferFiles = true;
  setParameter("FileTransfer",&_isTransferFiles);

}


//////////////////////////////////////////////////////////////////////////////

SubSysCouplerData::~SubSysCouplerData()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubSysCouplerData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CouplerData::configure(args);

  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  cf_assert(_isGeometryNonMatching.size() == _nonMatchingGeometryRotation.size());
  cf_assert(physModel->getDim() * _nonMatchingGeometryRotation.size() == _nonMatchingGeometryVector.size());

  ///Transform angle into radians
  for(CFuint j=0;j<_nonMatchingGeometryRotation.size();j++)
  {
  _nonMatchingGeometryRotation[j] *= (MathTools::MathConsts::CFrealPi()/180.);
  }

  ///Configure the Variable Transformers
  configureVariablesTransformers(args);

  ///Configure the coordinates type
  cf_assert(_coordTypeStr.size()!=0.);
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
  }
//   std::string provider = "Null";
//   if (_updateVarStr != "Null") {
//     provider = (physModel->getConvectiveName() != "Null") ?
//       (physModel->getConvectiveName() + _updateVarStr) :
//       (physModel->getDiffusiveName() + _updateVarStr) ;
//   }
//
//   CFLog(INFO, "SubSysCouplerData::_updateVarStr = " << provider
//         << "\n");
//
//   _updateVarSet.reset(Environment::Factory<VarSet>::getInstance().
//                       getProvider(provider)->create());
//
//   cf_assert(_updateVarSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void SubSysCouplerData::configureVariablesTransformers( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  SharedPtr<SubSysCouplerData> thisPtr(this);

  cf_assert(_preTransformersStr.size() > 0);
  _preTransformers.resize(_preTransformersStr.size());
  for(CFuint i=0;i<_preTransformersStr.size();i++)
  {
    Common::SafePtr<BaseMethodStrategyProvider<SubSysCouplerData,PreVariableTransformer> > prov =
      Environment::Factory<PreVariableTransformer>::getInstance().getProvider(_preTransformersStr[i]);

    cf_assert(prov.isNotNull());
    _preTransformers[i] = prov->create(_preTransformersStr[i],thisPtr);
    configureNested(_preTransformers[i].getPtr(), args);

    cf_assert(_preTransformers[i].isNotNull());
  }

  cf_assert(_postTransformersStr.size() > 0);
  _postTransformers.resize(_postTransformersStr.size());

  for(CFuint i=0;i<_postTransformersStr.size();i++)
  {
    Common::SafePtr<BaseMethodStrategyProvider<SubSysCouplerData,PostVariableTransformer> > prov =
      Environment::Factory<PostVariableTransformer>::getInstance().getProvider(_postTransformersStr[i]);

    cf_assert(prov.isNotNull());
    _postTransformers[i] = prov->create(_postTransformersStr[i],thisPtr);
    configureNested(_postTransformers[i].getPtr(), args);

    cf_assert(_postTransformers[i].isNotNull());
  }

}

//////////////////////////////////////////////////////////////////////////////

void SubSysCouplerData::setup()
{

  // set up the GeometricEntity builders
  _stdTrsGeoBuilder.setup();
  _faceTrsGeoBuilder.setup();

  //set up the pre variable transformers
  for(CFuint i=0;i<_preTransformers.size();i++)
  {
    _preTransformers[i]->setup();
  }

  //set up the post variable transformers
  for(CFuint i=0;i<_postTransformers.size();i++)
  {
    _postTransformers[i]->setup();
  }

}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> SubSysCouplerData::getThisCoupledCoordName(const std::string& interfaceName,
                                                      const std::string& trs)
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

std::string SubSysCouplerData::getThisCoupledCoordName(const std::string& interfaceName,
                                              const std::string& trs,
                                              const std::string& coordType)
{
  const std::string dataType = "COORD";
  return thisNameBuilder(interfaceName, trs, coordType,dataType);
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> SubSysCouplerData::getOtherCoupledCoordName(const std::string& interfaceName,
                                                       const std::string& trs,
                                                       const CFuint& iProcessor)
{
  const std::string dataType = "COORD";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize((_coupledSubSystemsCoordType[idx]).size());
  for(CFuint iType = 0; iType < (_coupledSubSystemsCoordType[idx]).size(); iType++)
  {
    _tempStringVector[iType] =
        otherNameBuilder(interfaceName, trs, (_coupledSubSystemsCoordType[idx])[iType],dataType);

    if(Common::PE::GetPE().IsParallel()) _tempStringVector[iType] += ".P" + StringOps::to_str(iProcessor);
  }

  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string SubSysCouplerData::getOtherCoupledCoordName(const std::string& interfaceName,
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

vector<std::string> SubSysCouplerData::getThisCoupledAcceptedName(const std::string& interfaceName,
                                                         const std::string& trs,
                                                         const CFuint& iProcessor)
{
  const std::string dataType = "ISACCEPTED";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize(_coordType[idx].size());
  for(CFuint iType = 0; iType < _coordType[idx].size(); iType++)
  {
    _tempStringVector[iType] = thisNameBuilder(interfaceName, trs, (_coordType[idx])[iType],dataType);

    if(Common::PE::GetPE().IsParallel())
    {
      const CFuint iProcLocal = Common::PE::GetPE().GetRank();
      //accepted by processor
      _tempStringVector[iType] += ".P" + StringOps::to_str(iProcessor);
      //data from processor
      _tempStringVector[iType] += "P" + StringOps::to_str(iProcLocal);
    }

  }

  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> SubSysCouplerData::getThisCoupledAcceptedName(const std::string& interfaceName,
                                                         const std::string& trs)
{
  const std::string dataType = "ISACCEPTED";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize(_coordType[idx].size());
  for(CFuint iType = 0; iType < _coordType[idx].size(); iType++)
  {
    _tempStringVector[iType] = thisNameBuilder(interfaceName, trs, (_coordType[idx])[iType],dataType);
  }

  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string SubSysCouplerData::getThisCoupledAcceptedName(const std::string& interfaceName,
                                                 const std::string& trs,
                                                 const std::string& coordType,
                                                 const CFuint& iProcessor)
{
  const std::string dataType = "ISACCEPTED";
  std::string name = thisNameBuilder(interfaceName, trs, coordType,dataType);
  if(Common::PE::GetPE().IsParallel())
  {
    const CFuint iProcLocal = Common::PE::GetPE().GetRank();
    //accepted by processor
    name += ".P" + StringOps::to_str(iProcessor);
    //data from processor
    name += "P" + StringOps::to_str(iProcLocal);
  }

  return name;
}

//////////////////////////////////////////////////////////////////////////////

std::string SubSysCouplerData::getThisCoupledAcceptedName(const std::string& interfaceName,
                                                       const std::string& trs,
                                                       const std::string& coordType)
{
  const std::string dataType = "ISACCEPTED";
  std::string name = thisNameBuilder(interfaceName, trs, coordType,dataType);

  return name;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> SubSysCouplerData::getOtherCoupledAcceptedName(const std::string& interfaceName,
                                                          const std::string& trs,
                                                          const CFuint& iProcessor)
{
  const std::string dataType = "ISACCEPTED";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize((_coupledSubSystemsCoordType[idx]).size());
  for(CFuint iType = 0; iType < (_coupledSubSystemsCoordType[idx]).size(); iType++)
  {
    _tempStringVector[iType] =
      otherNameBuilder(interfaceName, trs, (_coupledSubSystemsCoordType[idx])[iType],dataType);

    if(Common::PE::GetPE().IsParallel())
    {
      const CFuint iProcLocal = Common::PE::GetPE().GetRank();
      //accepted by processor
      _tempStringVector[iType] += ".P" + StringOps::to_str(iProcLocal);
      //data from processor
      _tempStringVector[iType] += "P" + StringOps::to_str(iProcessor);
    }
  }

  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string SubSysCouplerData::getOtherCoupledAcceptedName(const std::string& interfaceName,
                                                  const std::string& trs,
                                                  const std::string& coordType,
                                                  const CFuint& iProcessor)
{
  const std::string dataType = "ISACCEPTED";
  std::string name = otherNameBuilder(interfaceName, trs, coordType,dataType);
  if(Common::PE::GetPE().IsParallel())
  {
    const CFuint iProcLocal = Common::PE::GetPE().GetRank();
    //accepted by processor
    name += ".P" + StringOps::to_str(iProcLocal);
    //data from processor
    name += "P" + StringOps::to_str(iProcessor);
  }

  return name;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> SubSysCouplerData::getThisCoupledDataName(const std::string& interfaceName,
                                                           const std::string& trs)
{
  const std::string dataType = "DATA";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize(_coordType[idx].size());
  for(CFuint iType = 0; iType < _coordType[idx].size(); iType++)
  {
    _tempStringVector[iType] = thisNameBuilder(interfaceName, trs, (_coordType[idx])[iType],dataType);
  }

  return _tempStringVector;
}


//////////////////////////////////////////////////////////////////////////////

vector<std::string> SubSysCouplerData::getThisCoupledDataName(const std::string& interfaceName,
                                                     const std::string& trs,
                                                     const CFuint& iProcessor)
{
  const std::string dataType = "DATA";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize(_coordType[idx].size());
  for(CFuint iType = 0; iType < _coordType[idx].size(); iType++)
  {
    _tempStringVector[iType] = thisNameBuilder(interfaceName, trs, (_coordType[idx])[iType],dataType);

    if(Common::PE::GetPE().IsParallel())
    {
      const CFuint iProcLocal = Common::PE::GetPE().GetRank();
      //accepted by processor
      _tempStringVector[iType] += ".P" + StringOps::to_str(iProcessor);
      //data from processor
      _tempStringVector[iType] += "P" + StringOps::to_str(iProcLocal);
    }

  }

  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string SubSysCouplerData::getThisCoupledDataName(const std::string& interfaceName,
                                             const std::string& trs,
                                             const std::string& coordType,
                                             const CFuint& iProcessor)
{
  const std::string dataType = "DATA";
  std::string name = thisNameBuilder(interfaceName, trs, coordType,dataType);
  if(Common::PE::GetPE().IsParallel())
  {
    const CFuint iProcLocal = Common::PE::GetPE().GetRank();
    //accepted by processor
    name += ".P" + StringOps::to_str(iProcessor);
    //data from processor
    name += "P" + StringOps::to_str(iProcLocal);
  }

  return name;
}

//////////////////////////////////////////////////////////////////////////////

std::string SubSysCouplerData::getThisCoupledDataName(const std::string& interfaceName,
                                                   const std::string& trs,
                                                   const std::string& coordType)
{
  const std::string dataType = "DATA";
  std::string name = thisNameBuilder(interfaceName, trs, coordType,dataType);

  return name;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> SubSysCouplerData::getOtherCoupledDataName(const std::string& interfaceName,
                                                      const std::string& trs,
                                                      const CFuint& iProcessor)
{
  const std::string dataType = "DATA";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize((_coupledSubSystemsCoordType[idx]).size());
  for(CFuint iType = 0; iType < (_coupledSubSystemsCoordType[idx]).size(); iType++)
  {
    _tempStringVector[iType] = otherNameBuilder(interfaceName, trs, (_coupledSubSystemsCoordType[idx])[iType],dataType);
    if(Common::PE::GetPE().IsParallel())
    {
      const CFuint iProcLocal = Common::PE::GetPE().GetRank();
      //accepted by processor
      _tempStringVector[iType] += ".P" + StringOps::to_str(iProcLocal);
      //data from processor
      _tempStringVector[iType] += "P" + StringOps::to_str(iProcessor);
    }
  }

  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string SubSysCouplerData::getOtherCoupledDataName(const std::string& interfaceName,
                                              const std::string& trs,
                                              const std::string& coordType,
                                              const CFuint& iProcessor)
{
  const std::string dataType = "DATA";
  std::string name = otherNameBuilder(interfaceName, trs, coordType,dataType);
  if(Common::PE::GetPE().IsParallel())
  {
    const CFuint iProcLocal = Common::PE::GetPE().GetRank();
    //accepted by processor
    name += ".P" + StringOps::to_str(iProcLocal);
    //data from processor
    name += "P" + StringOps::to_str(iProcessor);
  }

  return name;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> SubSysCouplerData::getThisCoupledConnectedName(const std::string& interfaceName,
                                                          const std::string& trs)
{
  const std::string dataType = "CONNECTED";
  const CFuint idx = _coupledInterfacesMap.find(interfaceName);
  _tempStringVector.resize(_coordType[idx].size());
  for(CFuint iType = 0; iType < _coordType[idx].size(); iType++)
  {
    _tempStringVector[iType] = thisNameBuilder(interfaceName, trs, (_coordType[idx])[iType],dataType);
  }

  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string SubSysCouplerData::getThisCoupledConnectedName(const std::string& interfaceName,
                                                  const std::string& trs,
                                                  const std::string& coordType)
{
  const std::string dataType = "CONNECTED";
  return thisNameBuilder(interfaceName, trs, coordType,dataType);
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> SubSysCouplerData::getOtherCoupledConnectedName(const std::string& interfaceName,
                                                           const std::string& trs)
{

  throw Common::ShouldNotBeHereException (FromHere(),"getOtherCoupledConnectedName");

  return _tempStringVector;
}

//////////////////////////////////////////////////////////////////////////////

std::string SubSysCouplerData::getOtherCoupledConnectedName(const std::string& interfaceName,
                                                   const std::string& trs,
                                                   const std::string& coordType)
{

  throw Common::ShouldNotBeHereException (FromHere(),"getOtherCoupledConnectedName");

  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

std::string SubSysCouplerData::thisNameBuilder(const std::string& interface,
                                      const std::string& trs,
                                      const std::string& coordType,
                                      const std::string& type)
{
  const std::string nameSpace = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(nameSpace);
  Common::SafePtr<SubSystemStatus> subSys = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);
  const std::string subsystem = subSys->getSubSystemName();

  return "COUPLING_" + interface + "_" + trs + "_" + nameSpace + "_" + subsystem + "_" + coordType + "_" + type;
}

//////////////////////////////////////////////////////////////////////////////

std::string SubSysCouplerData::otherNameBuilder(const std::string& interface,
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
