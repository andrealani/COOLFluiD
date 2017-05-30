// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/DataProcessingData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/Framework.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<DataProcessingData>,
		      DataProcessingData, FrameworkLib>
nullDataProcessingComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void DataProcessingData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("updateVar","update variables set");
  options.addConfigOption< CFuint, Config::DynamicOption<> >("StartIter","Iteration at which postprocessing must start.");
}

//////////////////////////////////////////////////////////////////////////////

DataProcessingData::DataProcessingData(Common::SafePtr<Method> owner)
 : MethodData(owner),
   _spaceMtd(),
   _stdTrsGeoBuilder(),
   _faceTrsGeoBuilder(),
   _updateVarStr(),
   _updateVarSet()
{
  addConfigOptionsTo(this);
  _updateVarStr = "Null";
  setParameter("updateVar",&_updateVarStr);
  
  ///@info in order for this to be issued it must be called from the derived class through the command 
  //getMethodData().getStartIter()
  m_startIter = 0.;
  setParameter("StartIter",&m_startIter);
}
    
//////////////////////////////////////////////////////////////////////////////

DataProcessingData::~DataProcessingData()
{
}

//////////////////////////////////////////////////////////////////////////////

void DataProcessingData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  
  MethodData::configure(args);
  
  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = 
    PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  std::string provider = "Null";
  if (_updateVarStr != "Null") {
    provider = (physModel->getConvectiveName() != "Null") ?
      (physModel->getConvectiveName() + _updateVarStr) :
      (physModel->getDiffusiveName() + _updateVarStr) ;
  }
  
  CFLog(VERBOSE, "DataProcessingData : UpdateVarSet is " << provider << "\n");
  
  _updateVarSet.reset(FACTORY_GET_PROVIDER
    (getFactoryRegistry(), ConvectiveVarSet, provider)->
    create(physModel->getImplementor()->getConvectiveTerm()));
  
  cf_assert(_updateVarSet.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void DataProcessingData::setup()
{
  CFAUTOTRACE;

  CFLog ( VERBOSE, " +++ DataProcessingData::setup()\n" );  
  
  MethodData::setup();
 
  _stdTrsGeoBuilder.setup();
  _faceTrsGeoBuilder.setup();

}

//////////////////////////////////////////////////////////////////////////////

void DataProcessingData::unsetup()
{
  CFAUTOTRACE;
  
  CFLog ( VERBOSE, " +++ DataProcessingData::unsetup()\n" );
  
  _faceTrsGeoBuilder.unsetup();
  _stdTrsGeoBuilder.unsetup();
  
  MethodData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

