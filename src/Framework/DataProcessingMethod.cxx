// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/DataProcessingMethod.hh"
#include "Framework/SubSystemStatus.hh"
#include "Environment/CFEnv.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void DataProcessingMethod::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint, Config::DynamicOption<> >("ProcessRate","Rate to process the data."); 
  options.addConfigOption< CFuint, Config::DynamicOption<> >("StartIter","Iteration at which processing starts.");
  options.addConfigOption< CFuint, Config::DynamicOption<> >("StopIter", "Iteration to stop processing stops.");
  options.addConfigOption< bool >("NeedsInitialization","Flag telling if this processing needs initialization.");
}

//////////////////////////////////////////////////////////////////////////////

void DataProcessingMethod::build_dynamic_functions()
{
  add_dynamic_function("processData",&DataProcessingMethod::processData);
}

//////////////////////////////////////////////////////////////////////////////

void DataProcessingMethod::registActionListeners()
{
  Method::registActionListeners();
  
  //  const std::string ssname = SubSystemStatusStack::getCurrentName();  
  //  event_handler->addListener(event_handler->key(ssname, "CF_ON_DATAPROCESS_PROCESS"),
  //                              this,&DataProcessingMethod::solve);
}

//////////////////////////////////////////////////////////////////////////////

DataProcessingMethod::DataProcessingMethod(const std::string& name)
  : Method(name),
    m_convergenceMtd(),
    m_lssMtd()
{
  // define which functions might be called dynamic
  build_dynamic_functions();
  // regist which functions might be called by raising Events
  registActionListeners();
  // regist the configuration options
  addConfigOptionsTo(this);

  // by default the data processing is run once -> processRate=infinity
  m_processRate = numeric_limits<CFuint>::max();
  setParameter("ProcessRate",&m_processRate);

  m_startIter = 0;
  setParameter("StartIter",&m_startIter);

  m_stopIter = numeric_limits<CFuint>::max();
  setParameter("StopIter",&m_stopIter);  

  m_needsInitialization = false;
  setParameter("NeedsInitialization",&m_needsInitialization);
}
    
//////////////////////////////////////////////////////////////////////////////

DataProcessingMethod::~DataProcessingMethod()
{
}

//////////////////////////////////////////////////////////////////////////////

void DataProcessingMethod::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  Method::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DataProcessingMethod::processData()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();
  
  if (SubSystemStatusStack::getActive()->getNbIter() < m_stopIter 
      && SubSystemStatusStack::getActive()->getNbIter() >= m_startIter ) {
    CFLog(VERBOSE, "DataProcessingMethod::processData() for [" << getName() << "]\n");
    processDataImpl();
  }
  
  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void DataProcessingMethod::setMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void DataProcessingMethod::unsetMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
