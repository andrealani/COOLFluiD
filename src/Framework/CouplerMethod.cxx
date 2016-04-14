// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CouplerMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::build_dynamic_functions()
{
  add_dynamic_function("preProcessWrite",&CouplerMethod::preProcessWrite);
  add_dynamic_function("preProcessRead",&CouplerMethod::preProcessRead);
  add_dynamic_function("meshMatchingWrite",&CouplerMethod::meshMatchingWrite);
  add_dynamic_function("meshMatchingRead",&CouplerMethod::meshMatchingRead);
  add_dynamic_function("dataTransferRead",&CouplerMethod::dataTransferRead);
  add_dynamic_function("dataTransferWrite",&CouplerMethod::dataTransferWrite);
  add_dynamic_function("finalize",&CouplerMethod::finalize);
}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::registActionListeners()
{
  Method::registActionListeners();

  // const std::string ssname = SubSystemStatusStack::getCurrentName();   
  // event_handler->addListener(event_handler->key(ssname, "CF_ON_COUPLERMETHOD_MATCH"),
  //                            this,&CouplerMethod::match);
}

//////////////////////////////////////////////////////////////////////////////

CouplerMethod::CouplerMethod(const std::string& name)  : Method(name)
{
  // define which functions might be called dynamic
  build_dynamic_functions();
  // regist which functions might be called by raising Events
  registActionListeners();
  // regist the configuration options
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

CouplerMethod::~CouplerMethod()
{
}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  Method::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::setInfoToOtherSubSystem()
{
  CFAUTOTRACE;

}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::getInfoFromOtherSubSystem()
{
  CFAUTOTRACE;

}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::postConfigure( Config::ConfigArgs& args)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::preProcessWrite()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE,"-------------------------------------------------------------\n");
  CFLog(VERBOSE,"Writing coordinates of CouplerMethod [" << getName() << "]\n");
  cf_assert(isConfigured());
  cf_assert(isSetup());
  
  pushNamespace();
  
  preProcessWriteImpl();
  
  popNamespace();
  CFLog(VERBOSE,"-------------------------------------------------------------\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::preProcessRead()
{
  CFAUTOTRACE;

  CFLog(VERBOSE,"-------------------------------------------------------------\n");
  CFLog(VERBOSE,"Reading coordinates of CouplerMethod [" << getName() << "]\n");
  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  preProcessReadImpl();

  popNamespace();
  CFLog(VERBOSE,"-------------------------------------------------------------\n");
}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::meshMatchingWrite()
{
  CFAUTOTRACE;

  CFLog(VERBOSE,"-------------------------------------------------------------\n");
  CFLog(VERBOSE,"Mesh matching of CouplerMethod [" << getName() << "]\n");
  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  meshMatchingWriteImpl();

  popNamespace();
  CFLog(VERBOSE,"-------------------------------------------------------------\n");
}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::meshMatchingRead()
{
  CFAUTOTRACE;

  CFLog(VERBOSE,"-------------------------------------------------------------\n");
  CFLog(VERBOSE,"Reading matched mesh of CouplerMethod [" << getName() << "]\n");
  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  meshMatchingReadImpl();

  popNamespace();
  CFLog(VERBOSE,"-------------------------------------------------------------\n");
}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::dataTransferRead()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  dataTransferReadImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::dataTransferWrite()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  dataTransferWriteImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::finalize()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();
  
  finalizeImpl();
  
  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::setMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void CouplerMethod::unsetMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
