// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include <sstream>

#include "Common/PE.hh"

#include "Common/ProcessInfo.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Environment/CFEnv.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PathAppender.hh"
#include "Framework/DynamicBalancerMethod.hh"
//#include "Framework/DynamicBalancerMethodData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace boost::filesystem;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void DynamicBalancerMethod::defineConfigOptions(Config::OptionList& options)
{
  // some options
}

//////////////////////////////////////////////////////////////////////////////

void DynamicBalancerMethod::build_dynamic_functions()
{
}

//////////////////////////////////////////////////////////////////////////////

void DynamicBalancerMethod::registActionListeners()
{
  Method::registActionListeners();
  // const std::string ssname = SubSystemStatusStack::getCurrentName();   
  // event_handler->addListener(event_handler->key(ssname, "CF_ON_CONVERGENCE_TAKESTEP"),
  //                            this,&DynamicBalancerMethod::solve);
}

//////////////////////////////////////////////////////////////////////////////

DynamicBalancerMethod::DynamicBalancerMethod(const std::string& name)
  : Method(name)
{
  // define which functions might be called dynamic
  build_dynamic_functions();
  // regist which functions might be called by raising Events
  registActionListeners();
  // regist the configuration options
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

DynamicBalancerMethod::~DynamicBalancerMethod()
{
}

//////////////////////////////////////////////////////////////////////////////

void DynamicBalancerMethod::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  Method::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DynamicBalancerMethod::doDynamicBalance()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  doDynamicBalanceImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void DynamicBalancerMethod::setMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void DynamicBalancerMethod::unsetMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
