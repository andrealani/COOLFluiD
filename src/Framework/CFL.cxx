// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/CFL.hh"
#include "Environment/Factory.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

CFL::CFL() :
  OwnedObject(),
  ConfigObject("CFL"),
  m_value(1.0),
  m_computeCFLStr()
{
   addConfigOptionsTo(this);

  m_value = 1.0;
  setParameter("Value",&m_value);

  m_computeCFLStr = "Null";
  setParameter("ComputeCFL", &m_computeCFLStr);
}

//////////////////////////////////////////////////////////////////////////////

void CFL::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption<CFreal, Config::DynamicOption<> > ("Value","CFL number (it can be changed interactively).");
  options.addConfigOption<std::string>("ComputeCFL","CFL calculator.");
}

//////////////////////////////////////////////////////////////////////////////

CFL::~CFL()
{
}

//////////////////////////////////////////////////////////////////////////////

void CFL::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);

  // configure the CFL computer
  Common::SafePtr<ComputeCFL::PROVIDER> prov =
    FACTORY_GET_PROVIDER(getFactoryRegistry(), ComputeCFL, m_computeCFLStr);
  m_computeCFL.reset(prov->create(prov->getName()));
  
  configureNested ( m_computeCFL.getPtr(), args );
  
  (*m_computeCFL).setCFL(this);
}
    
//////////////////////////////////////////////////////////////////////////////

void CFL::update(Framework::ConvergenceStatus * cstatus)
{
  if (cstatus == CFNULL) {
    Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
    cstatus = new Framework::ConvergenceStatus;
    *cstatus = subSysStatus->getConvergenceStatus();
  }
  (*m_computeCFL)(*cstatus);
}

//////////////////////////////////////////////////////////////////////////////

void CFL::update()
{
   (*m_computeCFL)(SubSystemStatusStack::getActive()->getConvergenceStatus());
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

