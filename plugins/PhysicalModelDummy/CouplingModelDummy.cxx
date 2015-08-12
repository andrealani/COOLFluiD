// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "PhysicalModelDummy/Dummy.hh"
#include "PhysicalModelDummy/CouplingModelDummy.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace PhysicalModelDummy {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CouplingModelDummy,
			    PhysicalModelImpl,
			    DummyModule,
			    1 >
dummyCouplingModelProvider("CouplingModelDummy");

//////////////////////////////////////////////////////////////////////////////

CouplingModelDummy::CouplingModelDummy(const std::string& name) :
  PhysicalModelDummy(name)
{
  CFAUTOTRACE;
  
  addConfigOptionsTo(this);
  
  m_sendIDs = std::vector<CFuint>();
  setParameter("SendIDs", &m_sendIDs);
  
  m_recvIDs = std::vector<CFuint>();
  setParameter("RecvIDs", &m_recvIDs);
}
    
//////////////////////////////////////////////////////////////////////////////

CouplingModelDummy::~CouplingModelDummy()
{
}

//////////////////////////////////////////////////////////////////////////////

void CouplingModelDummy::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector< CFuint > >
    ("SendIDs","Vector storing IDs of the variables when sent.");
  
  options.addConfigOption< std::vector< CFuint > >
    ("RecvIDs","Vector storing IDs of the variables when received.");
}
    
//////////////////////////////////////////////////////////////////////////////

void CouplingModelDummy::configure ( Config::ConfigArgs& args )
{
  PhysicalModelDummy::configure(args);
 
  if (m_sendIDs.size() == 0) {
    CFLog(ERROR, "CouplingModelDummy::configure() => SendIDs not specified!");
    cf_assert(m_sendIDs.size() > 0);
  }
  
  if (m_recvIDs.size() == 0) {
    CFLog(ERROR, "CouplingModelDummy::configure() => RecvIDs not specified!");
    cf_assert(m_recvIDs.size() > 0);
  }

  if (m_sendIDs.size() != m_recvIDs.size()) {
    CFLog(ERROR, "CouplingModelDummy::configure() => SendIDs and RecvIDs must have same size!");
    cf_assert(m_sendIDs.size() == m_recvIDs.size());
  }
}
    
//////////////////////////////////////////////////////////////////////////////

  }  // namespace PhysicalModelDummy
}  // namespace COOLFluiD

