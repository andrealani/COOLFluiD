// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/ProcessInfo.hh"
#include "Common/OSystem.hh"
#include "Common/BadValueException.hh"

#include "Framework/MeshCreator.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/MethodData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void MeshCreator::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshCreator::build_dynamic_functions()
{
  add_dynamic_function("generateMeshData",&MeshCreator::generateMeshData);
  add_dynamic_function("processMeshData",&MeshCreator::processMeshData);
  add_dynamic_function("buildMeshData",&MeshCreator::buildMeshData);
}

//////////////////////////////////////////////////////////////////////////////

void MeshCreator::registActionListeners()
{
  Method::registActionListeners();

  // const std::string ssname = SubSystemStatusStack::getCurrentName(); 
  // event_handler->addListener(event_handler->key(ssname, "CF_ON_COUPLERMETHOD_MATCH"),
  //                            this,&MeshCreator::match);
}

//////////////////////////////////////////////////////////////////////////////

MeshCreator::MeshCreator(const std::string& name) : Method(name)
{
  // define which functions might be called dynamic
  build_dynamic_functions();
  // regist which functions might be called by raising Events
  registActionListeners();
  // regist the configuration options
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

MeshCreator::~MeshCreator()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshCreator::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  Method::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void MeshCreator::generateMeshData()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  CFLog(NOTICE,"-------------------------------------------------------------\n");
  CFLog(NOTICE,"MeshCreator [" << getName() << "] Generate or Read Mesh\n");
  CFLog(NOTICE,"\nMemory usage before building mesh: " << Common::OSystem::getInstance().getProcessInfo()->memoryUsage() << "\n\n");

  generateMeshDataImpl();

  CFLog(NOTICE,"\nMemory usage after building mesh: " << Common::OSystem::getInstance().getProcessInfo()->memoryUsage() << "\n");
  CFLog(NOTICE,"-------------------------------------------------------------\n");

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void MeshCreator::processMeshData()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  /// @todo this could be a post generation hook not directly accessible from
  ///       the interface, maybe controled by commands
  processMeshDataImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void MeshCreator::buildMeshData()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  CFLog(NOTICE,"-------------------------------------------------------------\n");
  CFLog(NOTICE,"MeshCreator [" << getName() << "] Building Mesh Data\n");

  buildMeshDataImpl();

  CFLog(NOTICE,"-------------------------------------------------------------\n");

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void MeshCreator::setMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshCreator::unsetMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

