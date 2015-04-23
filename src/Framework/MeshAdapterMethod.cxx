// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshAdapterMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void MeshAdapterMethod::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("AdaptRate","Rate to adapt the mesh.");
}

//////////////////////////////////////////////////////////////////////////////

void MeshAdapterMethod::build_dynamic_functions()
{
  add_dynamic_function("adaptMesh",&MeshAdapterMethod::adaptMesh);
  add_dynamic_function("remesh",&MeshAdapterMethod::remesh);
}

//////////////////////////////////////////////////////////////////////////////

void MeshAdapterMethod::registActionListeners()
{
  Method::registActionListeners();
  
  // const std::string ssname = SubSystemStatusStack::getCurrentName();   
  // event_handler->addListener(event_handler->key(ssname, "CF_ON_MESHADAOTER_ADAPT"),
  //                            this,&MeshAdapterMethod::solve);
}

//////////////////////////////////////////////////////////////////////////////

MeshAdapterMethod::MeshAdapterMethod(const std::string& name)
  : Method(name)
{
  // define which functions might be called dynamic
  build_dynamic_functions();
  // regist which functions might be called by raising Events
  registActionListeners();
  // regist the configuration options
  addConfigOptionsTo(this);

  _adaptRate = 1;
  setParameter("AdaptRate",&_adaptRate);
}

//////////////////////////////////////////////////////////////////////////////

MeshAdapterMethod::~MeshAdapterMethod()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshAdapterMethod::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  Method::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void MeshAdapterMethod::adaptMesh()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  adaptMeshImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void MeshAdapterMethod::remesh()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  remeshImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void MeshAdapterMethod::setMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshAdapterMethod::unsetMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
