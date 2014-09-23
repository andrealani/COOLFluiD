// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "SimpleMeshAdapterData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModelImpl.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

void SimpleMeshAdapterData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption<std::string>("OtherNamespace","Name of the other namespace.");
}

//////////////////////////////////////////////////////////////////////////////

SimpleMeshAdapterData::SimpleMeshAdapterData(Common::SafePtr<Framework::Method> owner)
  : MeshAdapterData(owner),
    _isRemeshNeeded(true)
{
   addConfigOptionsTo(this);

   setParameter("OtherNamespace",&_otherNamespace);
}

//////////////////////////////////////////////////////////////////////////////

SimpleMeshAdapterData::~SimpleMeshAdapterData()
{
}

//////////////////////////////////////////////////////////////////////////////

void SimpleMeshAdapterData::configure ( Config::ConfigArgs& args )
{
  MeshAdapterData::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void SimpleMeshAdapterData::setup()
{
  MeshAdapterData::setup();
}

//////////////////////////////////////////////////////////////////////////////

void SimpleMeshAdapterData::unsetup()
{
  MeshAdapterData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

