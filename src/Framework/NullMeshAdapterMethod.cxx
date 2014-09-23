// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NullMeshAdapterMethod.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullMeshAdapterMethod,
               MeshAdapterMethod,
               FrameworkLib,
               1>
nullMeshAdapterMethodProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullMeshAdapterMethod::NullMeshAdapterMethod(const std::string& name)
  : MeshAdapterMethod(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullMeshAdapterMethod::~NullMeshAdapterMethod()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> NullMeshAdapterMethod::getMethodData () const
{
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

void NullMeshAdapterMethod::adaptMeshImpl()
{
  CFLogDebugMed("NullMeshAdapterMethod::adaptMesh() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullMeshAdapterMethod::setMethodImpl()
{
  MeshAdapterMethod::setMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullMeshAdapterMethod::unsetMethodImpl()
{
  MeshAdapterMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
