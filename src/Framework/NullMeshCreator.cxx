// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NullMeshCreator.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullMeshCreator,
               MeshCreator,
               FrameworkLib,
               1>
nullMeshCreatorProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullMeshCreator::NullMeshCreator(const std::string& name)
  : MeshCreator(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullMeshCreator::~NullMeshCreator()
{
}

//////////////////////////////////////////////////////////////////////////////

void NullMeshCreator::setMethodImpl()
{
  MeshCreator::setMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void NullMeshCreator::unsetMethodImpl()
{
  MeshCreator::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> NullMeshCreator::getMethodData () const
{
  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

void NullMeshCreator::generateMeshDataImpl()
{
  CFLog(VERBOSE,"NullMeshCreator::generateMeshData() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullMeshCreator::processMeshDataImpl()
{
  CFLog(VERBOSE,"NullMeshCreator::processMeshData() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullMeshCreator::buildMeshDataImpl()
{
  CFLog(VERBOSE,"NullMeshCreator::buildMeshData() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullMeshCreator::modifyFileNameForRestart(const std::string filename)
{
  CFLog(VERBOSE,"NullMeshCreator::modifyFileNameForRestart() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullMeshCreator::modifyFileName(const std::string filename)
{
  CFLog(VERBOSE,"NullMeshCreator::modifyFileName() called!" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
