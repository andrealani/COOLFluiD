// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "DummyMeshGenerator.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<DummyMeshGenerator, SimpleMeshAdapterData, SimpleGlobalMeshAdapterModule> dummyMeshGeneratorProvider("DummyMeshGenerator");

//////////////////////////////////////////////////////////////////////////////

void DummyMeshGenerator::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("Filename","Name of the input file.");
}

//////////////////////////////////////////////////////////////////////////////

DummyMeshGenerator::DummyMeshGenerator(const std::string& name)  :
  SimpleMeshAdapterCom(name)
{
   addConfigOptionsTo(this);

   _filename = "";
   setParameter("Filename",&_filename);

}

//////////////////////////////////////////////////////////////////////////////

void DummyMeshGenerator::execute()
{
  CFAUTOTRACE;

  cf_assert(_filename != "");
  getMethodData().setAdaptedMeshFileName(_filename);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
