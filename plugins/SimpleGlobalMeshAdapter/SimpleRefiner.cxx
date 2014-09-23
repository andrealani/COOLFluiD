// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Common/SelfRegistPtr.hh"
#include "Environment/DirPaths.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshFormatConverter.hh"

#include "SimpleGlobalMeshAdapter/SimpleRefiner.hh"
#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SimpleRefiner, SimpleMeshAdapterData, SimpleGlobalMeshAdapterModule>
SimpleRefinerProvider("SimpleRefiner");

//////////////////////////////////////////////////////////////////////////////

void SimpleRefiner::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("Filename","Name of the input file.");
  options.addConfigOption< std::string >("Refiner","Name of the refiner object.");
}

//////////////////////////////////////////////////////////////////////////////

SimpleRefiner::SimpleRefiner(const std::string& name)  :
  SimpleMeshAdapterCom(name)
{
  addConfigOptionsTo(this);

  _filename = "";
  setParameter("Filename",&_filename);

  _refinerStr = "";
  setParameter("Refiner",&_refinerStr);
}

//////////////////////////////////////////////////////////////////////////////

SimpleRefiner::~SimpleRefiner() {}

//////////////////////////////////////////////////////////////////////////////

void SimpleRefiner::execute()
{
  CFAUTOTRACE;

  std::string outputFileName = "REFINED" + _filename;
  getMethodData().setAdaptedMeshFileName(outputFileName);

  // build the MeshFormatConverterConvertFrom
  SelfRegistPtr<MeshFormatConverter> refiner =
    Environment::Factory<MeshFormatConverter>::getInstance().getProvider(_refinerStr)->create(_refinerStr);
  configureNested ( refiner.getPtr(), m_stored_args );

  using namespace boost::filesystem;

  path fromfile = Environment::DirPaths::getInstance().getWorkingDir() / path ( _filename );
  path tofile   = Environment::DirPaths::getInstance().getWorkingDir() / path ( outputFileName );

  // refiner works serially on the processor with rank == 0
  if (PE::GetPE().IsParallel())
  {
    PE::GetPE().setBarrier();

    if (PE::GetPE().GetRank() == 0) {

      refiner->convert(fromfile, tofile);
    }

    PE::GetPE().setBarrier();
  }
  else
  {
    cf_assert(!PE::GetPE().IsParallel());

    refiner->convert(fromfile, tofile);
  }

}

//////////////////////////////////////////////////////////////////////////////

void SimpleRefiner::configure ( Config::ConfigArgs& args )
{
  SimpleMeshAdapterCom::configure ( args );

  m_stored_args = args; // store the configuration args to use later in MeshFormatConverter
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
