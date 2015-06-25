// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "GambitMeshGenerator.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/OSystem.hh"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<GambitMeshGenerator, SimpleMeshAdapterData, SimpleGlobalMeshAdapterModule> GambitMeshGeneratorProvider("GambitMeshGenerator");

//////////////////////////////////////////////////////////////////////////////

void GambitMeshGenerator::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("JournalFile","Name of the journal file.");
  options.addConfigOption< std::string >("OutputFile","Name of the output mesh.");
}

//////////////////////////////////////////////////////////////////////////////

GambitMeshGenerator::GambitMeshGenerator(const std::string& name)  :
  SimpleMeshAdapterCom(name)
{
   addConfigOptionsTo(this);

   _filename = "";
   setParameter("OutputFile",&_filename);

   _journalFile = "";
   setParameter("JournalFile",&_journalFile);


}

//////////////////////////////////////////////////////////////////////////////

void GambitMeshGenerator::execute()
{
  CFAUTOTRACE;

  //Set output filename
  cf_assert(_filename != "");
  getMethodData().setAdaptedMeshFileName(_filename);

  //Run mesh generation using the journal file created in the prepare phase
  std::string startGambit;
  startGambit = "gambit -inp " + _journalFile;
  
  const std::string nsp = getMethodData().getNamespace();
  
  if (PE::GetPE().IsParallel()) {

    PE::GetPE().setBarrier(nsp);

    if (PE::GetPE().GetRank (nsp) == 0) {
      Common::OSystem::getInstance().executeCommand(startGambit);
    }

    PE::GetPE().setBarrier(nsp);
  }
  else{
    Common::OSystem::getInstance().executeCommand(startGambit);
  }
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
