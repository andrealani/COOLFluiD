// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "DaedalusMeshGenerator.hh"
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

MethodCommandProvider<DaedalusMeshGenerator,
                      SimpleMeshAdapterData,
                      SimpleGlobalMeshAdapterModule>
                      daedalusMeshGeneratorProvider("DaedalusMeshGenerator");

//////////////////////////////////////////////////////////////////////////////

DaedalusMeshGenerator::DaedalusMeshGenerator(const std::string& name)  :
  SimpleMeshAdapterCom(name)
{
}

//////////////////////////////////////////////////////////////////////////////

void DaedalusMeshGenerator::execute()
{
  CFAUTOTRACE;

  //Run mesh generation using the startfile and cfg and pts files created in the prepare phase
  std::string startDaedalus;
  startDaedalus = "./dae";

  if (PE::GetPE().IsParallel()) {

    PE::GetPE().setBarrier();

    if (PE::GetPE().GetRank () == 0) {
      Common::OSystem::getInstance().executeCommand(startDaedalus);
    }

    PE::GetPE().setBarrier();
  }
  else{
    Common::OSystem::getInstance().executeCommand(startDaedalus);
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
