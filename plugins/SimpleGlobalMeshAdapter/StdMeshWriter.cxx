// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "StdMeshWriter.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/MeshCreator.hh"
#include "Framework/OutputFormatter.hh"
#include "Common/ProcessInfo.hh"
#include "Framework/SimulationStatus.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdMeshWriter, SimpleMeshAdapterData, SimpleGlobalMeshAdapterModule> StdMeshWriterProvider("StdMeshWriter");

//////////////////////////////////////////////////////////////////////////////

StdMeshWriter::StdMeshWriter(const std::string& name)  :
  SimpleMeshAdapterCom(name)
{
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshWriter::execute()
{
  CFAUTOTRACE;

  MultiMethodHandle<OutputFormatter> outputFormat = getMethodData().getOutputFormatter();
  cf_assert(outputFormat.size() == 1);

  const std::string subsystemName = SubSystemStatusStack::getActive()->getSubSystemName();
  std::string nspName = outputFormat[0]->getNamespace();
  std::string filename = "INTERPOLATED_"+ getMethodData().getAdaptedMeshFileName() ;

  outputFormat[0]->setMethod();

  outputFormat[0]->setOutputFileName(filename);
  outputFormat[0]->open ();
  outputFormat[0]->write();
  outputFormat[0]->close();

  std::string fullOutputName = SimulationStatus::getInstance().getLastOutputFile(subsystemName, getMethodData().getNamespace()).string();

///@todo here the other mesh should be in the same subsystem
  SimulationStatus::getInstance().setLastOutputFile(
        fullOutputName,getMethodData().getOtherNamespace(),subsystemName, true);

  outputFormat[0]->unsetMethod();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
