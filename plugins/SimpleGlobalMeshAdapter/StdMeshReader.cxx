// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "StdMeshReader.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/MeshCreator.hh"
#include "Common/ProcessInfo.hh"
#include "Common/OSystem.hh"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdMeshReader, SimpleMeshAdapterData, SimpleGlobalMeshAdapterModule> stdMeshReaderProvider("StdMeshReader");

//////////////////////////////////////////////////////////////////////////////

StdMeshReader::StdMeshReader(const std::string& name)  :
  SimpleMeshAdapterCom(name)
{
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshReader::execute()
{
  CFAUTOTRACE;

  MultiMethodHandle<MeshCreator> meshCreator = getMethodData().getMeshCreator();
  cf_assert(meshCreator.size() == 1);

  std::string subsystemName = getName();
  std::string nspName = meshCreator[0]->getNamespace();
  std::string filename = getMethodData().getAdaptedMeshFileName();

  meshCreator[0]->modifyFileName(filename);

  buildMeshData();
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshReader::buildMeshData()
{
//   CFLogInfo("Setting up MeshCreator\n");

  MultiMethodHandle<MeshCreator> meshCreator = getMethodData().getMeshCreator();
  cf_assert(meshCreator.size() == 1);

  meshCreator[0]->setMethod();
  
  CFLog(NOTICE,"Building MeshData's\n");
  CFLog(NOTICE,"Memory Usage before building mesh: " << Common::OSystem::getInstance().getProcessInfo()->memoryUsage() << "\n");
  
  //Create the CFMeshData for each Namespace
  meshCreator[0]->generateMeshData();
  
  //Process the CFMeshData to do:
  //  - renumbering
  //  - conversion FVM <-> FEM
  meshCreator[0]->processMeshData();
  
  //Use the CFMeshData to build the mesh
  meshCreator[0]->buildMeshData();
  
  CFLog(NOTICE,"Memory Usage after building mesh: " << Common::OSystem::getInstance().getProcessInfo()->memoryUsage() << "\n");

  vector <Common::SafePtr<MeshData> > meshDataVector = MeshDataStack::getInstance().getAllEntries();

  bool found = false;
  CFuint meshDataID(0);
  for(CFuint iMeshData = 0; iMeshData < meshDataVector.size(); iMeshData++)
  {
    if(meshCreator[0]->getNamespace() == meshDataVector[iMeshData]->getPrimaryNamespace())
    {
      meshDataID = iMeshData;
      found = true;
    }
  }
  cf_assert(found);

  CFLogNotice("Building TRS info for MeshData in Namespace " << meshDataVector[meshDataID]->getPrimaryNamespace() << "\n");

  // === conversion fix ===
  vector<std::string> & TotalTRSNames =
    meshDataVector[meshDataID]->getTotalTRSNames ();
  vector<vector<CFuint> > & TotalTRSInfo =
    meshDataVector[meshDataID]->getTotalTRSInfo ();

  if (TotalTRSInfo.empty ()) {
    cf_assert(TotalTRSNames.empty());
    vector< SafePtr<TopologicalRegionSet> > trsList =
      meshDataVector[meshDataID]->getTrsList();

    // count in advance the number of writable TRS
    CFuint sizeTRSToWrite = 0;
    for (CFuint i = 0; i < trsList.size(); ++i) {
      if (trsList[i]->hasTag("writable")) {
        sizeTRSToWrite++;
      }
    }
    cf_assert(sizeTRSToWrite < trsList.size());

    TotalTRSNames.resize(sizeTRSToWrite);
    TotalTRSInfo.resize(sizeTRSToWrite);
    CFuint counter = 0;
    vector< Common::SafePtr<TopologicalRegionSet> >::iterator it;
    for (it = trsList.begin(); it != trsList.end(); ++it)
    {
      if ((*it)->hasTag("writable"))
      {
        TotalTRSNames[counter] = (*it)->getName();
        const CFuint nbTRs = (*it)->getNbTRs();
        TotalTRSInfo[counter].resize(nbTRs);
        for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
          TotalTRSInfo[counter][iTR] = (*(*it))[iTR]->getLocalNbGeoEnts();
          CFLog(VERBOSE, "TRS : " << (*it)->getName()
                                << ", TR : "<< iTR
                                << ", nbGeos : "
                                << TotalTRSInfo[counter][iTR] << "\n");
        }
        counter++;
      }
    }

    cf_assert(counter == sizeTRSToWrite);
  }

  if (PE::GetPE().IsParallel()) {

      const std::string parStateVecName = meshDataVector[meshDataID]->getPrimaryNamespace() + "-states";
      const std::string parNodeVecName = meshDataVector[meshDataID]->getPrimaryNamespace() + "-nodes";

      DataHandle<State*, GLOBAL> states =
        meshDataVector[meshDataID]->getDataStorage()->getGlobalData<State*>(parStateVecName);

      DataHandle<Node*, GLOBAL> nodes =
        meshDataVector[meshDataID]->getDataStorage()->getGlobalData<Node*>(parNodeVecName);

      states.buildMap ();
      nodes.buildMap ();
    }

  meshCreator[0]->unsetMethod();
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
