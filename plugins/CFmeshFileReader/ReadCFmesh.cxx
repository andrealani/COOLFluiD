// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <set>

#include "Common/Table.hh"
#include "Common/Stopwatch.hh"
#include "Common/ProcessInfo.hh"

#include "Environment/ConcreteProvider.hh"
#include "Environment/DirPaths.hh"

#include "Common/BadValueException.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshDataBuilder.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/MethodData.hh"

#include "CFmeshFileReader/CFmeshFileReader.hh"
#include "CFmeshFileReader/ReadCFmesh.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ReadCFmesh, CFmeshReaderData, CFmeshFileReaderPlugin>
stdReadCFmeshProvider("StdReadCFmesh");

//////////////////////////////////////////////////////////////////////////////

ReadCFmesh::ReadCFmesh(const std::string& name) :
  ReadBase(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ReadCFmesh::~ReadCFmesh()
{
}

//////////////////////////////////////////////////////////////////////////////

void ReadCFmesh::configure ( Config::ConfigArgs& args )
{
  ReadBase::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ReadCFmesh::execute()
{
  CFAUTOTRACE;
  
  m_data->setDataSockets(socket_states, socket_nodes, &m_sockets);
  m_data->setExtraVarNamesAndTags(getMethodData().getExtraNVarSocketNamesAndTags(),
                                  getMethodData().getExtraSVarSocketNamesAndTags(),
                                  getMethodData().getExtraVarSocketNamesAndTags());
  SafePtr<CFmeshReaderSource> ptr = m_data.get();
  cf_assert(ptr.isNotNull());
  
  m_reader.setFactoryRegistry(getFactoryRegistry());
  m_reader.setReadData(ptr);
  m_reader.setStateInitValues(getMethodData().getUseInitValues(),
            getMethodData().getInitValues(),
            getMethodData().getInitValuesIDs());


  Stopwatch<WallTime> stp;
  stp.start();

  try
  {
    m_reader.readFromFile(Environment::DirPaths::getInstance().getWorkingDir()
      / getMethodData().getFileName());
  }
  catch (Common::Exception& e)
  {
    CFLog( ERROR, e.what() << CFendl );
    throw;
  }

  CFLog(INFO,"Reading data from "
  << getMethodData().getFileName().string()
  << " took " << stp.read() << "s\n");

  CFLog(NOTICE,"Memory Usage after mesh reading: " << Common::OSystem::getInstance().getProcessInfo()->memoryUsage() << "\n");

  if (m_data->getDimension() != PhysicalModelStack::getActive()->getDim())
    throw BadValueException (FromHere(),"Mesh has wrong dimensionality.");

  // if the file contained states that do not match the physical
  // model then we need to correct their size
  correctStates();

  // builder of the basic data in MeshData
  // dont forget to release in the end
  Common::SelfRegistPtr<MeshDataBuilder> meshDataBuilder =
    getMethodData().getCollaborator<SpaceMethod>()->createMeshBuilder();

  DataHandle<State*,GLOBAL> states = socket_states.getDataHandle();

  for (CFuint i = 0; i < states.size(); ++i){
    const CFuint globalID = states[i]->getLocalID();
    states[i]->setGlobalID(globalID);
  }

  /// @todo AL : hack to restore the conversion to hdf5
  /// this will not work correctly in parallel
  // Nodes are here adimensionalized with a reference length
    DataHandle<Node*,GLOBAL> nodes = socket_nodes.getDataHandle();
  for (CFuint i = 0; i < nodes.size(); ++i){
    const CFuint globalID = nodes[i]->getLocalID();
    nodes[i]->setGlobalID(globalID);
  }

  meshDataBuilder->setCFmeshData(m_data.get());

  applyScalings();
  applyTranslation();

  // get all useful global info, store them in MeshData for later use
  // the following are local values (not really global) :
  // this works only serial

  const CFuint totalNbNodes    = nodes.size();
  const CFuint totalNbStates   = states.size();
  const CFuint nbElemTypes     = m_data->getNbElementTypes();

  cf_assert(m_data->getNbElements() > 0);
  cf_assert(nbElemTypes > 0);
  
  SafePtr<vector<ElementTypeData> > eData = m_data->getElementTypeData();
  cf_assert(nbElemTypes == eData->size());
  
  // element type info
  vector<CFuint> elementTypeCount(nbElemTypes);
  
  // AL: should I reset element type data in MeshData  here??
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
    elementTypeCount[iType] = (*eData)[iType].getNbElems();
    (*eData)[iType].setNbTotalElems((*eData)[iType].getNbElems());
  }
  
  cf_assert(totalNbStates > 0);
  cf_assert(totalNbNodes > 0);
  cf_assert(elementTypeCount.size() > 0);

  CFLog(VERBOSE, "ReadCFmesh:: Total nodes: " << totalNbNodes
  << ", states: "  << totalNbStates
  << ", cells: "   << m_data->getNbElements() << "\n");

  // set the total number of nodes, states and elements in the MeshData
  // to make it available later while writing
  MeshDataStack::getActive()->setTotalNodeCount(totalNbNodes);
  MeshDataStack::getActive()->setTotalStateCount(totalNbStates);
  MeshDataStack::getActive()->setTotalElementCount(elementTypeCount);
  
  // mesh data builder creates data that are IO-dependent
  meshDataBuilder->computeGeoTypeInfo();
  meshDataBuilder->createTopologicalRegionSets();

  CFLog(VERBOSE,"Memory Usage after TRS: " << Common::OSystem::getInstance().getProcessInfo()->memoryUsage() << "\n");

  //set some global mesh values useful in many places
  meshDataBuilder->setMaxGlobalInfo();

  CFLog(INFO,"Building MeshData from "
  << getMethodData().getFileName().string()
  << " took "
  << stp.read()
  <<"s\n");

  // deallocate the unnecessary memory
  m_data->releaseMemory();
  meshDataBuilder->releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD
