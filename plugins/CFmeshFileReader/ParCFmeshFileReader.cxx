// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <numeric>

#include <boost/progress.hpp>

#include "Common/PE.hh"
#include "Common/CFPrintContainer.hh"
#include "Common/ProcessInfo.hh"
#include "Common/OSystem.hh"
#include "Common/StringOps.hh"
#include "Common/SwapEmpty.hh"
#include "Common/BadValueException.hh"

#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/NamespaceSwitcher.hh"
#include "Framework/MeshData.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/MeshPartitioner.hh"
#include "Framework/SubSystemStatus.hh"

#include "CFmeshFileReader/ParCFmeshFileReader.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

ParCFmeshFileReader::ParCFmeshFileReader() :
  FileReader(),
  ConfigObject("ParCFmeshFileReader"),
  m_mapString2Reader(),
  m_partitioner(),
  m_comm(0),
  m_local_elem(CFNULL),
  m_myRank(0),
  m_nbProc(0),
  m_startNodeList(0),
  m_startStateList(0),
  m_useInitValues(),
  m_initValues(),
  m_initValuesIDs(),
  m_readData(CFNULL),
  m_originalNbEqs(0),
  m_totNbElem(0),
  m_totNbElemTypes(0),
  m_totNbNodes(0),
  m_totNbStates(0),
  m_nbElemPerProc(0),
  m_partitionerOutData(),
  m_localNodeIDs(),
  m_localStateIDs(),
  m_ghostNodeIDs(),
  m_ghostStateIDs(),
  m_mapGlobToLocNodeID(),
  m_mapGlobToLocStateID(),
  m_mapNodeElemID(),
  m_localElemIDs(),
  m_hasPastNodes(false),
  m_hasPastStates(false),
  m_hasInterNodes(false),
  m_hasInterStates(false)
{
  addConfigOptionsTo(this);

  m_nbOverLayers = 1;
  setParameter("NbOverlapLayers",&m_nbOverLayers);

#ifdef CF_HAVE_PARMETIS
  m_partitionerName = "ParMetis";
#else
  m_partitionerName = "Dumb";
#endif
  setParameter("Partitioner",&m_partitionerName);

  m_merge_trs = vector<std::string>();
  setParameter("MergeTRS",&m_merge_trs);
  
  m_inputToUpdateVecStr = "Identity";
  setParameter("InputToUpdate",&m_inputToUpdateVecStr);
}

//////////////////////////////////////////////////////////////////////////////

ParCFmeshFileReader::~ParCFmeshFileReader()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("NbOverlapLayers", "Number of layers of overlap");

  options.addConfigOption< std::string >("Partitioner", "Mesh partitioner to use");

  options.addConfigOption< std::vector<std::string> > ("MergeTRS", "Topological regions sets to be merged");

  options.addConfigOption< std::string >("InputToUpdate", "Transformer from input to update variables");
}

/////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);

  configureTRSMerging();

  SafePtr<MeshPartitioner::PROVIDER> provider;
  try
  {
     provider = Environment::Factory<MeshPartitioner>::getInstance().
      getProvider(m_partitionerName);
  }
  catch (NoSuchValueException& e)
  {
	CFLog(VERBOSE, e.what() << "\n");
    std::string message = "Could not find mesh partitioner " + m_partitionerName;
    throw NoSuchValueException (FromHere(),message.c_str());
  }

  cf_assert(provider.isNotNull());
  m_partitioner = provider->create(m_partitionerName);
  cf_assert(m_partitioner.isNotNull());

  // configure the mesh partitioner
  configureNested ( m_partitioner.getPtr(), args );
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::configureTRSMerging()
{
  m_curr_trs = "";
  m_curr_nbtr = 0;
  m_trs_reduction = 0;
  m_rev_trs_merge.clear();
  m_trs_idxmap.clear();
  if (!m_merge_trs.empty())
  {
    for (CFuint i = 0; i < m_merge_trs.size(); ++i)
    {
      std::string opt = m_merge_trs[i];
      std::vector<std::string> wrds = Common::StringOps::getWords(opt,':');

      if(wrds.size() < 3)
        throw BadValueException (FromHere(),"MergeTRS badly defined, should be at least 3 TRS names separated by ':'\n");

      std::vector<std::string>* sum = new std::vector<std::string>();
      std::string newtrs = wrds[0];
      CFLog(INFO, "Merging into TRS " << newtrs << " following TRS's [ " );
      for (CFuint w = 1; w < wrds.size(); ++w)
      {
        m_rev_trs_merge[wrds[w]] = newtrs;
        sum->push_back(wrds[w]);
        m_trs_reduction++;
        CFLog(INFO, wrds[w] << " " );
      }
      CFLog(INFO, "]\n" );
      m_trs_reduction--; // one less to account for the one created
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::setup()
{
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  
  m_comm   = PE::GetPE().GetCommunicator(nsp);
  m_myRank = PE::GetPE().GetRank(nsp);
  m_nbProc = PE::GetPE().GetProcessorCount(nsp);
  
  // each processor allocates the elmdist array with size==m_nbProc+1
  m_nbElemPerProc.resize(m_nbProc);
  
  setMapString2Readers();
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::setMapString2Readers()
{
  m_mapString2Reader["!COOLFLUID_VERSION"]     = &ParCFmeshFileReader::readCFVersion;
  m_mapString2Reader["!COOLFLUID_SVNVERSION"]  = &ParCFmeshFileReader::readSvnVersion;
  m_mapString2Reader["!CFMESH_FORMAT_VERSION"] = &ParCFmeshFileReader::readCFmeshVersion;
  m_mapString2Reader["!NB_DIM"]             = &ParCFmeshFileReader::readDimension;
  m_mapString2Reader["!NB_EQ"]              = &ParCFmeshFileReader::readNbEquations;
  m_mapString2Reader["!NB_NODES"]           = &ParCFmeshFileReader::readNbNodes;
  m_mapString2Reader["!NB_STATES"]          = &ParCFmeshFileReader::readNbStates;
  m_mapString2Reader["!STORE_PASTSTATES"]   = &ParCFmeshFileReader::readStorePastStates;
  m_mapString2Reader["!STORE_PASTNODES"]    = &ParCFmeshFileReader::readStorePastNodes;
  m_mapString2Reader["!STORE_INTERSTATES"]   = &ParCFmeshFileReader::readStoreInterStates;
  m_mapString2Reader["!STORE_INTERNODES"]    = &ParCFmeshFileReader::readStoreInterNodes;
  m_mapString2Reader["!NB_EXTRA_SVARS"]     = &ParCFmeshFileReader::readNbExtraStateVars;
  m_mapString2Reader["!NB_EXTRA_NVARS"]     = &ParCFmeshFileReader::readNbExtraNodalVars;
  m_mapString2Reader["!NB_EXTRA_VARS"]     = &ParCFmeshFileReader::readNbExtraVars;
  m_mapString2Reader["!EXTRA_NVARS_NAMES"]  = &ParCFmeshFileReader::readExtraNodalVarNames;
  m_mapString2Reader["!EXTRA_SVARS_NAMES"]  = &ParCFmeshFileReader::readExtraStateVarNames;
  m_mapString2Reader["!EXTRA_VARS_NAMES"]  = &ParCFmeshFileReader::readExtraVarNames;
  m_mapString2Reader["!EXTRA_NVARS_STRIDES"]  = &ParCFmeshFileReader::readExtraNodalVarStrides;
  m_mapString2Reader["!EXTRA_SVARS_STRIDES"]  = &ParCFmeshFileReader::readExtraStateVarStrides;
  m_mapString2Reader["!EXTRA_VARS_STRIDES"]  = &ParCFmeshFileReader::readExtraVarStrides;
  m_mapString2Reader["!EXTRA_VARS"]         = &ParCFmeshFileReader::readExtraVars;
  m_mapString2Reader["!NB_ELEM"]            = &ParCFmeshFileReader::readNbElements;
  m_mapString2Reader["!NB_ELEM_TYPES"]      = &ParCFmeshFileReader::readNbElementTypes;
  m_mapString2Reader["!GEOM_POLYTYPE"]      = &ParCFmeshFileReader::readGeometricPolyType;
  m_mapString2Reader["!SOL_POLYTYPE"]       = &ParCFmeshFileReader::readSolutionPolyType;
  m_mapString2Reader["!GEOM_POLYORDER"]     = &ParCFmeshFileReader::readGeometricPolyOrder;
  m_mapString2Reader["!SOL_POLYORDER"]      = &ParCFmeshFileReader::readSolutionPolyOrder;
  m_mapString2Reader["!LIST_NODE"]          = &ParCFmeshFileReader::readNodeList;
  m_mapString2Reader["!LIST_STATE"]         = &ParCFmeshFileReader::readStateList;
  m_mapString2Reader["!NB_TRSs"]            = &ParCFmeshFileReader::readNbTRSs;
  m_mapString2Reader["!TRS_NAME"]           = &ParCFmeshFileReader::readTRSName;
  m_mapString2Reader["!NB_TRs"]             = &ParCFmeshFileReader::readNbTRs;
  m_mapString2Reader["!NB_GEOM_ENTS"]       = &ParCFmeshFileReader::readNbGeomEnts;
  m_mapString2Reader["!GEOM_TYPE"]          = &ParCFmeshFileReader::readGeomType;
  m_mapString2Reader["!LIST_GEOM_ENT"]      = &ParCFmeshFileReader::readGeomEntList;
  m_mapString2Reader["!ELEM_TYPES"]         = &ParCFmeshFileReader::readElementTypes;
  m_mapString2Reader["!NB_ELEM_PER_TYPE"]   = &ParCFmeshFileReader::readNbElementsPerType;
  m_mapString2Reader["!NB_NODES_PER_TYPE"]  = &ParCFmeshFileReader::readNbNodesPerType;
  m_mapString2Reader["!NB_STATES_PER_TYPE"] = &ParCFmeshFileReader::readNbStatesPerType;
  m_mapString2Reader["!NB_GROUPS"]          = &ParCFmeshFileReader::readNbGroups;
  m_mapString2Reader["!GROUP_NAME"]         = &ParCFmeshFileReader::readGroupName;
  m_mapString2Reader["!GROUP_ELEM_NB"]      = &ParCFmeshFileReader::readGroupElementNb;
  m_mapString2Reader["!GROUP_ELEM_LIST"]    = &ParCFmeshFileReader::readGroupElementList;
  m_mapString2Reader["!LIST_ELEM"]          = &ParCFmeshFileReader::readElementList;

//  std::string fileRoot = "envFile";
//  std::string fileProc = fileRoot +  StringOps::to_str(PE::GetPE().GetRank ());
//  std::string commandEnv =  "printenv >& " + fileProc;
//  Common::OSystem::getInstance().executeCommand(commandEnv);

}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readCFVersion(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readCFVersion() start\n");

  std::string version;
  fin >> version;
  CFLogDebugMin( "CF version : " << version << "\n");

  CFLogDebugMin( "ParCFmeshFileReader::readCFVersion() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readSvnVersion(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readSvnVersion() start\n");

  std::string version;
  fin >> version;
  CFLogDebugMin( "svn version : " << version << "\n");

  CFLogDebugMin( "ParCFmeshFileReader::readSvnVersion() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readCFmeshVersion(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readCFmeshVersion() start\n");

  std::string version;
  fin >> version;
  CFLogDebugMin( "CFmesh version : " << version << "\n");

  CFLogDebugMin( "ParCFmeshFileReader::readCFmeshVersion() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readDimension(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readDimension() start\n");

  CFint dimension = 0;
  fin >> dimension;

  getReadData().setDimension(static_cast<CFuint>(dimension));

  if(getReadData().getDimension() != DIM_1D &&
     getReadData().getDimension() != DIM_2D &&
     getReadData().getDimension() != DIM_3D)
  {
    throw BadFormatException (FromHere(),"Wrong dimension in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshFileReader::readDimension() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbEquations(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbEquations() start\n");

  if(m_useInitValues.size() == m_initValues.size()) {
    if (m_useInitValues.size() == m_initValuesIDs.size() && m_useInitValues.size() != 0) {
      throw BadValueException
  (FromHere(), "Only one between InitValues and InitValueIDs must have same size as UseInitValues vector");
    }
  }

  if(m_useInitValues.size() == m_initValuesIDs.size()) {
    if (m_useInitValues.size() == m_initValues.size() && m_useInitValues.size() != 0) {
      throw BadValueException
  (FromHere(), "Only one between InitValues and InitValueIDs must have same size as UseInitValues vector");
    }
  }

  fin >> m_originalNbEqs;
  CFuint nbEquations = m_originalNbEqs;
  cf_assert(m_originalNbEqs > 0);

  // the following is useful if restarting from a file with less equations
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  if (m_useInitValues.size() > static_cast<CFuint>(nbEquations)) {
    if (m_useInitValues.size() != nbEqs) {
      throw BadValueException (FromHere(),"UseInitValues vector has size different from nbEqs");
    }
  }

  if (m_originalNbEqs != nbEqs) {
    nbEquations = nbEqs;
  }

  getReadData().setNbEquations(static_cast<CFuint>(nbEquations));

  if(nbEquations < 1 || nbEquations != nbEqs) {
    throw BadFormatException (FromHere(),"Wrong number of equations in CFmesh");
  }

  CFLogNotice("Original NbEquations = " << m_originalNbEqs << "\n");
  CFLogNotice("Final NbEquations    = " << nbEquations << "\n");

  CFLogDebugMin( "ParCFmeshFileReader::readNbEquations() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbNodes(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbNodes() start\n");

  CFint nbNonUpdatableNodes = 0;
  fin >> m_totNbNodes >> nbNonUpdatableNodes;

  // set the total number of nodes in the MeshData
  MeshDataStack::getActive()->setTotalNodeCount(m_totNbNodes);

  getReadData().setNbUpdatableNodes
    (static_cast<CFuint>(m_totNbNodes));

  getReadData().setNbNonUpdatableNodes
    (static_cast<CFuint>(nbNonUpdatableNodes));

  if(m_totNbNodes < 1) {
    throw BadFormatException (FromHere(),"Number of nodes < 1 in CFmesh");
  }

  if(nbNonUpdatableNodes != 0) {
    throw BadFormatException (FromHere(),"Number of non updatable nodes != 0 in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshFileReader::readNbNodes() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbStates(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbStates() start\n");

  CFint nbNonUpdatableStates = 0;
  fin >> m_totNbStates >> nbNonUpdatableStates;

  // set the total number of states in the MeshData
  MeshDataStack::getActive()->setTotalStateCount(m_totNbStates);

  getReadData().setNbUpdatableStates
    (static_cast<CFuint>(m_totNbStates));

  getReadData().setNbNonUpdatableStates
    (static_cast<CFuint>(nbNonUpdatableStates));

  if(m_totNbStates < 1) {
    throw BadFormatException (FromHere(),"Number of states < 1 in CFmesh");
  }

  if(nbNonUpdatableStates != 0) {
    throw BadFormatException (FromHere(),"Number of non updatable states != 0 in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshFileReader::readNbStates() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readStorePastStates(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readStorePastStates() start\n");

  fin >> m_hasPastStates;

  if(getReadData().storePastStates())
  {
    if(m_hasPastStates == false){
      throw BadFormatException (FromHere(),"PastStates asked but missing in CFmesh");
    }
  }

  CFLogDebugMin( "ParCFmeshFileReader::readStorePastStates() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readStorePastNodes(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readStorePastNodes() start\n");

  fin >> m_hasPastNodes;

  if(getReadData().storePastNodes())
  {
    if(m_hasPastNodes == false){
      throw BadFormatException (FromHere(),"PastNodes asked but missing in CFmesh");
    }
  }

  CFLogDebugMin( "ParCFmeshFileReader::readStorePastNodes() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readStoreInterStates(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readStoreInterStates() start\n");

  fin >> m_hasInterStates;

  if(getReadData().storeInterStates())
  {
    if(m_hasInterStates == false){
      throw BadFormatException (FromHere(),"InterStates asked but missing in CFmesh");
    }
  }

  CFLogDebugMin( "ParCFmeshFileReader::readStoreInterStates() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readStoreInterNodes(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readStoreInterNodes() start\n");

  fin >> m_hasInterNodes;

  if(getReadData().storeInterNodes())
  {
    if(m_hasInterNodes == false){
      throw BadFormatException (FromHere(),"InterNodes asked but missing in CFmesh");
    }
  }

  CFLogDebugMin( "ParCFmeshFileReader::readStoreInterNodes() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbElements(ifstream& fin)
{
  CFLog(NOTICE,"Memory Usage before assembling connectivity: " << Common::OSystem::getInstance().getProcessInfo()->memoryUsage() << "\n");

  CFLogDebugMin( "ParCFmeshFileReader::readNbElements() start\n");

  fin >> m_totNbElem;

  getReadData().setNbElements(m_totNbElem);

  if(m_totNbElem < 1) {
    throw BadFormatException (FromHere(),"Number of elements < 1 in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshFileReader::readNbElements() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbElementTypes(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbElementTypes() start\n");

  fin >> m_totNbElemTypes;
  getReadData().setNbElementTypes(m_totNbElemTypes);

  if (getReadData().getDimension() == DIM_1D && (m_totNbElemTypes < 1)) {
    throw BadFormatException (FromHere(),"Bad number of element types in 1D CFmesh");
  }

  if (getReadData().getDimension() == DIM_2D &&
      (m_totNbElemTypes < 1 || m_totNbElemTypes > 2)) {
    throw BadFormatException (FromHere(),"Bad number of element types in 2D CFmesh");
  }

  if(getReadData().getDimension() == DIM_3D &&
     (m_totNbElemTypes < 1 || m_totNbElemTypes > 4)) {
    throw BadFormatException (FromHere(),"Bad number of element types in 3D CFmesh");
  }

  CFLogDebugMin( "ParCFmeshFileReader::readNbElementTypes() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readGeometricPolyOrder(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readGeometricPolyOrder() start\n");

  CFuint order = 0;
  fin >> order;
  getReadData().setGeometricPolyOrder(static_cast<CFPolyOrder::Type>(order));

  if(getReadData().getGeometricPolyOrder() >= CFPolyOrder::MAXORDER ||
     getReadData().getGeometricPolyOrder() < CFPolyOrder::ORDER0) {
    throw BadFormatException (FromHere(),"Bad polynomial order of geometry in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshFileReader::readGeometricPolyOrder() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readSolutionPolyOrder(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readSolutionPolyOrder() start\n");

  CFuint order = 0;
  fin >> order;
  getReadData().setSolutionPolyOrder(static_cast<CFPolyOrder::Type>(order));

  if(getReadData().getSolutionPolyOrder() >= CFPolyOrder::MAXORDER ||
     getReadData().getSolutionPolyOrder() < CFPolyOrder::ORDER0) {
    throw BadFormatException (FromHere(),"Bad polynomial order of solution in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshFileReader::readSolutionPolyOrder() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readGeometricPolyType(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readGeometricPolyType() start\n");

  CFuint Type = 0;
  fin >> Type;
  getReadData().setGeometricPolyType(static_cast<CFPolyForm::Type>(Type));

  /// @todo should add check here for validity
  CFLogDebugMin( "ParCFmeshFileReader::readGeometricPolyType() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readSolutionPolyType(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readSolutionPolyType() start\n");

  CFuint Type = 0;
  fin >> Type;
  getReadData().setSolutionPolyType(static_cast<CFPolyForm::Type>(Type));

  /// @todo should add check here for validity
  CFLogDebugMin( "ParCFmeshFileReader::readSolutionPolyType() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readElementTypes(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readElementTypes() start\n");
  const CFuint nbElementTypes = getReadData().getNbElementTypes();
  cf_assert(nbElementTypes > 0);

  SafePtr< vector<ElementTypeData> > elementType =
    getReadData().getElementTypeData();
  elementType->resize(nbElementTypes);

  for (CFuint i = 0; i < nbElementTypes; ++i) {
    std::string elementShape = "";
    fin >> elementShape;

    (*elementType)[i].setShape(elementShape);

    (*elementType)[i].setGeoShape( CFGeoShape::Convert::to_enum(elementShape) );
  }

  CFLogDebugMin( "ParCFmeshFileReader::readElementTypes() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbElementsPerType(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readElementTypes() start\n");

  const CFuint nbElementTypes = getReadData().getNbElementTypes();
  cf_assert(nbElementTypes > 0);

  SafePtr< vector<ElementTypeData> > elementType =
    getReadData().getElementTypeData();

  CFuint sumNbElems = 0;
  CFint nbElemPerType = 0;

  for (CFuint i = 0; i < nbElementTypes; ++i) {
    fin >> nbElemPerType;
    (*elementType)[i].setNbElems(static_cast<CFuint>(nbElemPerType));

    sumNbElems += nbElemPerType;

    if (nbElemPerType < 0 || static_cast<CFuint>(nbElemPerType) >
  getReadData().getNbElements()) {
      throw BadFormatException (FromHere(),"Bad number of elements per type in CFmesh");
    }
  }

  if (sumNbElems != getReadData().getNbElements()){
    throw BadFormatException
      (FromHere(), "In CFmesh sum of elements per type differs from number elements");
  }

  CFLogDebugMin( "ParCFmeshFileReader::readElementTypes() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbNodesPerType(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbNodesPerType() start\n");

  const CFuint nbElementTypes = getReadData().getNbElementTypes();
  cf_assert(nbElementTypes > 0);

  Common::SafePtr< vector<ElementTypeData> > elementType =
    getReadData().getElementTypeData();

  CFint nbNodesPerType = 0;
  for (CFuint i = 0; i < nbElementTypes; ++i) {
    fin >> nbNodesPerType;

    (*elementType)[i].setNbNodes(static_cast<CFuint>(nbNodesPerType));

    CFLogDebugMin( "nbNodesPerType = " << nbNodesPerType << "\n");

    if(nbNodesPerType < 0) {
      throw BadFormatException
  (FromHere(), "In CFmesh, negative number of nodes per type");
    }

    /// @todo try to avoid having these hardcoded
    ///       only working for P1 geometry elements
    if(getReadData().getDimension() == DIM_1D &&
       nbNodesPerType != 2) {
      throw BadFormatException
	(FromHere(), "In CFmesh, bad number of nodes per 1D type");
    }

    /// @todo try to avoid having these hardcoded
    ///       only working for P1 and P2 geometry elements
    if(getReadData().getDimension() == DIM_2D &&
       nbNodesPerType != 3 &&
       nbNodesPerType != 4 &&
       nbNodesPerType != 6 &&
       nbNodesPerType != 8 &&
       nbNodesPerType != 9 &&
       nbNodesPerType != 10) {
      throw BadFormatException
  (FromHere(), "In CFmesh, bad number of nodes per 2D type");
    }

    if(getReadData().getDimension() == DIM_3D &&
       nbNodesPerType != 4 &&
       nbNodesPerType != 5 &&
       nbNodesPerType != 6 &&
       nbNodesPerType != 8 &&
       nbNodesPerType != 10 &&
       nbNodesPerType != 13 &&
       nbNodesPerType != 14 &&
       nbNodesPerType != 15 &&
       nbNodesPerType != 18 &&
       nbNodesPerType != 20 &&
       nbNodesPerType != 27) {
      throw BadFormatException
  (FromHere(), "In CFmesh, bad number of nodes per 3D type");
    }
  }
  CFLogDebugMin( "ParCFmeshFileReader::readNbNodesPerType() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbStatesPerType(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbStatesPerType() start\n");

  const CFuint nbElementTypes = getReadData().getNbElementTypes();
  cf_assert(nbElementTypes > 0);

  Common::SafePtr< vector<ElementTypeData> > elementType =
    getReadData().getElementTypeData();

  CFint nbStatesPerType = 0;
  for (CFuint i = 0; i < nbElementTypes; ++i) {
    fin >> nbStatesPerType;

    (*elementType)[i].setNbStates
      (static_cast<CFuint>(nbStatesPerType));

    CFLogDebugMin( "nbStatesPerType = " << nbStatesPerType << "\n");

    if(nbStatesPerType < 0) {
      throw BadFormatException
  (FromHere(), "In CFmesh, negative number of states per type");
    }
  }

  CFLogDebugMin( "ParCFmeshFileReader::readNbStatesPerType() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNodeList(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNodeList() start\n");

  const CFuint nbLocalNodes = m_localNodeIDs.size() + m_ghostNodeIDs.size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  const std::string parNodeVecName = nsp + "_nodes";
  
  DataHandle<Node*,GLOBAL> nodes = MeshDataStack::getActive()->getDataStorage()->
    getGlobalData<Node*>(parNodeVecName);
  
  cf_assert(nbLocalNodes > 0);
  nodes.reserve(nbLocalNodes, dim*sizeof(CFreal), nsp);
  getReadData().resizeNodes(nbLocalNodes);
  
  sort(m_localNodeIDs.begin(), m_localNodeIDs.end());
  sort(m_ghostNodeIDs.begin(), m_ghostNodeIDs.end());
  
  RealVector tmpNode(0.0, dim);
  RealVector tmpPastNode(0.0, dim);
  RealVector tmpInterNode(0.0, dim);

  if(!m_hasPastNodes && getReadData().storePastNodes()){
    throw BadFormatException
      (FromHere(), "ParCFmeshFileReader => readPastNodes is asked but PastNodes are not present in the CFmesh");
  }

  if(m_hasInterNodes && !getReadData().storeInterNodes()){
    CFout << "ParCFmeshFileReader => InterNodes are present in the CFmesh but will not be read\n";
  }

  if(!m_hasInterNodes && getReadData().storeInterNodes()){
    throw BadFormatException
      (FromHere(), "ParCFmeshFileReader => readInterNodes is asked but InterNodes are not present in the CFmesh");
  }

  if(m_hasPastNodes && !getReadData().storePastNodes()){
    CFout << "ParCFmeshFileReader => PastNodes are present in the CFmesh but will not be read\n";
  }

  const CFuint nbExtraVars = getReadData().getNbExtraNodalVars();
  const vector<CFuint>& nodalExtraVarsStrides = *getReadData().getExtraNodalVarStrides();

  RealVector extraVars;
  if (nbExtraVars > 0) {
    const CFuint sizeExtraVars = std::accumulate(nodalExtraVarsStrides.begin(),
              nodalExtraVarsStrides.end(), 0);
    extraVars.resize(sizeExtraVars);
    extraVars = 0.0;
  }

  getReadData().prepareNodalExtraVars();

  CFuint countLocals = 0;
  for (CFuint iNode = 0; iNode < m_totNbNodes; ++iNode) {

    // read the node
    fin >> tmpNode;

    if (m_hasPastNodes) {
      fin >> tmpPastNode;
    }

    if (m_hasInterNodes) {
      fin >> tmpInterNode;
    }

    if (nbExtraVars > 0) {
      fin >> extraVars;
    }

    CFuint localID = 0;
    bool isGhost = false;
    bool isFound = false;
    if (hasEntry(m_localNodeIDs, iNode)) {
      countLocals++;
      localID = nodes.addLocalPoint (iNode);
      cf_assert(localID < nbLocalNodes);
      isFound = true;
    }
    else if (hasEntry(m_ghostNodeIDs, iNode)) {
      countLocals++;
      localID = nodes.addGhostPoint (iNode);
      cf_assert(localID < nbLocalNodes);
      isGhost = true;
      isFound = true;
    }

    if (isFound) {
      Node* newNode = getReadData().createNode
  (localID, nodes.getGlobalData(localID), tmpNode, !isGhost);
      newNode->setGlobalID(iNode);
      if (m_hasPastNodes) {
        getReadData().setPastNode(localID, tmpPastNode);
      }
     if (m_hasInterNodes) {
        getReadData().setInterNode(localID, tmpInterNode);
      }


      // set the nodal extra variable
      if (nbExtraVars > 0) {
  getReadData().setNodalExtraVar(localID, extraVars);
      }
    }
  }

  cf_assert(countLocals == nbLocalNodes);

  CFLogDebugMin("countLocals  = " << countLocals << "\n");
  CFLogDebugMin("nbLocalNodes = " << nbLocalNodes << "\n");
  CFLogDebugMin("m_localNodeIDs.size() = " << m_localNodeIDs.size() << "\n");
  CFLogDebugMin("m_ghostNodeIDs.size() = " << m_ghostNodeIDs.size() << "\n");

  CFLogDebugMin( "ParCFmeshFileReader::readNodeList() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::emptyNodeListRead(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::emptyNodeListRead() start" << "\n");

  cf_assert(m_totNbNodes > 0);

  RealVector node(0.0, getReadData().getDimension());
  RealVector tmpPastNode(0.0, getReadData().getDimension());
  RealVector tmpInterNode(0.0, getReadData().getDimension());

  if(!m_hasPastNodes && getReadData().storePastNodes()){
    throw BadFormatException
      (FromHere(), "ParCFmeshFileReader => readPastNodes is asked but PastNodes are not present in the CFmesh");
  }

  if(m_hasPastNodes && !getReadData().storePastNodes()){
    CFout << "ParCFmeshFileReader => PastNodes are present in the CFmesh but will not be read\n";
  }

  if(!m_hasInterNodes && getReadData().storeInterNodes()){
    throw BadFormatException
      (FromHere(), "ParCFmeshFileReader => readInterNodes is asked but InterNodes are not present in the CFmesh");
  }

  if(m_hasInterNodes && !getReadData().storeInterNodes()){
    CFout << "ParCFmeshFileReader => InterNodes are present in the CFmesh but will not be read\n";
  }

  const CFuint nbExtraVars = getReadData().getNbExtraNodalVars();
  const vector<CFuint>& nodalExtraVarsStrides = *getReadData().getExtraNodalVarStrides();
  RealVector extraVars;
  if (nbExtraVars > 0) {
    const CFuint sizeExtraVars = std::accumulate(nodalExtraVarsStrides.begin(),
              nodalExtraVarsStrides.end(), 0);
    extraVars.resize(sizeExtraVars);
    extraVars = 0.0;
  }

  getReadData().prepareNodalExtraVars();

  for (CFuint n = 0; n < m_totNbNodes; ++n) {
    fin >> node;

    if (m_hasPastNodes) {
      fin >> tmpPastNode;
    }

    if (m_hasInterNodes) {
      fin >> tmpInterNode;
    }


    if (nbExtraVars > 0) {
      fin >> extraVars;
    }
  }

  CFLogDebugMin( "ParCFmeshFileReader::emptyNodeListRead() end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readStateList(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readStateList() start\n");
  
  if ((m_initValues.size()    != m_useInitValues.size()) && (m_initValuesIDs.size() != m_useInitValues.size())) {
    CFLog(VERBOSE, "ParCFmeshFileReader::readStateList() => m_initValues.size()    = " << m_initValues.size() << "\n");
    CFLog(VERBOSE, "ParCFmeshFileReader::readStateList() => m_initValuesIDs.size() = " << m_initValuesIDs.size() << "\n");
    CFLog(VERBOSE, "ParCFmeshFileReader::readStateList() => m_useInitValues.size() = " << m_useInitValues.size() << "\n");
    
    throw BadFormatException
      (FromHere(),"ParCFmeshFileReader => m_initValues && m_initValuesIDs sizes != m_useInitValues.size()");
  }
  
  bool isWithSolution = false;
  fin >> isWithSolution;
  getReadData().setWithSolution(isWithSolution);

  const CFuint nbLocalStates = m_localStateIDs.size() + m_ghostStateIDs.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  const std::string parStateVecName = nsp + "_states";
  
  DataHandle<State*,GLOBAL> states = MeshDataStack::getActive()->getDataStorage()->
    getGlobalData<State*>(parStateVecName);
  
  cf_assert(nbLocalStates > 0);
  states.reserve(nbLocalStates, nbEqs*sizeof(CFreal), nsp);
  getReadData().resizeStates(nbLocalStates);
  
  sort(m_localStateIDs.begin(), m_localStateIDs.end());
  sort(m_ghostStateIDs.begin(), m_ghostStateIDs.end());

  State tmpState;
  State dummyReadState;
  RealVector tmpPastState(0.0, nbEqs);
  RealVector tmpInterState(0.0, nbEqs);
  RealVector readState(0.0, m_originalNbEqs);
  cf_assert(m_originalNbEqs > 0);

  if(!m_hasPastStates && getReadData().storePastStates()){
    throw BadFormatException
      (FromHere(),"ParCFmeshFileReader => readPastStates is asked but PastStates are not present in the CFmesh");
  }

  if(m_hasPastStates && !getReadData().storePastStates()){
    CFout  << "ParCFmeshFileReader => PastStates are present in the CFmesh but will not be read\n";
  }

 if(!m_hasInterStates && getReadData().storeInterStates()){
    throw BadFormatException
      (FromHere(),"ParCFmeshFileReader => readInterStates is asked but InterStates are not present in the CFmesh");
  }

  if(m_hasInterStates && !getReadData().storeInterStates()){
    CFout  << "ParCFmeshFileReader => InterStates are present in the CFmesh but will not be read\n";
  }

  const CFuint nbExtraVars = getReadData().getNbExtraStateVars();
  const vector<CFuint>& stateExtraVarsStrides = *getReadData().getExtraStateVarStrides();
  RealVector extraVars;
  if (nbExtraVars > 0) {
    const CFuint sizeExtraVars = std::accumulate(stateExtraVarsStrides.begin(),
              stateExtraVarsStrides.end(), 0);
    extraVars.resize(sizeExtraVars);
    extraVars = 0.0;
  }

  getReadData().prepareStateExtraVars();

  bool hasTransformer = false;

  // read the state
  if (isWithSolution) { // warn if nbeqs differs from original and we dont provide mapping of variable ids
    if (m_originalNbEqs != nbEqs && (m_useInitValues.size() == 0))
      {
        CFLog(WARN, "************************************************************\n");
        CFLog(WARN, "WARNING: Nb equations in CFmesh differs from Physical Model.\n");
        CFLog(WARN, "         User did not provide mapping IDs.\n");
        CFLog(WARN, "         May incur in wrong initialization of solution.\n");
        CFLog(WARN, "************************************************************\n");
      }

    // configure the variable transformer
    const std::string name = MeshDataStack::getActive()->getPrimaryNamespace();
    SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
      (SubSystemStatusStack::getCurrentName()).getNamespace(name);
    SafePtr<PhysicalModel> physModel =
      PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

    SafePtr<VarSetTransformer::PROVIDER> vecTransProv = CFNULL;

    CFLog(VERBOSE, "ParCFmeshFileReader::configure(): configuring " << m_inputToUpdateVecStr << "\n");
    try
    {
      vecTransProv = Environment::Factory<VarSetTransformer>::getInstance().getProvider(m_inputToUpdateVecStr);
    }
    catch (NoSuchValueException& except) {
      m_inputToUpdateVecStr = "Identity";

      CFLog(VERBOSE, except.what() << "\n");
      CFLog(VERBOSE, "Choosing IdentityVarSetTransformer instead ..." << "\n");
      vecTransProv = Environment::Factory<VarSetTransformer>::getInstance().getProvider(m_inputToUpdateVecStr);
    }

    cf_assert(vecTransProv.isNotNull());
    m_inputToUpdateVecTrans.reset(vecTransProv->create(physModel->getImplementor()));
    cf_assert(m_inputToUpdateVecTrans.getPtr() != CFNULL);

    if (m_inputToUpdateVecStr != "Identity") {
      hasTransformer = true;
    }

    m_inputToUpdateVecTrans->setup(1);
  }
  
  CFuint countLocals = 0;
  for (CFuint iState = 0; iState < m_totNbStates; ++iState)
  {
    // read the state
    if (isWithSolution) 
    {      
      // no init values were used
      if (m_useInitValues.size() == 0)
      {
	fin >> readState;

        if (m_hasPastStates) 
        {
          fin >> tmpPastState;
        }
	
	if (m_hasInterStates) {
          fin >> tmpInterState;
        }

        if (nbExtraVars > 0) {
          fin >> extraVars;
        }

        if (!hasTransformer) {
          const CFuint currNbEqs = std::min(nbEqs,m_originalNbEqs); // AL: why min????
          for (CFuint iEq = 0; iEq < currNbEqs; ++iEq) {
            tmpState[iEq] = readState[iEq];
          }
        }
        else {
          for (CFuint iEq = 0; iEq < m_originalNbEqs; ++iEq) {
            dummyReadState[iEq] = readState[iEq];
          }
          tmpState = *m_inputToUpdateVecTrans->transform(&dummyReadState);
        }
      }
      
      // using init values
      else {
	cf_assert(m_useInitValues.size() == nbEqs);
	fin >> readState;
	
	if (m_hasPastStates) {
	  fin >> tmpPastState;
	}
	
	if (m_hasInterStates) {
	  fin >> tmpInterState;
	}
	
	if (nbExtraVars > 0) {
	  fin >> extraVars;
	}
	
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  if (!m_useInitValues[iEq] && iEq < m_originalNbEqs) {
	    cf_assert(iEq < tmpState.size());
	    tmpState[iEq] = readState[iEq];
	  }
          else {
	    // user must specify either all initial values or values IDs, NOT BOTH
	    cf_assert(m_initValues.size() != m_initValuesIDs.size());
	    
	    if (m_initValues.size() > 0) {
              cf_assert(m_initValuesIDs.size() == 0);
	      cf_assert(iEq < tmpState.size());
	      cf_assert(iEq < m_initValues.size());
              tmpState[iEq] = m_initValues[iEq];
            }
	    
            if (m_initValuesIDs.size() > 0) {
	      cf_assert(m_initValues.size() == 0);
	      const CFuint currID = m_initValuesIDs[iEq];
	      // if the current ID is >= nbEqs set this variable to 0.0
	      tmpState[iEq] = (currID < m_originalNbEqs) ? readState[currID] : 0.0;
	    }
          }
        }
	
        // in case the original nb of equations in the file
        // is bigger than the current number of equations
        // we read the rest of the states and discard them
        if (m_originalNbEqs > nbEqs)
        {
	  for (CFuint iEq = nbEqs; iEq < m_originalNbEqs; ++iEq)
	  {
            fin >> readState[iEq];
          }
        }
      }
    }

    CFuint localID = 0;
    bool isGhost = false;
    bool isFound = false;
    if (hasEntry(m_localStateIDs, iState)) {
      countLocals++;
      localID = states.addLocalPoint (iState);
      cf_assert(localID < nbLocalStates);
      isFound = true;
    }
    else if (hasEntry(m_ghostStateIDs, iState)) {
      countLocals++;
      localID = states.addGhostPoint (iState);
      cf_assert(localID < nbLocalStates);
      isGhost = true;
      isFound = true;
    }

    if (isFound) {
      State* newState = getReadData().createState
  (localID, states.getGlobalData(localID), tmpState, !isGhost);
      newState->setGlobalID(iState);

      if (m_hasPastStates) {
        getReadData().setPastState(localID, tmpPastState);
      }

      if (m_hasInterStates) {
        getReadData().setInterState(localID, tmpInterState);
      }
      // set the nodal extra variable
      if (nbExtraVars > 0) {
  getReadData().setStateExtraVar(localID, extraVars);
      }

      // getReadData().setStateLocalToGlobal (localID, globalID);
      // getReadData().setLocalState(localID, !isGhost);
    }
  }

  cf_assert(countLocals == nbLocalStates);

  CFLogDebugMin( "ParCFmeshFileReader::readStateList() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::emptyStateListRead(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::emptyStateListRead() start\n");

  bool isWithSolution = false;
  fin >> isWithSolution;

  if ((m_initValues.size() != m_useInitValues.size()) &&
      (m_initValuesIDs.size() != m_useInitValues.size())) {
    throw BadFormatException
      (FromHere(), "ParCFmeshFileReader => m_initValues && m_initValuesIDs sizes != m_useInitValues.size()");
  }

  cf_assert(m_totNbStates > 0);
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  RealVector readState(0.0, m_originalNbEqs);
  RealVector tmpState(0.0, nbEqs);
  cf_assert(m_originalNbEqs > 0);

  const CFuint nbExtraVars = getReadData().getNbExtraStateVars();
  const vector<CFuint>& stateExtraVarsStrides = *getReadData().getExtraStateVarStrides();
  RealVector extraVars;
  if (nbExtraVars > 0) {
    const CFuint sizeExtraVars = std::accumulate(stateExtraVarsStrides.begin(),
              stateExtraVarsStrides.end(), 0);
    extraVars.resize(sizeExtraVars);
    extraVars = 0.0;
  }

  getReadData().prepareStateExtraVars();

  if (isWithSolution) {
    for (CFuint s = 0; s < m_totNbStates; ++s) {
      // read the state values if they exist
      if (m_useInitValues.size() == 0)
      {
  fin >> readState;
        if (nbExtraVars > 0) {
          fin >> extraVars;
        }

  for (CFuint iEq = 0; iEq < std::min(nbEqs,m_originalNbEqs); ++iEq)
    tmpState[iEq] = readState[iEq];
      }

      // using init values
      else {

        cf_assert(m_useInitValues.size() == nbEqs);
  fin >> readState;
        if (nbExtraVars > 0) {
          fin >> extraVars;
        }


  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    if (!m_useInitValues[iEq] && iEq < m_originalNbEqs) {
      tmpState[iEq] = readState[iEq];
    }
    else {
      cf_assert(m_initValues.size() != m_initValuesIDs.size());
      if (m_initValues.size() > 0) {
        cf_assert(m_initValuesIDs.size() == 0);
        tmpState[iEq] = m_initValues[iEq];
      }

      if (m_initValuesIDs.size() > 0) {
        const CFuint currID = m_initValuesIDs[iEq];
        // if the current ID is >= nbEqs set this variable to 0.0
        tmpState[iEq] = (currID < m_originalNbEqs) ? readState[currID] : 0.0;
      }
    }
  }

  // in case the original nb of equations in the file
        // is bigger than the current number of equations
        // we read the rest of the states and discard them
        if (m_originalNbEqs > nbEqs)
  {
    for (CFuint iEq = nbEqs; iEq < m_originalNbEqs; ++iEq)
    {
      fin >> readState[iEq];
    }
        }
      }
    }
  }

  CFLogDebugMin( "ParCFmeshFileReader::emptyStateListRead() end" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readElementList(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readElementList() start\n");

  // allocate the partitioner data
  PartitionerData pdata;

  // set the local coloring array data inside the partitioner data
  // for later usage
  pdata.part = &m_partitionerOutData;

  // set the element distribution array
  setElmDistArray(pdata.elmdist);
  setSizeElemVec(pdata.sizeElemNodeVec, pdata.sizeElemStateVec);

  const CFuint sizeNodeElem = pdata.sizeElemNodeVec[m_myRank];
  const CFuint sizeStateElem = pdata.sizeElemStateVec[m_myRank];
  pdata.elemNode.resize(sizeNodeElem);
  pdata.elemState.resize(sizeStateElem);
  pdata.eptrn.resize(m_nbElemPerProc[m_myRank] + 1);
  pdata.eptrs.resize(m_nbElemPerProc[m_myRank] + 1);
  
  readElemListRank(pdata, fin);
  pdata.ndim=(CFint)PhysicalModelStack::getActive()->getDim();
  
  // do the partitioning of the mesh
  // global element IDs local to each processor after the partitioning
  // will be placed in pdata.part
  cf_assert(m_partitioner.isNotNull());
  
  CFLog(DEBUG_MAX, CFPrintContainer<vector<PartitionerData::IndexT> >("pdata.eptrn    = ", &pdata.eptrn));
  CFLog(DEBUG_MAX, CFPrintContainer<vector<PartitionerData::IndexT> >("pdata.elemNode = ", &pdata.elemNode));
  CFLog(DEBUG_MAX, CFPrintContainer<vector<PartitionerData::IndexT> >("pdata.elmdist  = ", &pdata.elmdist));
  
  // avoid mesh partitioning if you have just one processor
  if (m_nbProc > 1)
  {
    if (PhysicalModelStack::getActive()->getDim() > DIM_1D) {
      m_partitioner->SetCommunicator(m_comm);
      CFLog(NOTICE, "Calling mesh partitioner\n");
      CFLog(NOTICE, "+++\n");
      m_partitioner->doPartition(pdata);
      CFLog(NOTICE, "+++\n");
    }
    else {
      pdata.part->resize(m_nbElemPerProc[m_myRank], m_myRank);
    }
  }
  else
  {
    cf_assert(m_myRank == 0);
    pdata.part->resize(m_nbElemPerProc[m_myRank], 0);
  }
  
  // move the elements data to the right processor
  // and build info about the overlap region
  m_local_elem = new ElementDataArray<0>;
  moveElementData(*m_local_elem, pdata);

  cf_assert(m_localNodeIDs.size() > 0);
  cf_assert(m_localStateIDs.size() > 0);

  // set mapping between global and local node/state IDs
  setMapGlobalToLocalID(m_localNodeIDs, m_ghostNodeIDs, m_mapGlobToLocNodeID);

  setMapGlobalToLocalID(m_localStateIDs, m_ghostStateIDs, m_mapGlobToLocStateID);

  // set the elements in the readData
  setElements(*m_local_elem);

  setMapNodeElemID(*m_local_elem);

  CFLogDebugMin( "ParCFmeshFileReader::readElementList() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readElemListRank(PartitionerData& pdata,
					   ifstream& fin)
{
  CFuint start = 0;
  for (CFuint rank = 0; rank < m_myRank; ++rank) {
    start += m_nbElemPerProc[rank];
  }
  const CFuint ne = m_nbElemPerProc[m_myRank];
  const CFuint end = start + ne;

  SafePtr< vector<ElementTypeData> > elementType =
    getReadData().getElementTypeData();
  
  vector<PartitionerData::IndexT>& eNode  = pdata.elemNode;
  vector<PartitionerData::IndexT>& eState = pdata.elemState;
  vector<PartitionerData::IndexT>& eptrn  = pdata.eptrn;
  vector<PartitionerData::IndexT>& eptrs  = pdata.eptrs;
  
  const CFuint nbEEminProc       = m_nbElemPerProc[m_myRank];
  const CFuint nbEEminProcMinus1 = nbEEminProc - 1;
  CFuint ncount = 0;
  CFuint scount = 0;
  CFuint ipos = 0;
  CFuint iElemBegin = 0;
  CFuint nodeID = 0;
  CFuint stateID = 0;
  
  for (CFuint iType = 0; iType < m_totNbElemTypes; ++iType) {
    const CFuint nbNodesInElem  = (*elementType)[iType].getNbNodes();
    const CFuint nbStatesInElem = (*elementType)[iType].getNbStates();
    const CFuint nbElementsPerType = (*elementType)[iType].getNbElems();
    const CFuint iElemEnd = iElemBegin + nbElementsPerType;
    
    // loop over the elements in this type
    for (CFuint iElem = iElemBegin; iElem < iElemEnd; ++iElem) {
      if (iElem < start || iElem >= end) {
	for (CFuint iNode = 0; iNode < nbNodesInElem; ++iNode) {
	  fin >> nodeID;
	  checkDofID("node", iElem, iNode, nodeID, m_totNbNodes);
	}
	for (CFuint iState = 0; iState < nbStatesInElem; ++iState) {
	  fin >> stateID;
	  checkDofID("state", iElem, iState, stateID, m_totNbStates);
	}
      }
      else {
	// if the element belongs to the current rank store it
	cf_assert(iElem >= start || iElem < end);
	
	eptrn[ipos] = ncount;
	eptrs[ipos] = scount;
	
	for (CFuint j = 0; j < nbNodesInElem; ++j, ++ncount) {
	  fin >> eNode[ncount];
	  checkDofID("node", iElem, j, eNode[ncount], m_totNbNodes);
	}
	for (CFuint j = 0; j < nbStatesInElem; ++j, ++scount) {
	  fin >> eState[scount];
	  checkDofID("state", iElem, j, eState[scount], m_totNbStates);
	}
	
	if (ipos == nbEEminProcMinus1) {
	  eptrn[nbEEminProc] = eptrn[ipos] + nbNodesInElem;
	  eptrs[nbEEminProc] = eptrs[ipos] + nbStatesInElem;
	}
	++ipos;
      }
    }
    
    iElemBegin +=  nbElementsPerType;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbTRSs(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbTRSs() start\n");

  CFuint nbTRSs = 0;
  fin >> nbTRSs;

  if (nbTRSs < 1)
  {
    throw BadFormatException (FromHere(),"Number of TRSs in file must be at least 1");
  }

  // account for the reduction of TRS's
  nbTRSs -= m_trs_reduction;

  if (nbTRSs < 1)
  {
    throw BadFormatException (FromHere(),"Number of TRSs after merging must be at least 1");
  }

  getReadData().setNbTRSs(nbTRSs);
  getReadData().resizeGeoConn(getReadData().getNbTRSs());
  getReadData().getNbGeomEntsPerTR()->resize(nbTRSs);
  getReadData().getNbTRs()->resize(nbTRSs);
  // Set some global info
  SwapEmpty(MeshDataStack::getActive()->getTotalTRSInfo());
  SwapEmpty(MeshDataStack::getActive()->getTotalTRSNames ());
  SwapEmpty(*MeshDataStack::getActive()->getGlobalTRSGeoIDs());

  MeshDataStack::getActive()->getTotalTRSInfo().resize(nbTRSs);
  MeshDataStack::getActive()->getTotalTRSNames().resize(nbTRSs);
  // AL: here resize() is used on purpose to allow direct subscripting afterwards
  MeshDataStack::getActive()->getGlobalTRSGeoIDs()->resize(nbTRSs);

  CFLogDebugMin( "ParCFmeshFileReader::readNbTRSs() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readTRSName(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readTRSName() start\n");

  std::string name = "";
  fin >> name;

  CFLogDebugMin( "Found TRS " + name + "\n");

  // check if TRS is to be merged
  m_mergedtrs = false;
  m_mergedtrs_just_added = false;
  if (m_rev_trs_merge.count(name))
  {
    // should be merged
    std::string sumtrs = m_rev_trs_merge[name];
    CFLogDebugMin( "TRS " << name << " should be merged with " << sumtrs << "\n");
    CFLogDebugMin( "TRS " << sumtrs << " was merged ntimes : " << m_trs_idxmap.count(sumtrs) << "\n");
    // if not added already, add the summed TRS
    if (!m_trs_idxmap.count(sumtrs))
    {
      CFLogDebugMin( "Creating merged TRS " << sumtrs << "\n");

      const CFuint idx = m_trs_idxmap.size();
      m_trs_idxmap[sumtrs] = idx;
      getReadData().getNameTRS()->push_back(sumtrs);
      MeshDataStack::getActive()->getTotalTRSNames()[idx] = sumtrs;
      m_mergedtrs_just_added = true;
    }
    // set the name to the merged TRS
    m_mergedtrs = true;
    name = sumtrs;
  }
  else // not to be merged
  {
    CFLogDebugMin( "Creating TRS " << name << "\n");

    const CFuint idx = m_trs_idxmap.size();
    m_trs_idxmap[name] = idx;
    getReadData().getNameTRS()->push_back(name);
    MeshDataStack::getActive()->getTotalTRSNames()[idx] = name;
  }

  // set the current TRS being read
  m_curr_trs = name;
  // index
  CFLogDebugMin( "Current TRS " << m_curr_trs << " has index " << m_trs_idxmap[m_curr_trs] << "\n");

  CFLogDebugMin( "ParCFmeshFileReader::readTRSName() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbTRs(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbTRs() start\n");

  CFuint nbTRsInTRS = 0;
  fin >> nbTRsInTRS;

  cf_assert(nbTRsInTRS > 0);

  const CFuint idx = m_trs_idxmap[m_curr_trs];
  (*getReadData().getNbTRs())[idx] += nbTRsInTRS;

  CFLogDebugMin( "TRS " << m_curr_trs << " + " << nbTRsInTRS << "TR, total " << (*getReadData().getNbTRs())[idx] << " TR\n");

  // set the current number of TR's being read
  m_curr_nbtr = nbTRsInTRS;

  CFLogDebugMin( "ParCFmeshFileReader::readNbTRs() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbGeomEnts(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbGeomEnts() start\n");

  const CFuint idx = m_trs_idxmap[m_curr_trs];
  for (CFuint iTr = 0; iTr < m_curr_nbtr; ++iTr)
  {
    CFint nbgeo = 0;
    fin >> nbgeo;
    cf_assert(nbgeo > 0);

    (*getReadData().getNbGeomEntsPerTR())[idx].push_back(nbgeo);
    MeshDataStack::getActive()->getTotalTRSInfo()[idx].push_back(nbgeo);
  }

  CFLogDebugMin( "ParCFmeshFileReader::readNbGeomEnts() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readGeomType(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readGeomType() start\n");

  std::string geomName = "";
  fin >> geomName;

  CFGeoEnt::Type trs_type = CFGeoEnt::Convert::to_enum ( geomName );

  // check if TRS is to be merged
  if (m_mergedtrs && !m_mergedtrs_just_added)
  {
    const CFuint idx = m_trs_idxmap[m_curr_trs];
    CFGeoEnt::Type sum_trs_type = (*getReadData().getGeomType())[idx];
    std::string sum_trs_type_name =
      CFGeoEnt::Convert::to_str(sum_trs_type);
    if ( sum_trs_type != trs_type )
      throw BadValueException (FromHere(), "Trying to merge TRS's of different types "
            + geomName + " and " + sum_trs_type_name);
  }
  else
  {
    getReadData().getGeomType()->push_back(trs_type);
  }

  CFLogDebugMin( "ParCFmeshFileReader::readGeomType() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readGeomEntList(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readGeomEntList() start\n");

  typedef CFMultiMap<CFuint,CFuint>::MapIterator MapItr;

  // load TRs data into memory for further use
  // this is actually only useful to be able to write file
  // without having constructed TRSs
  // get id of last TRS read from file

  SafePtr< vector<CFuint> > nbTRs = getReadData().getNbTRs();
  SafePtr< vector<vector<CFuint> > > nbGeomEntsPerTR =
    getReadData().getNbGeomEntsPerTR();

  const CFuint iTRS = m_trs_idxmap[m_curr_trs];
  const CFuint nbTRsAdded = (*nbTRs)[iTRS];

  getReadData().resizeGeoConn(iTRS, nbTRsAdded);

  CFLogDebugMin("Rank " << m_myRank << " CurrTRS = " << m_curr_trs << "\n");
  CFLogDebugMin("Rank " << m_myRank << " iTRS = " << iTRS << "\n");
  CFLogDebugMin("Rank " << m_myRank << " nbTRsAdded = " << nbTRsAdded << "\n");

  // Set some global info
  SafePtr<vector<vector<vector<CFuint> > > > trsGlobalIDs =
    MeshDataStack::getActive()->getGlobalTRSGeoIDs();
  (*trsGlobalIDs)[iTRS].resize(nbTRsAdded);

  pair<std::valarray<CFuint>, std::valarray<CFuint> > geoConLocal;

  // loop only in the new TRs, which have not been read yet
  for (CFuint iTR = nbTRsAdded - m_curr_nbtr; iTR < nbTRsAdded; ++iTR)
  {
    CFLogDebugMin("Rank " << m_myRank << " iTR = "<< iTR << "\n");
    const CFuint nbTRGeos = (*nbGeomEntsPerTR)[iTRS][iTR];
    CFLogDebugMin("Rank " << m_myRank << " nbTRGeos = "<< nbTRGeos << "\n");

    CFuint countGeos = 0;
    for (CFuint iGeo = 0; iGeo < nbTRGeos; ++iGeo)
    {
      CFuint nbNodesInGeo = 0;
      CFuint nbStatesInGeo = 0;
      fin >> nbNodesInGeo;
      fin >> nbStatesInGeo;

      geoConLocal.first.resize(nbNodesInGeo);
      geoConLocal.second.resize(nbStatesInGeo);

      for(CFuint n = 0; n < nbNodesInGeo; ++n)
      {
        fin >> geoConLocal.first[n];
        cf_assert(geoConLocal.first[n] < m_totNbNodes);
      }

      for(CFuint s = 0; s < nbStatesInGeo; ++s)
      {
        fin >> geoConLocal.second[s];
        cf_assert(geoConLocal.second[s] < m_totNbStates);
      }

      // check if the global ID of the first node of the
      // geometric entity is referenced by any local element
      bool nodeFound = false;
      pair<MapItr, MapItr> etr =
       m_mapNodeElemID.find(geoConLocal.first[0],nodeFound);

      // if the first one is found, check if all the other
      // GE nodes are referenced by one amongst all vertex-neighbor elements
      if (nodeFound)
      {
        bool exitLoop = false;
        for (MapItr etm = etr.first; (etm != etr.second) && (!exitLoop); ++etm) {
          const CFuint localElemID = etm->second;
          const CFuint nbENodes = getReadData().
          getNbNodesInElement(localElemID);
          CFuint counter = 1; // the first node already matches
          for (CFuint in = 1; in < nbNodesInGeo; ++in) {
            bool hasLocalID = false;
            const CFuint localNodeID = m_mapGlobToLocNodeID.
            find(geoConLocal.first[in], hasLocalID);
            // if the flag is false, this global node ID is not
            // referenced by local elements, then skip this GE
            if (!hasLocalID) {
              exitLoop = true;
              break;
            }

      // search the local node ID among the nodes of the current element
            for (CFuint jn = 0; jn < nbENodes; ++jn) {
              const CFuint nodeID = getReadData().
                getElementNode(localElemID, jn);
                if (nodeID == localNodeID)
                {
                  counter++;
                  break;
                }
            }
          }

          // if all the nodes of the current GE are included
    // in the given element, store the current GE as local
          if (counter == nbNodesInGeo) {
            // convert the global in local node/state IDs
            for(CFuint n = 0; n < nbNodesInGeo; ++n) {
              geoConLocal.first[n] = m_mapGlobToLocNodeID.
              find(geoConLocal.first[n]);
            }

            for(CFuint s = 0; s < nbStatesInGeo; ++s) {
              geoConLocal.second[s] = m_mapGlobToLocStateID.
              find(geoConLocal.second[s]);
            }

            getReadData().addGeoConn(iTRS, iTR, geoConLocal);

            // set the global ID of the geometric entity inside this TR and TRS
            (*trsGlobalIDs)[iTRS][iTR].push_back(iGeo);

            // increment the counter of the nb of GEs
            countGeos++;
            break;
          }
        }
      } // found node

    } // loop iGeo

    // reset the number of GEs in the current TR
    (*nbGeomEntsPerTR)[iTRS][iTR] = countGeos;

    cf_assert((*trsGlobalIDs)[iTRS][iTR].size() == countGeos);

    CFLogDebugMin("Rank " << m_myRank << ", iTR = " << iTR
      << ", countGeos = " << countGeos << "\n");
  }

  CFLogDebugMin( "ParCFmeshFileReader::readGeomEntList() end\n");
}

//////////////////////////////////////////////////////////////////////////////

bool ParCFmeshFileReader::readString(ifstream& file)
{
  std::string key = "";
  file >> key;

  CFLogDebugMin("CFmesh key = " << key << "\n");

  // check end of file
  if (key != getReaderTerminator()) {

    MapString2Reader::iterator key_pair = m_mapString2Reader.find(key);

    // check if key exists else ignore it
    if (key_pair != m_mapString2Reader.end()) {

      // if the elements are not built, ignore list of nodes and states
      // and keep on reading
      if(((key == "!LIST_NODE") && !areElementsBuild()) ||
   ((key == "!LIST_STATE") && !areElementsBuild())) {

  CFLog(WARN, "Warning: old CFmesh format file. Node and state lists will be built later.\n");

  if (key == "!LIST_NODE") {
    m_startNodeList = file.tellg();
    emptyNodeListRead(file);
  }

  if (key == "!LIST_STATE") {
    m_startStateList = file.tellg();
    emptyStateListRead(file);
  }

  // this can only occur during the first reading
  cf_assert(getReadCount() == 0);

  // at the end of the first reading
  setReadAgain(true);
  return true;
      }

      // if this is the second reading, only read list of nodes and states
      if (getReadCount() > 0) {
  // read node list
  file.seekg(m_startNodeList);
  readNodeList(file);

  // read state list
  file.seekg(m_startStateList);
  readStateList(file);
  setReadAgain(false);

  return false;
      }
      else {
  CFLogDebugMin( "Read CFmesh Key: " << key << "\n");
  ReaderFunction function = key_pair->second;
  cf_assert(function != CFNULL);
  (this->*function)(file);
      }
    }
    else {
      std::string msg = "Key in CFmesh is not valid:" + key;
      throw Common::NoSuchValueException (FromHere(),msg);
    }
    // keep reading
    return true;
  }

  // found end of file
  return false;
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::setElmDistArray(vector<PartitionerData::IndexT>& elmdist)
{
  elmdist.resize(m_nbProc + 1);
  const CFuint ne = m_totNbElem/m_nbProc;
  m_nbElemPerProc[0] = ne + m_totNbElem%m_nbProc;
  elmdist[0] = 0;
  CFuint elmd = m_nbElemPerProc[0];
  for (CFuint i = 1; i < m_nbProc; ++i) {
    m_nbElemPerProc[i] = ne;
    elmdist[i] = elmd;
    elmd += m_nbElemPerProc[i];
  }
  elmdist[m_nbProc] = elmd;

  cf_assert(std::accumulate(m_nbElemPerProc.begin(),
                    m_nbElemPerProc.end(),static_cast<CFuint>(0)) == m_totNbElem);
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::setSizeElemVec(vector<PartitionerData::IndexT>& sizeElemNodeVec,
					 vector<PartitionerData::IndexT>& sizeElemStateVec)
{
  sizeElemNodeVec.resize(m_nbProc, static_cast<CFuint>(0));
  sizeElemStateVec.resize(m_nbProc, static_cast<CFuint>(0));

  SafePtr< vector<ElementTypeData> > elementType =
    getReadData().getElementTypeData();

  CFuint start = 0;
  CFuint elemID = 0;
  for (CFuint ip = 0; ip < m_nbProc; ++ip) {
    const CFuint nbEP = m_nbElemPerProc[ip];
    const CFuint end = start + nbEP;
    for (CFuint ie = start; ie < end; ++ie, ++elemID) {
      const CFuint iType = getElementType(*elementType, elemID);
      sizeElemNodeVec[ip]  += (*elementType)[iType].getNbNodes();
      sizeElemStateVec[ip] += (*elementType)[iType].getNbStates();
    }
    start += nbEP;
  }
}

//////////////////////////////////////////////////////////////////////////////
      
void ParCFmeshFileReader::moveElementData(ElementDataArray<0>& localElem, 
					  PartitionerData& pdata )
{
  // set the global element IDs
  const CFuint partSize = pdata.part->size();
  vector<CFuint> globalElemIDs(partSize);
  
  // compute the starting element ID
  CFuint elemID = 0;
  for (CFuint i = 0; i <m_nbProc; ++i) {
    if (i == m_myRank) break;
    else {
      elemID += m_nbElemPerProc[i];
    }
  }
  
  const CFuint startElemID = elemID;
  CFLogDebugMin(m_myRank << " startElemID = " << startElemID << "\n");
  vector<CFuint> elementSize(partSize);

  for (CFuint i = 0; i < partSize; ++i, ++elemID) {
    globalElemIDs[i] = elemID;
    elementSize[i] = 3 + pdata.eptrn[i+1] + pdata.eptrs[i+1];
  }
  
  CFMultiMap<PartitionerData::IndexT, CFuint> mapProcToOldElemID(partSize);
  sortPartVec(globalElemIDs, *pdata.part, mapProcToOldElemID);
  
  // counts for data to send
  vector<int> sendCount(m_nbProc, static_cast<CFuint>(0));
  // counts for data to receive
  vector<int> recvCount(m_nbProc, static_cast<CFuint>(0));
  // displacements for data to send
  vector<int> sendDispl(m_nbProc, static_cast<CFuint>(0));
  // displacements for data to receive
  vector<int> recvDispl(m_nbProc, static_cast<CFuint>(0));

  // counts for data to send
  vector<int> sendPtrCount(m_nbProc, static_cast<CFuint>(0));
  // counts for data to receive
  vector<int> recvPtrCount(m_nbProc, static_cast<CFuint>(0));
  // displacements for data to send
  vector<int> sendPtrDispl(m_nbProc, static_cast<CFuint>(0));
  // displacements for data to receive
  vector<int> recvPtrDispl(m_nbProc, static_cast<CFuint>(0));
  
  ElementDataArray<0> tmpElem;

  // (global ID + nbNodesInElem + nbStatesInElem)*nbElements +
  // elemNode.sum() + elemState.sum()
  /// @todo nbScalarInfo will have to be changed if other element
  /// info (polyorder??) needs to be distributed
  const CFuint nbScalarInfo = 3;
  const CFuint totSendElemSize = nbScalarInfo*partSize + pdata.elemNode.size() + pdata.elemState.size();
  tmpElem.reserve(partSize, totSendElemSize);
  
  // data have to be put in the array tmpElem in sorted order
  sendDispl[(*pdata.part)[0]] = 0;
  sendPtrDispl[(*pdata.part)[0]] = 0;
  CFuint displ = 0;
  CFuint ptrDispl = 0;
  vector<CFuint> tmpElemSize(partSize, static_cast<CFuint>(0));
  
  // pdata.elemNode is ok, pdata.elemState is ok
  
  for (CFuint i = 0; i < partSize; ++i) {
    const CFuint globalElemID = mapProcToOldElemID[i];
    const CFuint localElemID = globalElemID - startElemID;
    addElement(pdata, globalElemID, localElemID, tmpElem, tmpElemSize[i]);
        
    sendCount[(*pdata.part)[i]] += tmpElemSize[i];
    cf_assert(sendCount[(*pdata.part)[i]] > 0);

    sendPtrCount[(*pdata.part)[i]]++;
    cf_assert(sendPtrCount[(*pdata.part)[i]] > 0);
    
    // case with i == 0 is already treated out of the loop
    if (i > 0) {
      if ((*pdata.part)[i] != (*pdata.part)[i-1]) {
	sendDispl[(*pdata.part)[i]] = displ;
        cf_assert(sendDispl[(*pdata.part)[i]] > 0);

	sendPtrDispl[(*pdata.part)[i]] = ptrDispl;
        cf_assert(sendPtrDispl[(*pdata.part)[i]] > 0);
      }
    }
    
    displ += tmpElemSize[i];
    ptrDispl++;
  }
  
  MPIError::getInstance().check
    ("MPI_Alltoall", "ParCFmeshFileReader::moveElementData()",
     MPI_Alltoall(&sendCount[0], 1, MPIStructDef::getMPIType(&sendCount[0]),
		  &recvCount[0], 1, MPIStructDef::getMPIType(&recvCount[0]), m_comm));
  
  MPIError::getInstance().check
    ("MPI_Alltoall", "ParCFmeshFileReader::moveElementData()",
     MPI_Alltoall(&sendPtrCount[0], 1, MPIStructDef::getMPIType(&sendPtrCount[0]),
		  &recvPtrCount[0], 1, MPIStructDef::getMPIType(&recvPtrCount[0]), m_comm));
  
  CFuint count = 0;
  CFuint countPtr = 0;
  recvDispl[0] = 0;
  recvPtrDispl[0] = 0;
  count += recvCount[0];
  countPtr += recvPtrCount[0];

  for (CFuint i = 1; i < m_nbProc; ++i) {
    if (recvCount[i] > 0) {
      recvDispl[i] = count;
      recvPtrDispl[i] = countPtr;
    }
    count += recvCount[i];
    countPtr += recvPtrCount[i];
  }
  
  const CFuint totRecvCount = std::accumulate(recvCount.begin(), recvCount.end(), 0);
  const CFuint totRecvPtrCount = std::accumulate(recvPtrCount.begin(), recvPtrCount.end(), 0);
  
  // update the array storing the number of elements per processor
  vector<CFuint> newNbElemPerProc(m_nbProc, static_cast<CFuint>(0));
  newNbElemPerProc[m_myRank] = totRecvPtrCount;
  
  MPIError::getInstance().check
    ("MPI_Allreduce", "ParCFmeshFileReader::moveElementData()",
     MPI_Allreduce(&newNbElemPerProc[0], &m_nbElemPerProc[0], m_nbProc,
		   MPIStructDef::getMPIType(&newNbElemPerProc[0]), MPI_SUM, m_comm));
  
  // deallocate PartitionerData
  SwapEmpty(pdata.elemNode);
  SwapEmpty(pdata.elemState);
  SwapEmpty(pdata.eptrn);
  SwapEmpty(pdata.eptrs);
  //  SwapEmpty(pdata.part);
  
  // allocate the storage for the local ElementDataArray<0>
  localElem.resize(totRecvPtrCount, totRecvCount);
  vector<CFuint> elemSize(totRecvPtrCount, static_cast<CFuint>(0));
  
  if (m_nbProc > 1) {
    MPIError::getInstance().check
      ("MPI_Alltoallv", "ParCFmeshFileReader::moveElementData()",
       MPI_Alltoallv(tmpElem.startData(), &sendCount[0],
		     &sendDispl[0], MPIStructDef::getMPIType(tmpElem.startData()), localElem.startData(),
		     &recvCount[0], &recvDispl[0], MPIStructDef::getMPIType(localElem.startData()), m_comm));
    
    MPIError::getInstance().check
      ("MPI_Alltoallv", "ParCFmeshFileReader::moveElementData()",
       MPI_Alltoallv(&tmpElemSize[0], &sendPtrCount[0],
		     &sendPtrDispl[0], MPIStructDef::getMPIType(&tmpElemSize[0]), &elemSize[0],
		     &recvPtrCount[0], &recvPtrDispl[0],
		     MPIStructDef::getMPIType(&elemSize[0]), m_comm));
  }
  else {
    for (CFuint i = 0; i < totRecvCount; ++i) {
      localElem.startData()[i] = tmpElem.startData()[i];
    }
    
    for (CFuint i = 0; i < tmpElemSize.size(); ++i) {
      elemSize[i] = tmpElemSize[i];
    }
  }
  
  // print some hardcore info in case of need for debugging
  CFLogDebugMin(CFPrintContainer<vector<int> >("sendCount  = ", &sendCount));
  CFLogDebugMin(CFPrintContainer<vector<int> >("recvCount  = ", &recvCount));
  CFLogDebugMin(CFPrintContainer<vector<int> >("sendDispl  = ", &sendDispl));
  CFLogDebugMin(CFPrintContainer<vector<int> >("recvDispl  = ", &recvDispl));
  CFLogDebugMin(CFPrintContainer<vector<int> >("sendPtrCount  = ", &sendPtrCount));
  CFLogDebugMin(CFPrintContainer<vector<int> >("recvPtrCount  = ", &recvPtrCount));
  CFLogDebugMin(CFPrintContainer<vector<int> >("sendPtrDispl  = ", &sendPtrDispl));
  CFLogDebugMin(CFPrintContainer<vector<int> >("recvPtrDispl  = ", &recvPtrDispl));
  
  // wipe clean the sent element array
  tmpElem.clear();
  SwapEmpty(tmpElemSize);

  // set the element ptr array  from the element sizes
  localElem.setEptrFromElemSize(elemSize);
  SwapEmpty(elemSize);

  // globally reduce the size of the element array
  vector<CFuint> tmpSizeElemArray(m_nbProc, static_cast<CFuint>(0));
  tmpSizeElemArray[m_myRank] = localElem.sizeData();
  
  vector<CFuint> sizeElemArray(m_nbProc, static_cast<CFuint>(0));
  MPI_Allreduce(&tmpSizeElemArray[0], &sizeElemArray[0], m_nbProc,
		MPIStructDef::getMPIType(&tmpSizeElemArray[0]), MPI_SUM, m_comm);
  
  /// @todo here is a possible place to do optimization
  /// if memory or speed problems arise
  set<CFuint> isLocalNode;
  set<CFuint> isLocalState;
  vector<CFuint> localNodeIDs;
  vector<CFuint> localStateIDs;
  
  // each processor builds an initial list of local node and state IDs
  setIsLocalNodeState(localElem, isLocalNode, isLocalState, localNodeIDs, localStateIDs);
    
  vector<CFuint> tmpLocalNodeIDs;
  vector<CFuint> tmpLocalStateIDs;
  vector<CFuint> localNodeIDsToRemove;
  vector<CFuint> localStateIDsToRemove;
  ElementDataArray<0> overlapElem;

  CFuint nbLocalNodeIDs = 0;
  CFuint nbLocalStateIDs = 0;

  vector<CFuint> ghostNodeIDs;
  vector<CFuint> ghostStateIDs;
  vector<bool> isOverlap;

  for (CFuint root = 0; root < m_nbProc; ++root) {
    // broadcast the element connectivity of each processor
    const CFuint rootElemSize = sizeElemArray[root];
    const CFuint rootNbElem = m_nbElemPerProc[root];
    tmpElem.resize(rootNbElem, rootElemSize);
    isOverlap.resize(rootNbElem);

    // copy the element data array, the node and state IDs list
    // of the root process
    if (m_myRank == root) {
      cf_assert(tmpElem.sizeData() == localElem.sizeData());
      cf_assert(tmpElem.getNbElements() == localElem.getNbElements());
      
      tmpElem.copy(localElem);
      nbLocalNodeIDs  = localNodeIDs.size();
      nbLocalStateIDs = localStateIDs.size();
    }
    
    MPIStruct ms;
    
    // only broadcast if you are running parallel, otherwise all data are already local
    if (m_nbProc > 1) {
      int ln[2];
      ln[0] = ln[1] = 1;
      
      MPIStructDef::buildMPIStruct(&nbLocalNodeIDs, &nbLocalStateIDs, ln, ms);
      MPI_Bcast(ms.start, 1, ms.type, root, m_comm);
    }
    
    tmpLocalNodeIDs.resize(nbLocalNodeIDs);
    tmpLocalStateIDs.resize(nbLocalStateIDs);
    
    if (m_myRank == root) {
      copy(localNodeIDs.begin(), localNodeIDs.end(), tmpLocalNodeIDs.begin());
      copy(localStateIDs.begin(), localStateIDs.end(), tmpLocalStateIDs.begin());
    }
    
    // only broadcast if you are running parallel, otherwise all data are already local
    if (m_nbProc > 1) {
      int le[4];
      le[0] = rootElemSize;
      le[1] = rootNbElem + 1;
      le[2] = nbLocalNodeIDs;
      le[3] = nbLocalStateIDs;
      
      MPIStructDef::buildMPIStruct<CFuint, CFuint, CFuint, CFuint>
	(tmpElem.startData(), tmpElem.startPtr(),
	 &tmpLocalNodeIDs[0], &tmpLocalStateIDs[0], le, ms);
      MPI_Bcast(ms.start, 1, ms.type, root, m_comm);
    }
    
    // AL: the number of overlap layers is set here
    m_nbOverLayers = (m_nbOverLayers < 2) ? 
      MeshDataStack::getActive()->getNbOverlapLayers() : m_nbOverLayers;
    
    if (m_myRank != root) {
      isOverlap.assign(isOverlap.size(),false);
      vector<CFuint> newLocalNodeIDs;
      vector<CFuint> newLocalStateIDs;
      
      for (CFuint iOver = 0; iOver < m_nbOverLayers; ++iOver) {
	updateIsLocalNodeState(root, tmpElem, overlapElem,
			       isLocalNode, isLocalState,
			       ghostNodeIDs, ghostStateIDs,
			       newLocalNodeIDs, newLocalStateIDs,
			       localNodeIDsToRemove,
			       localStateIDsToRemove,
			       isOverlap,
			       iOver+1);
      }
      
      const CFuint nbNewLocalNodeIDs = newLocalNodeIDs.size();
      for (CFuint in = 0; in < nbNewLocalNodeIDs; ++in) {
	isLocalNode.erase(newLocalNodeIDs[in]);
      }
      const CFuint nbNewLocalStateIDs = newLocalStateIDs.size();
      for (CFuint in = 0; in < nbNewLocalStateIDs; ++in) {
	isLocalState.erase(newLocalStateIDs[in]);
      }
    }
  }
  
  // add the storage of the overlap elements to the local elements
  localElem.add(overlapElem);

  // set the number of the local (= locally owned + ghost) elements
  getReadData().setNbElements(localElem.getNbElements());
  
  // remove duplicated ghost nodes
  sort(ghostNodeIDs.begin(), ghostNodeIDs.end());
  unique_copy(ghostNodeIDs.begin() , ghostNodeIDs.end(), back_inserter(m_ghostNodeIDs));

  // remove duplicated ghost states
  sort(ghostStateIDs.begin(), ghostStateIDs.end());
  unique_copy(ghostStateIDs.begin(), ghostStateIDs.end(), back_inserter(m_ghostStateIDs));

  vector<CFuint> nodeIDsToRemove;
  vector<CFuint> stateIDsToRemove;

  sort(localNodeIDsToRemove.begin(), localNodeIDsToRemove.end());
  unique_copy(localNodeIDsToRemove.begin() , localNodeIDsToRemove.end(), back_inserter(nodeIDsToRemove));

  sort(localStateIDsToRemove.begin(), localStateIDsToRemove.end());
  unique_copy(localStateIDsToRemove.begin() , localStateIDsToRemove.end(), back_inserter(stateIDsToRemove));

  const CFuint lnodesSize = localNodeIDs.size();
  m_localNodeIDs.reserve(lnodesSize - nodeIDsToRemove.size()); // this is a safe maximum size
  for (CFuint i = 0; i < lnodesSize; ++i) {
    const CFuint nodeID = localNodeIDs[i];
    if (!binary_search(localNodeIDsToRemove.begin(), localNodeIDsToRemove.end(), nodeID)) {
      m_localNodeIDs.push_back(nodeID);
    }
    else {
      isLocalNode.erase(nodeID);
    }
  }

  const CFuint lstatesSize = localStateIDs.size();
  m_localStateIDs.reserve(lstatesSize - stateIDsToRemove.size()); // this is a safe maximum size
  for (CFuint i = 0; i < lstatesSize; ++i) {
    const CFuint stateID = localStateIDs[i];
    if (!binary_search(localStateIDsToRemove.begin(), localStateIDsToRemove.end(), stateID)) {
      m_localStateIDs.push_back(stateID);
    }
    else {
      isLocalState.erase(stateID);
    }
  }
    
  CFLogDebugMax(CFPrintContainer<vector<CFuint> >("localNodeIDs  = ", &m_localNodeIDs));
  CFLogDebugMax(CFPrintContainer<vector<CFuint> >("localStateIDs = ", &m_localStateIDs));
  CFLogDebugMax(CFPrintContainer<vector<CFuint> >("ghostNodeIDs  = ", &m_ghostNodeIDs));
  CFLogDebugMax(CFPrintContainer<vector<CFuint> >("ghostStateIDs = ", &m_ghostStateIDs));
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::sortPartVec(const vector<CFuint>& globalElemIDs,
				      vector<PartitionerData::IndexT>& part,
				      CFMultiMap<PartitionerData::IndexT,CFuint>& mapProcToOldElemID)
{
  const CFuint partSize = part.size();
  
  // map the processorID with the corresponding global element IDs
  for (CFuint i = 0; i < partSize; ++i) {
    mapProcToOldElemID.insert(part[i], globalElemIDs[i]);
  }
  mapProcToOldElemID.sortKeys();
  
  // override the old pdata.part with sorted keys
  for (CFuint i = 0; i < partSize; ++i) {
    part[i] = mapProcToOldElemID.getKey(i);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::setIsLocalNodeState(ElementDataArray<0>& elem,
					      set<CFuint>& isLocalNode,
					      set<CFuint>& isLocalState,
					      vector<CFuint>& localNodeIDs,
					      vector<CFuint>& localStateIDs)
{
  CFAUTOTRACE;
  
  ElementDataArray<0>::Itr itr = elem.begin();
  ElementDataArray<0>::Itr end = elem.end();
  
  CFuint countn = 0; // counter of local nodes
  CFuint counts = 0; // counter of local states
  CFuint localElemID = 0;
  
  // counters for sanity check
  vector<CFuint> ncounter;
  if (m_nbProc == 1) {ncounter.resize(m_totNbNodes, (CFuint)0);}
  vector<CFuint> scounter;
  if (m_nbProc == 1) {scounter.resize(m_totNbStates, (CFuint)0);}
  
  for (; itr != end; ++itr, ++localElemID) {
    const CFuint nbNodesInElem = itr.get(ElementDataArray<0>::NB_NODES);
    
    // if (itr.getState(0) == 0) cout << "E[" <<localElemID <<"] has stateID=0 and nodeIDs=(";
    
    for (CFuint in = 0; in < nbNodesInElem; ++in) { 
      const CFuint nodeID = itr.getNode(in);
      cf_assert(nodeID < m_totNbNodes);
      if (m_nbProc == 1) {ncounter[nodeID]++;}
      // if (itr.getState(0) == 0) cout << nodeID << " ";   
      if (isLocalNode.count(nodeID) == 0) {
	isLocalNode.insert(nodeID);
	countn++;
      }
    }
    // if (itr.getState(0) == 0) cout << ")\n";
    
    const CFuint nbStatesInElem = itr.get(ElementDataArray<0>::NB_STATES);
    for (CFuint is = 0; is < nbStatesInElem; ++is) {
      const CFuint stateID = itr.getState(is); 
      if (m_nbProc == 1) {scounter[stateID]++;}
      cf_assert(stateID < m_totNbStates);
      if (isLocalState.count(stateID) == 0) {
	isLocalState.insert(stateID);
	counts++;
      }
    }
  }
  
  // sanity checks (to be removed)
  /*bool exitLoop = false;
  if (m_nbProc == 1) {
    bool exitLoop = false;
    for (CFuint i = 0; i < m_totNbStates; ++i) {
      if (scounter[i] == 0) {
	static CFuint countStates = 0;
	CFLog(ERROR,"ParCFmeshFileReader::setIsLocalNodeState() => " << ++countStates <<  " state [" << i << "] not found\n");
	exitLoop = true;
      }
      if (scounter[i] > 1) {
	CFLog(ERROR,"ParCFmeshFileReader::setIsLocalNodeState() => state [" << i << "] found #" << scounter[i] << " times\n");
	exitLoop = true;
      }
    }  
  
    for (CFuint i = 0; i < m_totNbNodes; ++i) {	
      if (ncounter[i] == 0) {
	CFLog(ERROR,"ParCFmeshFileReader::setIsLocalNodeState() => node [" << i << "] not found\n");
	exitLoop = true;
      }
    }
    
    if (exitLoop) exit(1);
  }*/
  
  for (set<CFuint>::const_iterator it=isLocalNode.begin(); 
       it != isLocalNode.end(); ++it) {
    localNodeIDs.push_back(*it);
  }
  for (set<CFuint>::const_iterator it=isLocalState.begin(); 
       it != isLocalState.end(); ++it) {
    localStateIDs.push_back(*it);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::updateIsLocalNodeState
(CFuint root,
 ElementDataArray<0>& elem,
 ElementDataArray<0>& overlapElem,
 set<CFuint>& isLocalNode,
 set<CFuint>& isLocalState,
 vector<CFuint>& ghostNodeIDs,
 vector<CFuint>& ghostStateIDs,
 vector<CFuint>& newLocalNodeIDs,
 vector<CFuint>& newLocalStateIDs,
 vector<CFuint>& localNodeIDsToRemove,
 vector<CFuint>& localStateIDsToRemove,
 vector<bool>& isOverlap,
 CFuint nOverlap)
{
  ElementDataArray<0>::Itr itr = elem.begin();
  ElementDataArray<0>::Itr end = elem.end();
  CFuint iElem = 0;

  for (; itr != end; ++itr, ++iElem) {
    if (!isOverlap[iElem]) {
      const CFuint nbNodesInElem = itr.get(ElementDataArray<0>::NB_NODES);
      bool sharesNode = false;
      
      for (CFuint in = 0; in < nbNodesInElem; ++in) {
	const CFuint nodeID = itr.getNode(in);
	// if the node is marked as local but is broadcast by
	// a process of lower ranking, the node is marked as
	// ghost in this process
	if (isLocalNode.count(nodeID) > 0) {
	  if (root < m_myRank && nOverlap == 1) {
	    cf_assert(nOverlap == 1);
	    ghostNodeIDs.push_back(nodeID);
	    localNodeIDsToRemove.push_back(nodeID);
	  }
	  
	  // check if the given elements have at least 1 node in common with
	  // the local elements in the current partition
	  if (!sharesNode) {
	    sharesNode = true;
	  }
	}
      }
      
      const CFuint nbStatesInElem = itr.get(ElementDataArray<0>::NB_STATES);
      for (CFuint is = 0; is < nbStatesInElem; ++is) {
	const CFuint stateID = itr.getState(is);
	// if the state is marked as local but is broadcast by
	// a process of lower ranking, the state is marked as
	// ghost in this process
	if (isLocalState.count(stateID) > 0) {
	  if (root < m_myRank && nOverlap == 1) {
	    cf_assert(nOverlap == 1);
	    ghostStateIDs.push_back(stateID);
	    localStateIDsToRemove.push_back(stateID);
	  }
	}
      }
      
      // if the current element shares at least one node with the
      // local elements in the current partition, store the current element
      // in the overlap region
      if (sharesNode) {
	overlapElem.addElement(itr);
        isOverlap[iElem] = true;
	
	// update the list of ghost nodes
	const CFuint nbNodesInElem = itr.get(ElementDataArray<0>::NB_NODES);
	for (CFuint in = 0; in < nbNodesInElem; ++in) {
	  const CFuint nodeID = itr.getNode(in);
	  if (isLocalNode.count(nodeID) == 0) {
	    ghostNodeIDs.push_back(nodeID);
	    
	    // if the number of layers of overlap is not reached
	    // consider all the current nodes as local
	    if (nOverlap < m_nbOverLayers) {
	      newLocalNodeIDs.push_back(nodeID);
	    }
	  }
	}
	
	// update the list of ghost states
	const CFuint nbStatesInElem = itr.get(ElementDataArray<0>::NB_STATES);
	for (CFuint is = 0; is < nbStatesInElem; ++is) {
	  const CFuint stateID = itr.getState(is);
	  if (isLocalState.count(stateID) == 0) {
	    ghostStateIDs.push_back(stateID);
	    
	    // if the number of layers of overlap is not reached
	    // consider all the current states as local
	    if (nOverlap < m_nbOverLayers) {
	      newLocalStateIDs.push_back(stateID);
	    }
	  }
	}
      }
    }
  }
  
  const CFuint nbNewLocalNodeIDs = newLocalNodeIDs.size();
  for (CFuint in = 0; in < nbNewLocalNodeIDs; ++in) {
    isLocalNode.insert(newLocalNodeIDs[in]);
  }
  
  const CFuint nbNewLocalStateIDs = newLocalStateIDs.size();
  for (CFuint in = 0; in < nbNewLocalStateIDs; ++in) {
    isLocalState.insert(newLocalStateIDs[in]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::addElement(const PartitionerData& pdata,
				     const CFuint globalElemID,
				     const CFuint localElemID,
				     ElementDataArray<0>& tmpElem,
				     CFuint& elemSize)
{
  cf_assert(globalElemID < m_totNbElem);
  
  // reset to 0 the element size
  elemSize = 0;
  
  tmpElem.setBeginEptr();
  const CFuint localElemIDPlus1 = localElemID + 1;
  tmpElem.addElemDataEntry(globalElemID);
  
  // AL: the following must be set to be consistent 
  // local ID (to be modified later)
  tmpElem.addElemDataEntry(globalElemID);
  // entity type ID (to be modified later)
  tmpElem.addElemDataEntry(0);
  
  const CFuint nbNodesInElem = pdata.eptrn[localElemIDPlus1]
    - pdata.eptrn[localElemID];
  cf_assert(nbNodesInElem > 1);
  tmpElem.addElemDataEntry(nbNodesInElem);

  const CFuint nbStatesInElem = pdata.eptrs[localElemIDPlus1] -
    pdata.eptrs[localElemID];
  cf_assert(nbStatesInElem > 0);
  tmpElem.addElemDataEntry(nbStatesInElem);

  // element-node connectivity
  CFuint idEN = pdata.eptrn[localElemID];
  for (CFuint in = 0; in < nbNodesInElem; ++in, ++idEN) {
    tmpElem.addElemDataEntry(pdata.elemNode[idEN]);
  }

  // element-state connectivity
  CFuint idES = pdata.eptrs[localElemID];
  for (CFuint in = 0; in < nbStatesInElem; ++in, ++idES) {
    tmpElem.addElemDataEntry(pdata.elemState[idES]);
  }
  tmpElem.setEndEptr();
  
  // set the element size
  elemSize = ElementDataArray<0>::getElementSize(nbNodesInElem, nbStatesInElem);
}
      
//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::setElements(ElementDataArray<0>& localElem)
{
  const CFuint nbLocalElems = localElem.getNbElements();
  getReadData().setNbElements(nbLocalElems);

  Common::SafePtr< vector<ElementTypeData> > elementType =
    getReadData().getElementTypeData();

  for (CFuint i = 0; i < m_totNbElemTypes; ++i) {
    (*elementType)[i].setGeoOrder(getReadData().getGeometricPolyOrder());
    (*elementType)[i].setSolOrder(getReadData().getSolutionPolyOrder());
  }

  // element type info
  vector<CFuint> elementTypeCount(m_totNbElemTypes);
  for (CFuint iType = 0; iType < m_totNbElemTypes; ++iType) {
    elementTypeCount[iType] = (*elementType)[iType].getNbElems();  
    (*elementType)[iType].setNbTotalElems((*elementType)[iType].getNbElems()); 
  }
  
  // set the total number of nodes, states and elements in the MeshData
  // to make it available later while writing
  MeshDataStack::getActive()->setTotalElementCount(elementTypeCount);
  
  // here elements are reordered by type
  vector< vector<CFuint> > elemIDPerType(m_totNbElemTypes);
  ElementDataArray<0>::Itr it;
  CFuint ne = 0;
  for (it = localElem.begin(); it != localElem.end(); ++it, ++ne) {
    const CFuint etp = getElementType(*elementType, it.get(ElementDataArray<0>::GLOBAL_ID));
    elemIDPerType[etp].push_back(ne);
  }
  cf_assert(ne == nbLocalElems);

  CFuint startIdx = 0;
  for (CFuint i = 0; i < m_totNbElemTypes; ++i) {
    (*elementType)[i].setStartIdx(startIdx);
    (*elementType)[i].setNbElems(elemIDPerType[i].size());
    startIdx += elemIDPerType[i].size();
  }
  cf_assert(startIdx == nbLocalElems);

  // set the new local element ID, taking into account the ordering by element type
  CFuint newID = 0;
  vector<CFuint> newLocalElemID(nbLocalElems);
  for (CFuint iType = 0; iType < m_totNbElemTypes; ++iType) {
    const CFuint nbEInType = elemIDPerType[iType].size();
    for (CFuint i = 0; i < nbEInType; ++i, ++newID) {
      const CFuint oldID = elemIDPerType[iType][i];
      cf_assert(oldID < nbLocalElems);
      newLocalElemID[oldID] = newID;
    }
  }
  cf_assert(newID == nbLocalElems);
  SwapEmpty(elemIDPerType);

  std::valarray<CFuint> nbCols(nbLocalElems);
  // resize element-node
  ne = 0;
  for (it = localElem.begin(); it != localElem.end(); ++it, ++ne) {
    nbCols[newLocalElemID[ne]] = it.get(ElementDataArray<0>::NB_NODES);
  }
  getReadData().resizeElementNode(nbCols);

  // resize element-state
  ne = 0;
  for (it = localElem.begin(); it != localElem.end(); ++it, ++ne) {
    nbCols[newLocalElemID[ne]] = it.get(ElementDataArray<0>::NB_STATES);
  }
  getReadData().resizeElementState(nbCols);
  
  Common::SafePtr< vector<CFuint> > globalElementIDs =
    MeshDataStack::getActive()->getGlobalElementIDs();
  globalElementIDs->resize(nbLocalElems);
  
  // fill in element-node and element-state
  m_localElemIDs.resize(nbLocalElems);
  ne = 0;
  boost::progress_display* progressBar = NULL;
  const CFuint globalRank = PE::GetPE().GetRank("Default");
  if (globalRank == 0) {
    progressBar = new boost::progress_display( localElem.getNbElements() );
  }
  
  for (it = localElem.begin(); it != localElem.end(); ++it, ++ne)
  {
    if (globalRank == 0) {
      ++(*progressBar);
    }
    
    const CFuint nbStatesInElem = it.get(ElementDataArray<0>::NB_STATES);
    cf_assert(ne < newLocalElemID.size());
    CFuint localElemID = newLocalElemID[ne];

    // if there is a single state in the element (FVM case)
    // impose that the local element ID == local state ID
    if (nbStatesInElem == 1)
    {
      localElemID = m_mapGlobToLocStateID.find(it.getState(0));
      getReadData().setElementState(localElemID, 0, localElemID);
    }
    else
    {
      for (CFuint i = 0; i < nbStatesInElem; ++i)
      {
          const CFuint stateID = m_mapGlobToLocStateID.find(it.getState(i));
          getReadData().setElementState(localElemID, i, stateID);
      }
    }
    m_localElemIDs[ne] = localElemID;
    
    // set the global element IDs list
    (*globalElementIDs)[localElemID] = getNewGlobalElementID(elementType, it.get(ElementDataArray<0>::GLOBAL_ID));
    
    const CFuint nbNodesInElem = it.get(ElementDataArray<0>::NB_NODES);
    for (CFuint i = 0; i < nbNodesInElem; ++i)
    {
      const CFuint nodeID = m_mapGlobToLocNodeID.find(it.getNode(i));
      getReadData().setElementNode(localElemID, i, nodeID);
    }
  } //  loop local elements
  
  if (m_myRank == 0) {delete progressBar;}
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::setMapGlobalToLocalID(const vector<CFuint>& localIDs,
            const vector<CFuint>& ghostIDs,
            CFMap<CFuint,CFuint>& m)
{
  const CFuint totCount = localIDs.size() + ghostIDs.size();
  vector<CFuint> allIDs;
  allIDs.reserve(totCount);

  // sort the full list of IDs
  vector<CFuint>::const_iterator it;
  for (it = localIDs.begin(); it != localIDs.end(); ++it) {
    allIDs.push_back(*it);
  }

  vector<CFuint>::const_iterator itg;
  for (itg = ghostIDs.begin(); itg != ghostIDs.end(); ++itg) {
    allIDs.push_back(*itg);
  }
  sort(allIDs.begin(), allIDs.end());

  m.reserve(totCount);
  for (CFuint i = 0; i < totCount; ++i) {
    m.insert(allIDs[i], i);
  }
  m.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::setMapNodeElemID(ElementDataArray<0>& localElem)
{
  // calculate the size of the map to be able to preallocate
  // the exact memory
  CFuint mapSize = 0; // counter of local nodes
  ElementDataArray<0>::Itr itr;
  for (itr = localElem.begin(); itr != localElem.end(); ++itr) {
    mapSize += itr.get(ElementDataArray<0>::NB_NODES);
  }
  m_mapNodeElemID.reserve(mapSize);

  CFuint elemID = 0;
  for (itr = localElem.begin(); itr != localElem.end(); ++itr, ++elemID) {
    const CFuint nbNodesInElem = itr.get(ElementDataArray<0>::NB_NODES);
    for (CFuint in = 0; in < nbNodesInElem; ++in) {
      const CFuint nodeID = itr.getNode(in);
      cf_assert(nodeID < m_totNbNodes);
      m_mapNodeElemID.insert(nodeID, m_localElemIDs[elemID]);
    }
  }
  m_mapNodeElemID.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::finish()
{
  SwapEmpty(m_localNodeIDs);
  SwapEmpty(m_localStateIDs);
  SwapEmpty(m_ghostNodeIDs);
  SwapEmpty(m_ghostStateIDs);
  m_mapGlobToLocNodeID.clear();
  m_mapGlobToLocStateID.clear();
  m_mapNodeElemID.clear();
  SwapEmpty(m_localElemIDs);

  //We dont need this anymore
  deletePtr(m_local_elem);

  //set default values
  if(getReadData().getNbGroups() == 0)
  {
    const CFuint nbGroups = 1;
    getReadData().setNbGroups(nbGroups);
    getReadData().getGroupNames()->resize(nbGroups);
    (*getReadData().getGroupNames())[0] = "InnerCells";
    getReadData().getGroupSizes()->resize(nbGroups);
    (*getReadData().getGroupSizes())[0] = getReadData().getNbElements();
    getReadData().getGroupElementLists()->resize(nbGroups);
    (*getReadData().getGroupElementLists())[0].resize(getReadData().getNbElements());

    for(CFuint iElem=0; iElem < getReadData().getNbElements(); ++iElem){
      (*getReadData().getGroupElementLists())[0][iElem] = iElem;
    }
  }
  else{
    CFout << "Nb Groups: " << getReadData().getNbGroups() << "\n";
    CFout << "Nb Groups Names: " << getReadData().getGroupNames()->size() << "\n";
    CFout << "Nb Groups Sizes: " << getReadData().getGroupSizes()->size() << "\n";
    cf_assert(getReadData().getGroupNames()->size() == getReadData().getNbGroups());
    cf_assert(getReadData().getGroupSizes()->size() == getReadData().getNbGroups());
    cf_assert(getReadData().getGroupElementLists()->size() == getReadData().getNbGroups());

    for (CFuint iGroup=0; iGroup < getReadData().getNbGroups(); ++iGroup){
      cf_assert(((*getReadData().getGroupElementLists())[iGroup]).size() == (*getReadData().getGroupSizes())[iGroup]);
    }
  }


}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbExtraVars(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbExtraVars() start\n");
  
  CFint nbExtraVars = 0;
  fin >> nbExtraVars;
  
  getReadData().setNbExtraVars
  (static_cast<CFuint>(nbExtraVars));
  
  if(nbExtraVars < 0) {
    throw BadFormatException (FromHere(),"Negative number of extra nodal variables in CFmesh");
  }
  
  CFLogDebugMin( "ParCFmeshFileReader::readNbExtraVars() end\n");
}
      
      
//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbExtraNodalVars(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbExtraNodalVars() start\n");

  CFint nbExtraNodalVars = 0;
  fin >> nbExtraNodalVars;

  getReadData().setNbExtraNodalVars
    (static_cast<CFuint>(nbExtraNodalVars));

  if(nbExtraNodalVars < 0) {
    throw BadFormatException (FromHere(),"Negative number of extra nodal variables in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshFileReader::readNbExtraNodalVars() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbExtraStateVars(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbExtraStateVars() start\n");

  CFint nbExtraStateVars = 0;
  fin >> nbExtraStateVars;

  getReadData().setNbExtraStateVars
    (static_cast<CFuint>(nbExtraStateVars));

  if(nbExtraStateVars < 0) {
    throw BadFormatException (FromHere(),"Negative number of extra state variables in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshFileReader::readNbExtraStateVars() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readExtraVarNames(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readExtraVarNames() start\n");
  
  const CFuint nbExtraVars = getReadData().getNbExtraVars();
  vector<std::string> extraVarNames(nbExtraVars);
  
  for (CFuint i = 0; i < nbExtraVars; ++i) {
    extraVarNames[i] = "";
    fin >> extraVarNames[i];
  }
  
  getReadData().setExtraVarNames(extraVarNames);
  
  CFLogDebugMin( "ParCFmeshFileReader::readExtraVarNames() end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readExtraStateVarNames(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readExtraStateVarNames() start\n");

  const CFuint nbExtraStateVars = getReadData().getNbExtraStateVars();
  vector<std::string> extraStateVarNames(nbExtraStateVars);

  for (CFuint i = 0; i < nbExtraStateVars; ++i) {
    extraStateVarNames[i] = "";
    fin >> extraStateVarNames[i];
 }

  getReadData().setExtraStateVarNames(extraStateVarNames);

  CFLogDebugMin( "ParCFmeshFileReader::readExtraStateVarNames() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readExtraNodalVarNames(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readExtraNodalVarNames() start\n");

  const CFuint nbExtraNodalVars = getReadData().getNbExtraNodalVars();
  vector<std::string> extraNodalVarNames(nbExtraNodalVars);

  for (CFuint i = 0; i < nbExtraNodalVars; ++i) {
    extraNodalVarNames[i] = "";
    fin >> extraNodalVarNames[i];
  }

  getReadData().setExtraNodalVarNames(extraNodalVarNames);

  CFLogDebugMin( "ParCFmeshFileReader::readExtraNodalVarNames() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readExtraVarStrides(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readExtraVarStrides() start\n");
  
  const CFuint nbExtraVars = getReadData().getNbExtraVars();
  vector<CFuint> extraVarStrides(nbExtraVars);
  
  for (CFuint i = 0; i < nbExtraVars; ++i) {
    fin >> extraVarStrides[i];
  }
  
  getReadData().setExtraVarStrides(extraVarStrides);
  
  CFLogDebugMin( "ParCFmeshFileReader::readExtraVarStrides() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readExtraStateVarStrides(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readExtraStateVarStrides() start\n");

  const CFuint nbExtraStateVars = getReadData().getNbExtraStateVars();
  vector<CFuint> extraStateVarStrides(nbExtraStateVars);

  for (CFuint i = 0; i < nbExtraStateVars; ++i) {
    fin >> extraStateVarStrides[i];
  }

  getReadData().setExtraStateVarStrides(extraStateVarStrides);

  CFLogDebugMin( "ParCFmeshFileReader::readExtraStateVarStrides() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readExtraNodalVarStrides(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readExtraNodalVarStrides() start\n");

  const CFuint nbExtraNodalVars = getReadData().getNbExtraNodalVars();
  vector<CFuint> extraNodalVarStrides(nbExtraNodalVars);

  for (CFuint i = 0; i < nbExtraNodalVars; ++i) {
    fin >> extraNodalVarStrides[i];
  }

  getReadData().setExtraNodalVarStrides(extraNodalVarStrides);

  CFLogDebugMin( "ParCFmeshFileReader::readExtraNodalVarStrides() end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readExtraVars(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readExtraVars() start\n");
  getReadData().resizeExtraVars();
  const CFuint nbExtraVars = getReadData().getNbExtraVars();
  const vector<CFuint>& extraVarStrides = *getReadData().getExtraVarStrides();
  RealVector extraVars;
  if (nbExtraVars > 0) {
    const CFuint sizeExtraVars = std::accumulate(extraVarStrides.begin(),
                                                 extraVarStrides.end(), 0);
    extraVars.resize(sizeExtraVars);
    extraVars = 0.0;
    
    getReadData().prepareExtraVars();
    
    for (CFuint i = 0; i < sizeExtraVars; ++i) {
      fin >> extraVars[i];
    }  
    
    getReadData().setExtraVar(extraVars);
    
  }
  

  
  CFLogDebugMin( "ParCFmeshFileReader::readExtraNodalVarStrides() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readNbGroups(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readNbGroups() start\n");

  CFuint nbGroups = 0;
  fin >> nbGroups;

  if (nbGroups < 1)
  {
    throw BadFormatException (FromHere(),"Number of nbGroups in file must be at least 1");
  }

  getReadData().setNbGroups(nbGroups);

  CFLogDebugMin( "ParCFmeshFileReader::readNbTRSs() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readGroupName(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readGroupName() start\n");

  std::string name = "";
  fin >> name;

  if(getReadData().getNbGroups() == 1) name = "InnerCells";
  CFLogDebugMin( "Creating Group " + name + "\n");

  const CFuint idx = m_groups_idxmap.size();
  m_groups_idxmap[name] = idx;
  getReadData().getGroupNames()->push_back(name);

  // set the current TRS being read
  m_curr_group = name;
  // index
  CFLogDebugMin( "Current TRS " << m_curr_group << " has index " << m_groups_idxmap[m_curr_group] << "\n");

  CFLogDebugMin( "ParCFmeshFileReader::readGroupName() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readGroupElementNb(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readGroupElementNb() start\n");

  CFint nbgeo = 0;
  fin >> nbgeo;

  cf_assert(nbgeo > 0);

  getReadData().getGroupSizes()->push_back(nbgeo);

  CFLogDebugMin( "ParCFmeshFileReader::readGroupElementNb() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshFileReader::readGroupElementList(ifstream& fin)
{
  CFLogDebugMin( "ParCFmeshFileReader::readGroupElementList() start\n");

  getReadData().getGroupElementLists()->resize(getReadData().getNbGroups());

  const CFuint nbLocalElements = m_localElemIDs.size();
  const CFuint idx = m_groups_idxmap[m_curr_group];
  const CFuint totalNbGroupElems = (*getReadData().getGroupSizes())[idx];

  //global to local mapping in a map
  Common::CFMap<CFuint,CFuint> global2LocalID;
  global2LocalID.reserve(nbLocalElements);

  ElementDataArray<0>::Itr itr = m_local_elem->begin();
  ElementDataArray<0>::Itr end = m_local_elem->end();

  CFuint localElemID = 0;
  for (; itr != end; ++itr, ++localElemID) {
    global2LocalID.insert(itr.get(ElementDataArray<0>::GLOBAL_ID),localElemID);
  }
  global2LocalID.sortKeys();

  // for each element of the group,
  // check if the local processor has the element
  CFuint globalElementID;
  for(CFuint iElem = 0 ; iElem < totalNbGroupElems; ++iElem){
    fin >> globalElementID;

    //our elements start by 0 while Gambit starts with 1
    globalElementID -= 1;
    //if yes, store it - otherwise continue reading
    if(global2LocalID.exists(globalElementID)){
      const CFuint localElemID = global2LocalID.find(globalElementID);
      ((*getReadData().getGroupElementLists())[idx]).push_back(localElemID);
    }
  }
  const CFuint nbGroupLocalElements = ((*getReadData().getGroupElementLists())[idx]).size();

  (*getReadData().getGroupSizes())[idx] = nbGroupLocalElements;

  CFLogDebugMin( "ParCFmeshFileReader::readGroupElementList() end\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
