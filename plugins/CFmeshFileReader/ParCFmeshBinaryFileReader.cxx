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
#include "Common/MPI/MPIIOFunctions.hh"

#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/NamespaceSwitcher.hh"
#include "Framework/MeshData.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/MeshPartitioner.hh"

#include "CFmeshFileReader/ParCFmeshBinaryFileReader.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

ParCFmeshBinaryFileReader::ParCFmeshBinaryFileReader() :
  ParCFmeshFileReader(), // here you have to pass the name of this object to configure
  m_mapString2ReaderFun(),
  m_fh(),
  m_status()
{
  addConfigOptionsTo(this);
  
  m_maxBuffSize = 2147479200; // (CFuint) std::numeric_limits<int>::max();
  setParameter("MaxBuffSize",&m_maxBuffSize);
}

//////////////////////////////////////////////////////////////////////////////

ParCFmeshBinaryFileReader::~ParCFmeshBinaryFileReader()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< int >("MaxBuffSize", "Maximum buffer size for MPI I/O");
}
 
/////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::setup()
{
  ParCFmeshFileReader::setup();
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::setMapString2Readers()
{
  m_mapString2ReaderFun["!COOLFLUID_VERSION"]     = &ParCFmeshBinaryFileReader::readCFVersion;
  m_mapString2ReaderFun["!COOLFLUID_SVNVERSION"]  = &ParCFmeshBinaryFileReader::readSvnVersion;
  m_mapString2ReaderFun["!CFMESH_FORMAT_VERSION"] = &ParCFmeshBinaryFileReader::readCFmeshVersion;
  m_mapString2ReaderFun["!NB_DIM"]             = &ParCFmeshBinaryFileReader::readDimension;
  m_mapString2ReaderFun["!NB_EQ"]              = &ParCFmeshBinaryFileReader::readNbEquations;
  m_mapString2ReaderFun["!NB_NODES"]           = &ParCFmeshBinaryFileReader::readNbNodes;
  m_mapString2ReaderFun["!NB_STATES"]          = &ParCFmeshBinaryFileReader::readNbStates;
  m_mapString2ReaderFun["!STORE_PASTSTATES"]   = &ParCFmeshBinaryFileReader::readStorePastStates;
  m_mapString2ReaderFun["!STORE_PASTNODES"]    = &ParCFmeshBinaryFileReader::readStorePastNodes;
  m_mapString2ReaderFun["!STORE_INTERSTATES"]   = &ParCFmeshBinaryFileReader::readStoreInterStates;
  m_mapString2ReaderFun["!STORE_INTERNODES"]    = &ParCFmeshBinaryFileReader::readStoreInterNodes;
  m_mapString2ReaderFun["!NB_EXTRA_SVARS"]     = &ParCFmeshBinaryFileReader::readNbExtraStateVars;
  m_mapString2ReaderFun["!NB_EXTRA_NVARS"]     = &ParCFmeshBinaryFileReader::readNbExtraNodalVars;
  m_mapString2ReaderFun["!NB_EXTRA_VARS"]     = &ParCFmeshBinaryFileReader::readNbExtraVars;
  m_mapString2ReaderFun["!EXTRA_NVARS_NAMES"]  = &ParCFmeshBinaryFileReader::readExtraNodalVarNames;
  m_mapString2ReaderFun["!EXTRA_SVARS_NAMES"]  = &ParCFmeshBinaryFileReader::readExtraStateVarNames;
  m_mapString2ReaderFun["!EXTRA_VARS_NAMES"]  = &ParCFmeshBinaryFileReader::readExtraVarNames;
  m_mapString2ReaderFun["!EXTRA_NVARS_STRIDES"]  = &ParCFmeshBinaryFileReader::readExtraNodalVarStrides;
  m_mapString2ReaderFun["!EXTRA_SVARS_STRIDES"]  = &ParCFmeshBinaryFileReader::readExtraStateVarStrides;
  m_mapString2ReaderFun["!EXTRA_VARS_STRIDES"]  = &ParCFmeshBinaryFileReader::readExtraVarStrides;
  m_mapString2ReaderFun["!EXTRA_VARS"]         = &ParCFmeshBinaryFileReader::readExtraVars;
  m_mapString2ReaderFun["!NB_ELEM"]            = &ParCFmeshBinaryFileReader::readNbElements;
  m_mapString2ReaderFun["!NB_ELEM_TYPES"]      = &ParCFmeshBinaryFileReader::readNbElementTypes;
  m_mapString2ReaderFun["!GEOM_POLYTYPE"]      = &ParCFmeshBinaryFileReader::readGeometricPolyType;
  m_mapString2ReaderFun["!SOL_POLYTYPE"]       = &ParCFmeshBinaryFileReader::readSolutionPolyType;
  m_mapString2ReaderFun["!GEOM_POLYORDER"]     = &ParCFmeshBinaryFileReader::readGeometricPolyOrder;
  m_mapString2ReaderFun["!SOL_POLYORDER"]      = &ParCFmeshBinaryFileReader::readSolutionPolyOrder;
  m_mapString2ReaderFun["!LIST_NODE"]          = &ParCFmeshBinaryFileReader::readNodeList;
  m_mapString2ReaderFun["!LIST_STATE"]         = &ParCFmeshBinaryFileReader::readStateList;
  m_mapString2ReaderFun["!NB_TRSs"]            = &ParCFmeshBinaryFileReader::readNbTRSs;
  m_mapString2ReaderFun["!TRS_NAME"]           = &ParCFmeshBinaryFileReader::readTRSName;
  m_mapString2ReaderFun["!NB_TRs"]             = &ParCFmeshBinaryFileReader::readNbTRs;
  m_mapString2ReaderFun["!NB_GEOM_ENTS"]       = &ParCFmeshBinaryFileReader::readNbGeomEnts;
  m_mapString2ReaderFun["!GEOM_TYPE"]          = &ParCFmeshBinaryFileReader::readGeomType;
  m_mapString2ReaderFun["!LIST_GEOM_ENT"]      = &ParCFmeshBinaryFileReader::readGeomEntList;
  m_mapString2ReaderFun["!ELEM_TYPES"]         = &ParCFmeshBinaryFileReader::readElementTypes;
  m_mapString2ReaderFun["!NB_ELEM_PER_TYPE"]   = &ParCFmeshBinaryFileReader::readNbElementsPerType;
  m_mapString2ReaderFun["!NB_NODES_PER_TYPE"]  = &ParCFmeshBinaryFileReader::readNbNodesPerType;
  m_mapString2ReaderFun["!NB_STATES_PER_TYPE"] = &ParCFmeshBinaryFileReader::readNbStatesPerType;
  m_mapString2ReaderFun["!NB_GROUPS"]          = &ParCFmeshBinaryFileReader::readNbGroups;
  m_mapString2ReaderFun["!GROUP_NAME"]         = &ParCFmeshBinaryFileReader::readGroupName;
  m_mapString2ReaderFun["!GROUP_ELEM_NB"]      = &ParCFmeshBinaryFileReader::readGroupElementNb;
  m_mapString2ReaderFun["!GROUP_ELEM_LIST"]    = &ParCFmeshBinaryFileReader::readGroupElementList;
  m_mapString2ReaderFun["!LIST_ELEM"]          = &ParCFmeshBinaryFileReader::readElementList;
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readCFVersion(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readCFVersion() start\n");
  
  string version = MPIIOFunctions::readAndTrimString(fh);
  CFLog(INFO, "CF version : " << version << "\n");
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readCFVersion() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readSvnVersion(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readSvnVersion() start\n");
  
  string version = MPIIOFunctions::readAndTrimString(fh);
  CFLog(INFO, "SVN version : " << version << "\n");
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readSvnVersion() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readCFmeshVersion(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readCFmeshVersion() start\n");
  
  string version = MPIIOFunctions::readAndTrimString(fh);
  CFLog(INFO, "CFmesh version : " << version << "\n");
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readCFmeshVersion() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readDimension(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readDimension() start\n");

  CFint dimension = 0;
  MPIIOFunctions::readScalar(fh, dimension);
  
  getReadData().setDimension(static_cast<CFuint>(dimension));
  
  if(getReadData().getDimension() != DIM_1D &&
     getReadData().getDimension() != DIM_2D &&
     getReadData().getDimension() != DIM_3D)
  {
    throw BadFormatException (FromHere(),"Wrong dimension in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshBinaryFileReader::readDimension() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbEquations(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbEquations() start\n");

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
  
  MPIIOFunctions::readScalar(fh, m_originalNbEqs);
  cf_assert(m_originalNbEqs > 0);
  CFuint nbEquations = m_originalNbEqs;
  
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
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbEquations() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbNodes(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbNodes() start\n");
  
  MPIIOFunctions::readScalar(fh, m_totNbNodes);
  CFint nbNonUpdatableNodes = 0;
  MPIIOFunctions::readScalar(fh, nbNonUpdatableNodes);
  
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

  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbNodes() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbStates(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbStates() start\n");
  
  MPIIOFunctions::readScalar(fh, m_totNbStates);
  CFint nbNonUpdatableStates = 0;
  MPIIOFunctions::readScalar(fh, nbNonUpdatableStates);
  
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
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbStates() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readStorePastStates(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readStorePastStates() start\n");
  
  CFint flag = 0;
  MPIIOFunctions::readScalar(fh, flag);
  m_hasPastStates = (bool)flag;
  
  if(getReadData().storePastStates())
  {
    if(!m_hasPastStates){
      throw BadFormatException (FromHere(),"PastStates asked but missing in CFmesh");
    }
  }

  CFLogDebugMin( "ParCFmeshBinaryFileReader::readStorePastStates() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readStorePastNodes(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readStorePastNodes() start\n");
 
  CFint flag = 0;
  MPIIOFunctions::readScalar(fh, flag);
  m_hasPastNodes = (bool)flag;

  if(getReadData().storePastNodes())
  {
    if(!m_hasPastNodes){
      throw BadFormatException (FromHere(),"PastNodes asked but missing in CFmesh");
    }
  }

  CFLogDebugMin( "ParCFmeshBinaryFileReader::readStorePastNodes() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readStoreInterStates(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readStoreInterStates() start\n");
  
  CFint flag = 0;
  MPIIOFunctions::readScalar(fh, flag);
  m_hasInterStates = (bool)flag;
  
  if(getReadData().storeInterStates())
  {
    if(!m_hasInterStates){
      throw BadFormatException (FromHere(),"InterStates asked but missing in CFmesh");
    }
  }

  CFLogDebugMin( "ParCFmeshBinaryFileReader::readStoreInterStates() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readStoreInterNodes(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readStoreInterNodes() start\n");
  
  CFint flag = 0;
  MPIIOFunctions::readScalar(fh, flag);
  m_hasInterNodes = (bool)flag;
  
  if(getReadData().storeInterNodes())
  {
    if(!m_hasInterNodes){
      throw BadFormatException (FromHere(),"InterNodes asked but missing in CFmesh");
    }
  }

  CFLogDebugMin( "ParCFmeshBinaryFileReader::readStoreInterNodes() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbElements(MPI_File* fh)
{
  CFLog(NOTICE,"Memory Usage before assembling connectivity: " << Common::OSystem::getInstance().getProcessInfo()->memoryUsage() << "\n");

  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbElements() start\n");
  
  MPIIOFunctions::readScalar(fh, m_totNbElem);
  
  getReadData().setNbElements(m_totNbElem);
  
  if(m_totNbElem < 1) {
    throw BadFormatException (FromHere(),"Number of elements < 1 in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbElements() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbElementTypes(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbElementTypes() start\n");
 
  MPIIOFunctions::readScalar(fh, m_totNbElemTypes);
  
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

  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbElementTypes() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readGeometricPolyOrder(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGeometricPolyOrder() start\n");
  
  CFuint order = 0;
  MPIIOFunctions::readScalar(fh, order);
  getReadData().setGeometricPolyOrder(static_cast<CFPolyOrder::Type>(order));
  
  if(getReadData().getGeometricPolyOrder() >= CFPolyOrder::MAXORDER ||
     getReadData().getGeometricPolyOrder() < CFPolyOrder::ORDER0) {
    throw BadFormatException (FromHere(),"Bad polynomial order of geometry in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGeometricPolyOrder() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readSolutionPolyOrder(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readSolutionPolyOrder() start\n");
  
  CFuint order = 0;
  MPIIOFunctions::readScalar(fh, order);
  getReadData().setSolutionPolyOrder(static_cast<CFPolyOrder::Type>(order));
  
  if(getReadData().getSolutionPolyOrder() >= CFPolyOrder::MAXORDER ||
     getReadData().getSolutionPolyOrder() < CFPolyOrder::ORDER0) {
    throw BadFormatException (FromHere(),"Bad polynomial order of solution in CFmesh");
  }
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readSolutionPolyOrder() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readGeometricPolyType(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGeometricPolyType() start\n");
  
  CFuint Type = 0;
  MPIIOFunctions::readScalar(fh, Type);
  getReadData().setGeometricPolyType(static_cast<CFPolyForm::Type>(Type));
  
  /// @todo should add check here for validity
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGeometricPolyType() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readSolutionPolyType(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readSolutionPolyType() start\n");

  CFuint Type = 0;
  MPIIOFunctions::readScalar(fh, Type);
  getReadData().setSolutionPolyType(static_cast<CFPolyForm::Type>(Type));
  
  /// @todo should add check here for validity
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readSolutionPolyType() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readElementTypes(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readElementTypes() start\n");
  
  const CFuint nbElementTypes = getReadData().getNbElementTypes();
  cf_assert(nbElementTypes > 0);
  
  SafePtr< vector<ElementTypeData> > elementType = 
    getReadData().getElementTypeData();
  elementType->resize(nbElementTypes);
  
  for (CFuint i = 0; i < nbElementTypes; ++i) {
    const std::string elementShape = MPIIOFunctions::readAndTrimString(fh);
    (*elementType)[i].setShape(elementShape);
    (*elementType)[i].setGeoShape( CFGeoShape::Convert::to_enum(elementShape) );
  }
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readElementTypes() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbElementsPerType(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readElementTypes() start\n");

  const CFuint nbElementTypes = getReadData().getNbElementTypes();
  cf_assert(nbElementTypes > 0);
  
  SafePtr< vector<ElementTypeData> > elementType = 
    getReadData().getElementTypeData();
  
  CFuint sumNbElems = 0;
  CFint nbElemPerType = 0;
  
  for (CFuint i = 0; i < nbElementTypes; ++i) {
    MPIIOFunctions::readScalar(fh, nbElemPerType);
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
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readElementTypes() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbNodesPerType(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbNodesPerType() start\n");

  const CFuint nbElementTypes = getReadData().getNbElementTypes();
  cf_assert(nbElementTypes > 0);
  
  Common::SafePtr< vector<ElementTypeData> > elementType =
    getReadData().getElementTypeData();
  
  CFint nbNodesPerType = 0;
  for (CFuint i = 0; i < nbElementTypes; ++i) {
    MPIIOFunctions::readScalar(fh, nbNodesPerType);
    
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
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbNodesPerType() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbStatesPerType(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbStatesPerType() start\n");

  const CFuint nbElementTypes = getReadData().getNbElementTypes();
  cf_assert(nbElementTypes > 0);
  
  Common::SafePtr< vector<ElementTypeData> > elementType =
    getReadData().getElementTypeData();
  
  CFint nbStatesPerType = 0;
  for (CFuint i = 0; i < nbElementTypes; ++i) {
    MPIIOFunctions::readScalar(fh, nbStatesPerType);
    
    (*elementType)[i].setNbStates
      (static_cast<CFuint>(nbStatesPerType));
    
    CFLogDebugMin( "nbStatesPerType = " << nbStatesPerType << "\n");
    
    if(nbStatesPerType < 0) {
      throw BadFormatException
	(FromHere(), "In CFmesh, negative number of states per type");
    }
  }
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbStatesPerType() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readElementList(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readElementList() start\n");
  
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
  
  readElemListRank(pdata, fh);
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
    
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readElementList() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readElemListRank(PartitionerData& pdata,
						 MPI_File* fh)
{
  CFuint start = 0;
  for (CFuint rank = 0; rank < m_myRank; ++rank) {
    start += m_nbElemPerProc[rank];
  }
  const CFuint ne = m_nbElemPerProc[m_myRank];
  const CFuint end = start + ne;
  
  vector<PartitionerData::IndexT>& eNode  = pdata.elemNode;
  vector<PartitionerData::IndexT>& eState = pdata.elemState;
  vector<PartitionerData::IndexT>& eptrn  = pdata.eptrn;
  vector<PartitionerData::IndexT>& eptrs  = pdata.eptrs;
  
  SafePtr< vector<ElementTypeData> > elementType =
    getReadData().getElementTypeData();
  
  const CFuint nbEEminProc       = m_nbElemPerProc[m_myRank];
  const CFuint nbEEminProcMinus1 = nbEEminProc - 1;
  CFuint ncount = 0;
  CFuint scount = 0;
  CFuint ipos = 0;
  CFuint iElemBegin = 0;
  CFuint localElemSize = 0;
  CFuint elemOffset = 0;
  CFuint elemMaxOffset = 0;
  
  for (CFuint iType = 0; iType < m_totNbElemTypes; ++iType) {
    const CFuint nbNodesInElem  = (*elementType)[iType].getNbNodes();
    const CFuint nbStatesInElem = (*elementType)[iType].getNbStates();
    const CFuint nbNodesStates  = nbNodesInElem + nbStatesInElem;
    const CFuint nbElementsPerType = (*elementType)[iType].getNbElems();
    const CFuint iElemEnd = iElemBegin + nbElementsPerType;
    elemMaxOffset += nbElementsPerType*nbNodesStates;
    if (iElemBegin > end) continue;
    
    // loop over the elements in this type
    for (CFuint iElem = iElemBegin; iElem < iElemEnd; ++iElem) {
      if (iElem >= start && iElem < end) {
	localElemSize += nbNodesStates;
	eptrn[ipos] = ncount;
	eptrs[ipos] = scount;
	
	if (ipos == nbEEminProcMinus1) {
	  eptrn[nbEEminProc] = eptrn[ipos] + nbNodesInElem;
	  eptrs[nbEEminProc] = eptrs[ipos] + nbStatesInElem;
	}
	
	++ipos;
	ncount += nbNodesInElem;
	scount += nbStatesInElem;
      }
      else if (iElem < start) {
	elemOffset += nbNodesStates;
      }
      else if (iElem >= end) {
	break;
      }
    }
    
    iElemBegin += nbElementsPerType;
  }
  
  cf_assert(elemOffset < elemMaxOffset);
  elemOffset    *= sizeof(CFuint);
  elemMaxOffset *= sizeof(CFuint);
  
  CFLog(VERBOSE, "P["<< m_myRank << "] => [start, end] = [" << start << ", " << end << "]\n");
  cf_assert(pdata.elemNode.size() + pdata.elemState.size() == localElemSize);
  
  vector<CFuint> buf(localElemSize);
  MPI_Offset offset;
  MPI_File_get_position(*fh, &offset);
  MPI_Offset startPos = offset + elemOffset + 1;    // the "1" is for the character "\n" 
  MPI_Offset endPos   = offset + elemMaxOffset + 1; // the "1" is for the character "\n" 
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileReader::readElemListRank() => elements read in position [" << offset + 1 << ", " << endPos << "]\n");
  CFLog(VERBOSE, "P[" << m_myRank << "] reads " << localElemSize  << " elements starting from position " << startPos << "\n");
    
  MPIIOFunctions::readAll("ParCFmeshBinaryFileReader::readElemListRank()", fh, startPos, &buf[0], (CFuint)buf.size(), m_maxBuffSize);
  
  CFLog(DEBUG_MAX, CFPrintContainer<vector<CFuint> >("buf = ", &buf));
  
  CFuint counter = 0;
  const CFuint nbLocalElem = eptrn.size()-1;
  for (CFuint i = 0; i < nbLocalElem; ++i) {
    const CFuint nncount = eptrn[i];
    const CFuint sscount = eptrs[i];
    const CFuint nbNodesInElem  = eptrn[i+1]-nncount;
    const CFuint nbStatesInElem = eptrs[i+1]-sscount;
    const CFuint nbNodeStates   = nbNodesInElem + nbStatesInElem; 
    
    // store element-node connectivity
    CFLog(DEBUG_MAX, "Element[" << i << "] => element-nodeIDs = ");
    for (CFuint iNode = 0; iNode < nbNodesInElem; ++iNode, ++counter) {
      cf_assert(counter < buf.size()); 
      const CFuint nodeID = buf[counter];
      CFLog(DEBUG_MAX, " " << nodeID << " ");
      checkDofID("node", i, iNode, nodeID, m_totNbNodes);
      cf_assert(nncount + iNode < eNode.size()); 
      eNode[nncount + iNode] = nodeID;
    }
    
    CFLog(DEBUG_MAX, ", element-stateIDs = ");
    // store element-state connectivity
    for (CFuint iState = 0; iState < nbStatesInElem; ++iState, ++counter) {
      cf_assert(counter < buf.size());
      const CFuint stateID = buf[counter];
      CFLog(DEBUG_MAX, " " << stateID << " ");
      checkDofID("state", i, iState, stateID, m_totNbStates);
      cf_assert(sscount + iState < eState.size()); 
      eState[sscount + iState] = stateID;
    }
    CFLog(DEBUG_MAX, "\n"); 
  }
  cf_assert(counter == buf.size());
  
  MPI_Barrier(m_comm);
  MPI_File_seek(*fh, endPos, MPI_SEEK_SET);
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbTRSs(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbTRSs() start\n");
  
  CFuint nbTRSs = 0;
  MPIIOFunctions::readScalar(fh, nbTRSs);
  
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
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbTRSs() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readTRSName(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readTRSName() start\n");
  
  std::string name = MPIIOFunctions::readAndTrimString(fh);
  CFLog(VERBOSE, "Found TRS " + name + "\n");
  
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
    
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readTRSName() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbTRs(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbTRs() start\n");

  CFuint nbTRsInTRS = 0;
  MPIIOFunctions::readScalar(fh, nbTRsInTRS);
  
  cf_assert(nbTRsInTRS > 0);

  const CFuint idx = m_trs_idxmap[m_curr_trs];
  (*getReadData().getNbTRs())[idx] += nbTRsInTRS;
  
  CFLog(VERBOSE, "TRS " << m_curr_trs << " + " << nbTRsInTRS << "TR, total " << (*getReadData().getNbTRs())[idx] << " TR\n");
  
  // set the current number of TR's being read
  m_curr_nbtr = nbTRsInTRS;
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbTRs() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbGeomEnts(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbGeomEnts() start\n");
  
  const CFuint idx = m_trs_idxmap[m_curr_trs];
  for (CFuint iTr = 0; iTr < m_curr_nbtr; ++iTr) {
    CFint nbgeo = 0;
    MPIIOFunctions::readScalar(fh, nbgeo);
    cf_assert(nbgeo > 0);
    
    (*getReadData().getNbGeomEntsPerTR())[idx].push_back(nbgeo);
    MeshDataStack::getActive()->getTotalTRSInfo()[idx].push_back(nbgeo);
  }
    
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbGeomEnts() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readGeomType(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGeomType() start\n");
  
  const std::string geomName = MPIIOFunctions::readAndTrimString(fh);
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
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGeomType() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readGeomEntList(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGeomEntList() start\n");
 
  SafePtr< vector<CFuint> > nbTRs = getReadData().getNbTRs();
  
  // need to read maximum number of nodes and states in GEO for each TR in
  // this TRS for preallocation purposes
  // must be written down by the writer
  // AL: this is an important format change compared to ASCII format
  
  // read in the maximum number of nodes/states for each TR geos
  CFMat<CFuint> nbNodesStatesInTRGeo(m_curr_nbtr, 2, static_cast<CFuint>(0));
  MPIIOFunctions::readArray(fh, &nbNodesStatesInTRGeo[0], nbNodesStatesInTRGeo.size());
  char c; MPIIOFunctions::readScalar(fh,c);  // reads a "\n"
    
  typedef CFMultiMap<CFuint,CFuint>::MapIterator MapItr;
  
  // load TRs data into memory for further use
  // this is actually only useful to be able to write file
  // without having constructed TRSs
  // get id of last TRS read from file
  SafePtr< vector<vector<CFuint> > > nbGeomEntsPerTR =
    getReadData().getNbGeomEntsPerTR();
  
  const CFuint iTRS = m_trs_idxmap[m_curr_trs];
  const CFuint nbTRsAdded = (*nbTRs)[iTRS];
  
  getReadData().resizeGeoConn(iTRS, nbTRsAdded);
  
  CFLogDebugMin("Rank " << m_myRank << " CurrTRS    = " << m_curr_trs << "\n");
  CFLogDebugMin("Rank " << m_myRank << " iTRS       = " << iTRS << "\n");
  CFLogDebugMin("Rank " << m_myRank << " nbTRsAdded = " << nbTRsAdded << "\n");
  
  // Set some global info
  SafePtr<vector<vector<vector<CFuint> > > > trsGlobalIDs =
    MeshDataStack::getActive()->getGlobalTRSGeoIDs();
  (*trsGlobalIDs)[iTRS].resize(nbTRsAdded);
  
  CFuint sizeBuf = 0;
  // loop only in the new TRs, which have not been read yet
  for (CFuint iTR = nbTRsAdded - m_curr_nbtr; iTR < nbTRsAdded; ++iTR) {
    const CFuint nbTRGeos = (*nbGeomEntsPerTR)[iTRS][iTR];
    sizeBuf += nbTRGeos*(2 + nbNodesStatesInTRGeo(iTR,0) + nbNodesStatesInTRGeo(iTR,1));
  }
  
  // read the full list of geometric entities
  vector<CFint> buf(sizeBuf);
  MPIIOFunctions::readArray(fh, &buf[0], sizeBuf);
  
  pair<std::valarray<CFuint>, std::valarray<CFuint> > geoConLocal;
  
  CFuint counter = 0;
  // loop only in the new TRs, which have not been read yet
  for (CFuint iTR = nbTRsAdded - m_curr_nbtr; iTR < nbTRsAdded; ++iTR)
  {
    const CFuint nbTRGeos = (*nbGeomEntsPerTR)[iTRS][iTR];
    CFLogDebugMin("Rank " << m_myRank << " iTR = " << iTR << ", nbTRGeos = "<< nbTRGeos << "\n");
    const CFuint stride = 2 + nbNodesStatesInTRGeo(iTR,0) + nbNodesStatesInTRGeo(iTR,1);
    
    CFuint countGeos = 0;
    for (CFuint iGeo = 0; iGeo < nbTRGeos; ++iGeo, counter+=stride)
    {
      const CFuint nbNodesInGeo  = buf[counter];
      cf_assert(counter < buf.size());
      const CFuint nbStatesInGeo = buf[counter+1];
      cf_assert(counter+1 < buf.size());
      cf_assert(nbNodesInGeo  == 2 || nbNodesInGeo == 3 || nbNodesInGeo == 4);
      cf_assert(nbStatesInGeo == 1);
      
      geoConLocal.first.resize(nbNodesInGeo);
      geoConLocal.second.resize(nbStatesInGeo);
      
      const CFuint startn = counter + 2;
      for(CFuint n = 0; n < nbNodesInGeo; ++n)
      {
	cf_assert(startn+n < buf.size());
	geoConLocal.first[n] = buf[startn+n];
	if (geoConLocal.first[n] >= m_totNbNodes) {
	  CFLog(ERROR, "ParCFmeshBinaryFileReader::readGeomEntList() => TRD nodeID " 
		<< geoConLocal.first[n] << " >= " << m_totNbNodes << "\n");
	  cf_assert(geoConLocal.first[n] < m_totNbNodes);
	}
      }

      const CFuint starts = startn + nbNodesInGeo;
      for(CFuint s = 0; s < nbStatesInGeo; ++s)
      {
	cf_assert(starts+s < buf.size());
	geoConLocal.second[s] = buf[starts + s];
	
	if (geoConLocal.second[s] >= m_totNbStates) {
	  CFLog(ERROR, "ParCFmeshBinaryFileReader::readGeomEntList() => TRD stateID " 
		<< geoConLocal.second[s] << " >= " << m_totNbStates << "\n");
	  cf_assert(geoConLocal.second[s] < m_totNbStates);
	}
      }
      
      // check if the global ID of the first node of the
      // geometric entity is referenced by any local element
      bool nodeFound = false;
      pair<MapItr, MapItr> etr = m_mapNodeElemID.find(geoConLocal.first[0],nodeFound);
      
      // if the first one is found, check if all the other
      // GE nodes are referenced by one amongst all vertex-neighbor elements
      if (nodeFound) {
        bool exitLoop = false;
        for (MapItr etm = etr.first; (etm != etr.second) && (!exitLoop); ++etm) {
          const CFuint localElemID = etm->second;
          const CFuint nbENodes = getReadData().getNbNodesInElement(localElemID);
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
              const CFuint nodeID = getReadData().getElementNode(localElemID, jn);
	      if (nodeID == localNodeID) {
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
    
    CFLog(VERBOSE, "Rank " << m_myRank << ", iTR = " << iTR << ", countGeos = " << countGeos << "\n");
  }
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGeomEntList() end\n");
}

//////////////////////////////////////////////////////////////////////////////

bool ParCFmeshBinaryFileReader::readString(MPI_File* fh)
{
  string key = MPIIOFunctions::readAndTrimString(fh);
  
  // check end of file
  if (key != getReaderTerminator()) {
    CFLog(VERBOSE, "ParCFmeshBinaryFileReader::readString() => key " << key << " before\n");
    MapString2ReaderP::iterator key_pair = m_mapString2ReaderFun.find(key);
    CFLog(VERBOSE, "ParCFmeshBinaryFileReader::readString() => key " << key << " after\n");
    
    // check if key exists else ignore it
    if (key_pair != m_mapString2ReaderFun.end()) {
      CFLogDebugMin( "Read CFmesh Key: " << key << "\n");
      ReaderFun function = key_pair->second;
      cf_assert(function != CFNULL);
      (this->*function)(fh);
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
    
void ParCFmeshBinaryFileReader::finish()
{
  ParCFmeshFileReader::finish();
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbExtraVars(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbExtraVars() start\n");
 
  CFint nbExtraVars = 0;
  MPIIOFunctions::readScalar(fh, nbExtraVars);
  
  getReadData().setNbExtraVars(static_cast<CFuint>(nbExtraVars));
  
  if(nbExtraVars < 0) {
    throw BadFormatException (FromHere(),"Negative number of extra nodal variables in CFmesh");
  }
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbExtraVars() end\n");
}
      
      
//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbExtraNodalVars(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbExtraNodalVars() start\n");

  CFint nbExtraNodalVars = 0;
  MPIIOFunctions::readScalar(fh, nbExtraNodalVars);
  
  getReadData().setNbExtraNodalVars(static_cast<CFuint>(nbExtraNodalVars));
  
  if(nbExtraNodalVars < 0) {
    throw BadFormatException (FromHere(),"Negative number of extra nodal variables in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbExtraNodalVars() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbExtraStateVars(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbExtraStateVars() start\n");
  
  CFint nbExtraStateVars = 0;
  MPIIOFunctions::readScalar(fh, nbExtraStateVars);
  
  getReadData().setNbExtraStateVars(static_cast<CFuint>(nbExtraStateVars));

  if(nbExtraStateVars < 0) {
    throw BadFormatException (FromHere(),"Negative number of extra state variables in CFmesh");
  }

  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbExtraStateVars() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readExtraVarNames(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraVarNames() start\n");
  
  const CFuint nbExtraVars = getReadData().getNbExtraVars();
  vector<std::string> extraVarNames(nbExtraVars);
  
  for (CFuint i = 0; i < nbExtraVars; ++i) {
    extraVarNames[i] = MPIIOFunctions::readAndTrimString(fh);
  }
  
  getReadData().setExtraVarNames(extraVarNames);
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraVarNames() end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readExtraStateVarNames(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraStateVarNames() start\n");
  
  const CFuint nbExtraStateVars = getReadData().getNbExtraStateVars();
  vector<std::string> extraStateVarNames(nbExtraStateVars);
  
  for (CFuint i = 0; i < nbExtraStateVars; ++i) {
    extraStateVarNames[i] = MPIIOFunctions::readAndTrimString(fh);
  }

  getReadData().setExtraStateVarNames(extraStateVarNames);
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraStateVarNames() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readExtraNodalVarNames(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraNodalVarNames() start\n");

  const CFuint nbExtraNodalVars = getReadData().getNbExtraNodalVars();
  vector<std::string> extraNodalVarNames(nbExtraNodalVars);

  for (CFuint i = 0; i < nbExtraNodalVars; ++i) {
    extraNodalVarNames[i] = MPIIOFunctions::readAndTrimString(fh);
  }

  getReadData().setExtraNodalVarNames(extraNodalVarNames);
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraNodalVarNames() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readExtraVarStrides(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraVarStrides() start\n");
  
  const CFuint nbExtraVars = getReadData().getNbExtraVars();
  vector<CFuint> extraVarStrides(nbExtraVars);
  
  MPIIOFunctions::readArray(fh, &extraVarStrides[0], nbExtraVars);
  getReadData().setExtraVarStrides(extraVarStrides);
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraVarStrides() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readExtraStateVarStrides(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraStateVarStrides() start\n");

  const CFuint nbExtraStateVars = getReadData().getNbExtraStateVars();
  vector<CFuint> extraStateVarStrides(nbExtraStateVars);

  MPIIOFunctions::readArray(fh, &extraStateVarStrides[0], nbExtraStateVars);
  getReadData().setExtraStateVarStrides(extraStateVarStrides);
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraStateVarStrides() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readExtraNodalVarStrides(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraNodalVarStrides() start\n");

  const CFuint nbExtraNodalVars = getReadData().getNbExtraNodalVars();
  vector<CFuint> extraNodalVarStrides(nbExtraNodalVars);
  
  MPIIOFunctions::readArray(fh, &extraNodalVarStrides[0], nbExtraNodalVars);
  getReadData().setExtraNodalVarStrides(extraNodalVarStrides);
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraNodalVarStrides() end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readExtraVars(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraVars() start\n");
  
  getReadData().resizeExtraVars();
  const CFuint nbExtraVars = getReadData().getNbExtraVars();
  const vector<CFuint>& extraVarStrides = *getReadData().getExtraVarStrides();
  
  RealVector extraVars;
  if (nbExtraVars > 0) {
    const CFuint sizeExtraVars = std::accumulate(extraVarStrides.begin(),
                                                 extraVarStrides.end(), 0);
    extraVars.resize(sizeExtraVars, 0.);
    getReadData().prepareExtraVars();
    MPIIOFunctions::readArray(fh, &extraVars[0], sizeExtraVars);
    getReadData().setExtraVar(extraVars);
  }
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readExtraNodalVarStrides() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNbGroups(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbGroups() start\n");

  CFuint nbGroups = 0;
  MPIIOFunctions::readScalar(fh, nbGroups);
  
  if (nbGroups < 1) {
    throw BadFormatException (FromHere(),"Number of nbGroups in file must be at least 1");
  }
  
  getReadData().setNbGroups(nbGroups);
 
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNbTRSs() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readGroupName(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGroupName() start\n");
  
  std::string name = MPIIOFunctions::readAndTrimString(fh);
  
  if(getReadData().getNbGroups() == 1) name = "InnerCells";
  CFLogDebugMin( "Creating Group " + name + "\n");
  
  const CFuint idx = m_groups_idxmap.size();
  m_groups_idxmap[name] = idx;
  getReadData().getGroupNames()->push_back(name);

  // set the current TRS being read
  m_curr_group = name;
  // index
  CFLogDebugMin( "Current TRS " << m_curr_group << " has index " << m_groups_idxmap[m_curr_group] << "\n");
    
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGroupName() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readGroupElementNb(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGroupElementNb() start\n");
  
  CFint nbgeo = 0;
  MPIIOFunctions::readScalar(fh, nbgeo);
  cf_assert(nbgeo > 0);
  getReadData().getGroupSizes()->push_back(nbgeo);
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGroupElementNb() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readGroupElementList(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGroupElementList() start\n");

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
  
  // AL: quick and dirty, this need to be reconsidered
  vector<CFuint> buf(totalNbGroupElems);
  MPIIOFunctions::readArray(fh, &buf[0], buf.size());
  
  // for each element of the group,
  // check if the local processor has the element
  CFuint globalElementID;
  for(CFuint iElem = 0 ; iElem < totalNbGroupElems; ++iElem){
    globalElementID = buf[iElem];
    
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
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readGroupElementList() end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readFromFile(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;
  
  char* fileName = const_cast<char*>(filepath.string().c_str());

  // open the file in parallel
  MPI_File_open(m_comm, fileName, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &m_fh); 
  
  bool keepOnReading = true;
  do {
    keepOnReading = readString(&m_fh);
  } while (keepOnReading);
  
  // close the file in parallel
  MPI_File_close(&m_fh);
  
  finish();
}
      
//////////////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readNodeList(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNodeList() start\n");
  
  const CFuint nbLocalNodes = m_localNodeIDs.size() + m_ghostNodeIDs.size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const std::string parNodeVecName = MeshDataStack::getActive()->
    getPrimaryNamespace() +"_nodes";
  
  DataHandle<Node*,GLOBAL> nodes = MeshDataStack::getActive()->getDataStorage()->
    getGlobalData<Node*>(parNodeVecName);
  
  cf_assert(nbLocalNodes > 0);
  nodes.reserve(nbLocalNodes, dim*sizeof(CFreal));
  getReadData().resizeNodes(nbLocalNodes);
  
  sort(m_localNodeIDs.begin(), m_localNodeIDs.end());
  sort(m_ghostNodeIDs.begin(), m_ghostNodeIDs.end());
  
  if (!m_hasPastNodes && getReadData().storePastNodes()){
    throw BadFormatException
      (FromHere(), "ParCFmeshBinaryFileReader => readPastNodes is asked but PastNodes are not present in the CFmesh");
  }
  
  if(m_hasInterNodes && !getReadData().storeInterNodes()){
    CFLog(INFO, "ParCFmeshBinaryFileReader => InterNodes are present in the CFmesh but will not be read\n");
  }
  
  if(!m_hasInterNodes && getReadData().storeInterNodes()){
    throw BadFormatException
      (FromHere(), "ParCFmeshBinaryFileReader => readInterNodes is asked but InterNodes are not present in the CFmesh");
  }
  
  if(m_hasPastNodes && !getReadData().storePastNodes()){
    CFLog(INFO, "ParCFmeshBinaryFileReader => PastNodes are present in the CFmesh but will not be read\n");
  }
  
  const CFuint nbExtraVars = getReadData().getNbExtraNodalVars();
  const vector<CFuint>& nodalExtraVarsStrides = *getReadData().getExtraNodalVarStrides();
  
  RealVector extraVars;
  if (nbExtraVars > 0) {
    const CFuint sizeExtraVars = std::accumulate(nodalExtraVarsStrides.begin(), nodalExtraVarsStrides.end(), 0);
    extraVars.resize(sizeExtraVars, 0.);
  }
  
  getReadData().prepareNodalExtraVars();
  
  // let each rank read a part of the nodes
  // do a MPI_All_to_allv to gather all the local node data
  
  CFuint nodeSize = dim;
  if (m_hasPastNodes)  {nodeSize += dim;}
  if (m_hasInterNodes) {nodeSize += dim;}
  if (nbExtraVars > 0) {nodeSize += extraVars.size();}
  
  // read "\n" character
  char c; MPIIOFunctions::readScalar(fh, c); 
  
  // get the position of the current pointer in the file
  MPI_Offset startListOffset;
  MPI_File_get_position(*fh, &startListOffset);
  
  // set the number of nodes to read in each processor
  vector<CFuint> nbNodesPerProc(m_nbProc);
  // min-max IDs of node data read by each process 
  vector<pair<CFuint, CFuint> > ranges(m_nbProc);
  setReadingRanges(m_totNbNodes, ranges, nbNodesPerProc);
  const CFuint nbNodesUpToRank = accumulate
    (&nbNodesPerProc[0], &nbNodesPerProc[0] + m_myRank, 0);
  
  MPI_Offset startPos = startListOffset + nbNodesUpToRank*nodeSize*sizeof(CFreal);
  const CFuint sizeRead = nbNodesPerProc[m_myRank]*nodeSize;
  vector<CFreal> buf(nbNodesPerProc[0]*nodeSize); // buffer is oversized 
  
  // each processor reads the portion of nodes that is associated to its rank
  CFLog(VERBOSE, "ParCFmeshBinaryFileReader::readNodeList() => nodes read in position [" << startPos << 
	", " << startPos + sizeRead*sizeof(CFreal) << "]\n");
  
  MPIIOFunctions::readAll("ParCFmeshBinaryFileReader::readNodeList()", fh, startPos, &buf[0], (CFuint)sizeRead, m_maxBuffSize);
  
  vector<CFreal> localNodesData(m_localNodeIDs.size()*nodeSize);
  getLocalData(buf, ranges, m_localNodeIDs, nodeSize, localNodesData);
  
  vector<CFreal> ghostNodesData;
  if (m_ghostNodeIDs.size() > 0) {
    ghostNodesData.resize(m_ghostNodeIDs.size()*nodeSize);
    getLocalData(buf, ranges, m_ghostNodeIDs, nodeSize, ghostNodesData);
  }
  
  createNodesAll(localNodesData, ghostNodesData, nodes);
  
  MPI_Barrier(m_comm);
  MPI_Offset endPos = startListOffset + m_totNbNodes*nodeSize*sizeof(CFreal);
  MPI_File_seek(*fh, endPos, MPI_SEEK_SET);
  
  CFLogDebugMin("m_localNodeIDs.size() = " << m_localNodeIDs.size() << "\n");
  CFLogDebugMin("m_ghostNodeIDs.size() = " << m_ghostNodeIDs.size() << "\n");
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readNodeList() end\n");
}
      
//////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::createNodesAll(const vector<CFreal>& localNodesData, 
					       const vector<CFreal>& ghostNodesData, 
					       DataHandle<Node*,GLOBAL> nodes)
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  RealVector tmpNode(0.0, dim);
  RealVector tmpPastNode(0.0, dim);
  RealVector tmpInterNode(0.0, dim);
  const CFuint nbLocalNodes = m_localNodeIDs.size() + m_ghostNodeIDs.size();
  const CFuint nbExtraVars = getReadData().getNbExtraNodalVars();
  const vector<CFuint>& nodalExtraVarsStrides = *getReadData().getExtraNodalVarStrides();
  
  RealVector extraVars;
  if (nbExtraVars > 0) {
    const CFuint sizeExtraVars = std::accumulate(nodalExtraVarsStrides.begin(), 
						 nodalExtraVarsStrides.end(), 0);
    extraVars.resize(sizeExtraVars, 0.);
  }
  
  // create a sorted single list of all global IDs locally present 
  vector<CFuint> globalIDs; globalIDs.reserve(nbLocalNodes);
  for (CFuint i = 0; i < m_localNodeIDs.size(); ++i) {
    globalIDs.push_back(m_localNodeIDs[i]);
  }
  for (CFuint i = 0; i < m_ghostNodeIDs.size(); ++i) {
    globalIDs.push_back(m_ghostNodeIDs[i]);
  }
  sort(globalIDs.begin(), globalIDs.end());
  
  CFuint countBufLocal = 0;
  CFuint countBufGhost = 0;
  CFuint countLocals = 0;
  for (CFuint iNode = 0; iNode < nbLocalNodes; ++iNode) {
    CFuint localID = 0;
    bool isGhost = false;
    bool isFound = false;
    const CFuint globalID = globalIDs[iNode];
    const vector<CFreal>* nodesData = NULL;
    CFuint* countBuf = NULL; 
    if (hasEntry(m_localNodeIDs, globalID)) {
      countLocals++;
      localID = nodes.addLocalPoint (globalID);
      cf_assert(localID < nbLocalNodes);
      isFound = true;
      nodesData = &localNodesData;
      countBuf  = &countBufLocal;
    }
    else if (hasEntry(m_ghostNodeIDs, globalID)) {
      countLocals++;
      localID = nodes.addGhostPoint (globalID);
      cf_assert(localID < nbLocalNodes);
      isGhost = true;
      isFound = true;
      nodesData = &ghostNodesData;
      countBuf  = &countBufGhost;
    }
    cf_assert(isFound);
    
    // read nodes
    for (CFuint d = 0; d < dim; ++d) {tmpNode[d] = (*nodesData)[(*countBuf)++];}
    // read past nodes
    if (m_hasPastNodes) {
      for (CFuint d = 0; d < dim; ++d) {tmpPastNode[d] = (*nodesData)[(*countBuf)++];}
    }
    // read inter nodes
    if (m_hasInterNodes) {
      for (CFuint d = 0; d < dim; ++d) {tmpInterNode[d] = (*nodesData)[(*countBuf)++];}
    }
    // read extra variables
    if (nbExtraVars  > 0) {
      for (CFuint d = 0; d < extraVars.size(); ++d) {
	extraVars[d] = (*nodesData)[(*countBuf)++];
      }
    }
    
    Node* newNode = getReadData().createNode
      (localID, nodes.getGlobalData(localID), tmpNode, !isGhost);
    newNode->setGlobalID(globalID);
    
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
  
  cf_assert(countLocals == nbLocalNodes);
}
      
//////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::createStatesAll(const vector<CFreal>& localStatesData, 
						const vector<CFreal>& ghostStatesData, 
						DataHandle<State*,GLOBAL> states)
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  State tmpState;
  State dummyReadState;
  RealVector tmpPastState(0.0, nbEqs);
  RealVector tmpInterState(0.0, nbEqs);
  RealVector readState(0.0, m_originalNbEqs);
  cf_assert(m_originalNbEqs > 0);
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileReader::createStatesAll() => nbLocalStates = " << m_localStateIDs.size() << "\n");
  CFLog(VERBOSE, "ParCFmeshBinaryFileReader::createStatesAll() => nbGhostStates = " << m_ghostStateIDs.size() << "\n");
  
  const CFuint nbLocalStates = m_localStateIDs.size() + m_ghostStateIDs.size();
  const CFuint nbExtraVars = getReadData().getNbExtraStateVars();
  const vector<CFuint>& stateExtraVarsStrides = *getReadData().getExtraStateVarStrides();
  
  RealVector extraVars;
  if (nbExtraVars > 0) {
    const CFuint sizeExtraVars = std::accumulate(stateExtraVarsStrides.begin(), 
						 stateExtraVarsStrides.end(), 0);
    extraVars.resize(sizeExtraVars, 0.);
  }
  
  // create a sorted single list of all global IDs locally present 
  vector<CFuint> globalIDs; globalIDs.reserve(nbLocalStates);
  for (CFuint i = 0; i < m_localStateIDs.size(); ++i) {
    globalIDs.push_back(m_localStateIDs[i]);
  }
  for (CFuint i = 0; i < m_ghostStateIDs.size(); ++i) {
    globalIDs.push_back(m_ghostStateIDs[i]);
  }
  sort(globalIDs.begin(), globalIDs.end());
  
  bool hasTransformer = false;
  if (m_inputToUpdateVecStr != "Identity") {
    hasTransformer = true;
  }
  
  CFuint countBufLocal = 0;
  CFuint countBufGhost = 0;
  CFuint countLocals = 0;
  for (CFuint iState = 0; iState < nbLocalStates; ++iState) {
    CFuint localID = 0;
    bool isGhost = false;
    bool isFound = false;
    const CFuint globalID = globalIDs[iState];
    const vector<CFreal>* statesData = NULL;
    CFuint* countBuf = NULL; 
    if (hasEntry(m_localStateIDs, globalID)) {
      countLocals++;
      localID = states.addLocalPoint (globalID);
      cf_assert(localID < nbLocalStates);
      isFound = true;
      statesData = &localStatesData;
      countBuf  = &countBufLocal;
    }
    else if (hasEntry(m_ghostStateIDs, globalID)) {
      countLocals++;
      localID = states.addGhostPoint (globalID);
      cf_assert(localID < nbLocalStates);
      isGhost = true;
      isFound = true;
      statesData = &ghostStatesData;
      countBuf  = &countBufGhost;
    }
    cf_assert(isFound);
    
    // no init values were used
    if (m_useInitValues.size() == 0) {
      // read states
      for (CFuint d = 0; d < readState.size(); ++d) {readState[d] = (*statesData)[(*countBuf)++];}
      // read past states
      if (m_hasPastStates) {
	for (CFuint d = 0; d < tmpPastState.size(); ++d) {tmpPastState[d] = (*statesData)[(*countBuf)++];}
      }
      // read inter states
      if (m_hasInterStates) {
	for (CFuint d = 0; d < tmpInterState.size(); ++d) {tmpInterState[d] = (*statesData)[(*countBuf)++];}
      }
      // read extra variables
      if (nbExtraVars  > 0) {
	for (CFuint d = 0; d < extraVars.size(); ++d) {
	  extraVars[d] = (*statesData)[(*countBuf)++];
	}
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
      
      // read states
      for (CFuint d = 0; d < readState.size(); ++d) {
	readState[d] = (*statesData)[(*countBuf)++];
      }
      
      // read past states
      if (m_hasPastStates) {
	for (CFuint d = 0; d < tmpPastState.size(); ++d) {
	  tmpPastState[d] = (*statesData)[(*countBuf)++];
	}
      }
      
      // read inter states
      if (m_hasInterStates) {
	for (CFuint d = 0; d < tmpInterState.size(); ++d) {
	  tmpInterState[d] = (*statesData)[(*countBuf)++];
	}
      }
      
      // read extra variables
      if (nbExtraVars > 0) {
	for (CFuint d = 0; d < extraVars.size(); ++d) {
	  extraVars[d] = (*statesData)[(*countBuf)++];
	}
      }
      
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	if (!m_useInitValues[iEq]) {
	  tmpState[iEq] = readState[iEq];
	}
	else {
	  // here we use all initial values or all initial values IDs
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
    }
    
    State* newState = getReadData().createState
      (localID, states.getGlobalData(localID), tmpState, !isGhost);
    newState->setGlobalID(globalID);
    
    if (m_hasPastStates) {
      getReadData().setPastState(localID, tmpPastState);
    }
    
    if (m_hasInterStates) {
      getReadData().setInterState(localID, tmpInterState);
    }
    
    // set the state extra variable
    if (nbExtraVars > 0) {
      getReadData().setStateExtraVar(localID, extraVars);
    }
    
    // getReadData().setStateLocalToGlobal (localID, globalID);
    // getReadData().setLocalState(localID, !isGhost);
  }
  
  cf_assert(countLocals == nbLocalStates);
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readStateList() end\n");
}
      
//////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::setSendCountDispl(const std::vector<CFuint>& ranks,
						  const CFuint stride,
						  std::vector<int>& sendCount,
						  std::vector<int>& sendDispl)
{
  sendDispl[ranks[0]] = 0;
  CFuint displ = 0;
  for (CFuint i = 0; i < ranks.size(); ++i) {
    sendCount[ranks[i]]++;
    // case with i == 0 is already treated out of the loop
    if (i > 0) {
      if (ranks[i] != ranks[i-1]) {
	sendDispl[ranks[i]] = displ;
      }
    }
    displ += stride;
  }
}
  
//////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::setRecvDispl(const std::vector<int>& recvCount,
					     std::vector<int>& recvDispl)
{
  recvDispl[0] = 0;
  CFuint count = recvCount[0];
  for (CFuint i = 1; i < m_nbProc; ++i) {
    if (recvCount[i] > 0) {
      recvDispl[i] = count;
    }
    count += recvCount[i];
  }
}

//////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::getLocalData(const vector<CFreal>& buf, 
					     const vector<pair<CFuint, CFuint> > ranges,
					     const vector<CFuint>& listIDs, 
					     const CFuint nodeSize, 
					     vector<CFreal>& recvBuf)
{
  // each process requests data from processes who have read nodes whose
  // global IDs are associated to their rank 
  vector<CFuint> donorRanks(listIDs.size());
  vector<CFuint> localIDInDonor(listIDs.size());
  for (CFuint i = 0; i < listIDs.size(); ++i) {
    // rank that has read node data corresponding to locally stored (global) node ID 
    setRankLocalID(ranges, listIDs[i], donorRanks[i], localIDInDonor[i]);
  }
    
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
  
  // two total exchanges for IDs list
  // 1) process I send process J the IDs of the data that I needs from J 
  // 2) process J send to I data as agreed
  
  setSendCountDispl(donorRanks, 1, sendCount, sendDispl);
  
  MPIError::getInstance().check
    ("MPI_Alltoall", "ParCFmeshBinaryFileReader::getLocalData()", 
     MPI_Alltoall(&sendCount[0], 1, MPIStructDef::getMPIType(&sendCount[0]), 
		  &recvCount[0], 1, MPIStructDef::getMPIType(&recvCount[0]), m_comm));
  
  // print some hardcore info in case of need for debugging
  // cout << m_myRank << CFPrintContainer<vector<CFint> >("sendCount  = ", &sendCount);
  // cout << m_myRank << CFPrintContainer<vector<CFint> >("recvCount  = ", &recvCount);
  
  setRecvDispl(recvCount, recvDispl);
  
  // CFLogDebugMin(CFPrintContainer<vector<CFint> >("sendDispl  = ", &sendDispl));
  // CFLogDebugMin(CFPrintContainer<vector<CFint> >("recvDispl  = ", &recvDispl));
  
  const CFuint totRecvCount = std::accumulate(recvCount.begin(), recvCount.end(), 0);
  vector<CFuint> localIDToSend(totRecvCount, 0);
  
  CFLog(DEBUG_MIN, CFPrintContainer<vector<CFuint> >("donorRanks      = ", &donorRanks) << "\n");
  
  // rank X receives from rank Y (in position corresponding to rank X in the recv buffer) 
  // the list of its own locally read nodes (by local within the list of read nodes) to send to Y 
  MPIError::getInstance().check
    ("MPI_Alltoallv", "ParCFmeshBinaryFileReader::getLocalData()", 
     MPI_Alltoallv(&localIDInDonor[0], &sendCount[0], &sendDispl[0], 
		   MPIStructDef::getMPIType(&localIDInDonor[0]), 
		   &localIDToSend[0], &recvCount[0], &recvDispl[0], 
		   MPIStructDef::getMPIType(&localIDToSend[0]), m_comm));
  
  // if (m_myRank == 0) {
  //   ofstream fout("file0");
  //   fout << "P0 RECEIVES IDs: "; 
  //   for (CFuint i = 0; i < m_localNodeIDs.size(); ++i) {
  //     if (donorRanks[i] == 1) {
  // 	fout << m_localNodeIDs[i] << " ";
  // 	counter++;
  //     }
  //   }
  //   fout.close();
  // }
  
  // if (m_myRank == 1) {
  //   ofstream fout("file1");
  //   fout << "P1 SENDS ID: ";
  //   for (CFuint i = 0; i < recvCount[0]; ++i) {
  //     fout << ranges[1].first + localIDToSend[i] << " ";
  //   }
  //   fout.close();
  // }
    
  // CFLog(INFO, CFPrintContainer<vector<CFuint> >("localIDInDonor  = ", &localIDInDonor) << "\n");
  // CFLog(INFO, CFPrintContainer<vector<CFuint> >("localIDToSend   = ", &localIDToSend) << "\n");
  
  for (CFuint i = 0; i < m_nbProc; ++i) {
    sendCount[i] = recvCount[i]*nodeSize;
    sendDispl[i] = recvDispl[i]*nodeSize;
  }
    
  MPIError::getInstance().check
    ("MPI_Alltoall", "ParCFmeshBinaryFileReader::getLocalData()", 
     MPI_Alltoall(&sendCount[0], 1, MPIStructDef::getMPIType(&sendCount[0]), 
		  &recvCount[0], 1, MPIStructDef::getMPIType(&recvCount[0]), m_comm));
    
  setRecvDispl(recvCount, recvDispl);
  
  const CFuint totRecvCountNodes = std::accumulate(recvCount.begin(), recvCount.end(), 0);
  vector<CFreal> sendBuf(localIDToSend.size()*nodeSize);
  cf_assert(totRecvCountNodes == listIDs.size()*nodeSize);
  
  CFuint index = 0;
  for (CFuint i = 0; i < localIDToSend.size(); ++i) {
    const CFuint start = localIDToSend[i]*nodeSize;
    for (CFuint d = 0; d < nodeSize; ++d, ++index) {
      cf_assert(index < sendBuf.size());
      cf_assert(start+d < buf.size());
      sendBuf[index] = buf[start+d];
    }
  }
    
  MPIError::getInstance().check
    ("MPI_Alltoallv", "ParCFmeshBinaryFileReader::getLocalData()", 
     MPI_Alltoallv(&sendBuf[0], &sendCount[0], &sendDispl[0], 
		   MPIStructDef::getMPIType(&sendBuf[0]), 
		   &recvBuf[0], &recvCount[0], &recvDispl[0], 
		   MPIStructDef::getMPIType(&recvBuf[0]), m_comm));
  
  // CFLog(INFO, CFPrintContainer<vector<CFreal> >("sendBuf  = ", &sendBuf, nodeSize) << "\n");
  // CFLog(INFO, CFPrintContainer<vector<CFreal> >("recvBuf  = ", &recvBuf, nodeSize) << "\n");
}
      
//////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::setReadingRanges
(const CFuint total, 
 vector<pair<CFuint, CFuint> >& ranges, 
 vector<CFuint>& nbNodesPerProc)
{
  const CFuint nbReaders = nbNodesPerProc.size();
  const CFuint ne = total/nbReaders;
  nbNodesPerProc[0] = ne + total%nbReaders;
  ranges[0].first  = 0; 
  ranges[0].second = nbNodesPerProc[0]-1;
  for (CFuint i = 1; i < nbReaders; ++i) {
    nbNodesPerProc[i] = ne;
    ranges[i].first  = ranges[i-1].second + 1;
    ranges[i].second = ranges[i].first + ne - 1;
  }
  cf_assert(ranges[nbReaders-1].second == total-1);
}
      
//////////////////////////////////////////////////////////////////////

void ParCFmeshBinaryFileReader::readStateList(MPI_File* fh)
{
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readStateList() start\n");
  
  if ((m_initValues.size() != m_useInitValues.size()) && 
      (m_initValuesIDs.size() != m_useInitValues.size())) {
    CFLog(VERBOSE, "ParCFmeshBinaryFileReader::readStateList() => m_initValues.size()    = " << m_initValues.size() << "\n");
    CFLog(VERBOSE, "ParCFmeshBinaryFileReader::readStateList() => m_initValuesIDs.size() = " << m_initValuesIDs.size() << "\n");
    CFLog(VERBOSE, "ParCFmeshBinaryFileReader::readStateList() => m_useInitValues.size() = " << m_useInitValues.size() << "\n");
    
    throw BadFormatException
      (FromHere(),"ParCFmeshBinaryFileReader => m_initValues && m_initValuesIDs sizes != m_useInitValues.size()");
  }
  
  // read flag telling if there is a solution
  CFuint flag = 0;
  MPIIOFunctions::readScalar(fh, flag);
  bool isWithSolution = (bool) flag;
  getReadData().setWithSolution(isWithSolution);
  
  const CFuint nbLocalStates = m_localStateIDs.size() + m_ghostStateIDs.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const std::string parStateVecName = MeshDataStack::getActive()->
    getPrimaryNamespace() + "_states";
  
  DataHandle<State*,GLOBAL> states = MeshDataStack::getActive()->getDataStorage()->
    getGlobalData<State*>(parStateVecName);
  
  cf_assert(nbLocalStates > 0);
  states.reserve(nbLocalStates, nbEqs*sizeof(CFreal));
  getReadData().resizeStates(nbLocalStates);
  
  sort(m_localStateIDs.begin(), m_localStateIDs.end());
  sort(m_ghostStateIDs.begin(), m_ghostStateIDs.end());
  
  if(!m_hasPastStates && getReadData().storePastStates()){
    throw BadFormatException
      (FromHere(),"ParCFmeshBinaryFileReader => readPastStates is asked but PastStates are not present in the CFmesh");
  }
  
  if(m_hasPastStates && !getReadData().storePastStates()){
    CFLog(INFO, "ParCFmeshBinaryFileReader => PastStates are present in the CFmesh but will not be read\n");
  }
  
  if(!m_hasInterStates && getReadData().storeInterStates()){
    throw BadFormatException
      (FromHere(),"ParCFmeshBinaryFileReader => readInterStates is asked but InterStates are not present in the CFmesh");
  }
  
  if(m_hasInterStates && !getReadData().storeInterStates()){
    CFLog(INFO, "ParCFmeshBinaryFileReader => InterStates are present in the CFmesh but will not be read\n");
  }
  
  const CFuint nbExtraVars = getReadData().getNbExtraStateVars();
  const vector<CFuint>& stateExtraVarsStrides = *getReadData().getExtraStateVarStrides();
  RealVector extraVars;
  if (nbExtraVars > 0) {
    const CFuint sizeExtraVars = std::accumulate(stateExtraVarsStrides.begin(), stateExtraVarsStrides.end(), 0);
    extraVars.resize(sizeExtraVars, 0.);
  }
  
  getReadData().prepareStateExtraVars();
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileReader::readStateList() => nbExtraVars = " << nbExtraVars << "\n");
  
  bool hasTransformer = false;
  
  // read the state
  if (isWithSolution) { // warn if nbeqs differs from original and we dont provide mapping of variable ids
    if (m_originalNbEqs != nbEqs && (m_useInitValues.size() == 0)) {
      CFLog(WARN, "************************************************************\n");
      CFLog(WARN, "WARNING: Nb equations in CFmesh differs from Physical Model.\n");
      CFLog(WARN, "         User did not provide mapping IDs.\n");
      CFLog(WARN, "         May incur in wrong initialization of solution.\n");
      CFLog(WARN, "************************************************************\n");
    }
    
    // configure the variable transformer
    std::string name = MeshDataStack::getActive()->getPrimaryNamespace();
    SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
    SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
    SafePtr<VarSetTransformer::PROVIDER> vecTransProv = CFNULL;
    
    CFLog(VERBOSE, "ParCFmeshBinaryFileReader::configure(): configuring " << m_inputToUpdateVecStr << "\n");
    try {
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
  
  CFuint stateSize = m_originalNbEqs;
  cf_assert(m_originalNbEqs > 0);
  if (m_hasPastStates) {stateSize += nbEqs;}
  if (m_hasInterStates) {stateSize += nbEqs;}
  if (nbExtraVars > 0) {stateSize += extraVars.size();}
  
  CFLog(VERBOSE, "ParCFmeshBinaryFileReader::readStateList() => stateSize = " << stateSize << "\n");
  
  // read "\n" character 
  char c; MPIIOFunctions::readScalar(fh, c);
  
  // get the position of the current pointer in the file
  MPI_Offset startListOffset;
  MPI_File_get_position(*fh, &startListOffset);
  
  vector<CFreal> localStatesData(m_localStateIDs.size()*stateSize);
  vector<CFreal> ghostStatesData;
  if (m_ghostStateIDs.size() > 0) {
    ghostStatesData.resize(m_ghostStateIDs.size()*stateSize);
  }
  
  if (isWithSolution)  { 
    // set the number of states to read in each processor
    vector<CFuint> nbStatesPerProc(m_nbProc);
    // min-max IDs of node data read by each process 
    vector<pair<CFuint, CFuint> > ranges(m_nbProc);
    setReadingRanges(m_totNbStates, ranges, nbStatesPerProc);
    const CFuint nbStatesUpToRank = accumulate
      (&nbStatesPerProc[0], &nbStatesPerProc[0] + m_myRank, 0);
    
    MPI_Offset startPos = startListOffset + nbStatesUpToRank*stateSize*sizeof(CFreal);
    const CFuint sizeRead = nbStatesPerProc[m_myRank]*stateSize;
    vector<CFreal> buf(nbStatesPerProc[0]*stateSize); // buffer is oversized 
    
    // each processor reads the portion of nodes that is associated to its rank
    CFLog(VERBOSE, "ParCFmeshBinaryFileReader::readStateList() => states read in position [" << startPos << 
	  ", " << startPos + sizeRead*sizeof(CFreal) << "]\n");
    
    MPIIOFunctions::readAll("ParCFmeshBinaryFileReader::readStateList()", fh, startPos, &buf[0], (CFuint)sizeRead, m_maxBuffSize);
    getLocalData(buf, ranges, m_localStateIDs, stateSize, localStatesData);
    
    if (m_ghostStateIDs.size() > 0) {
      getLocalData(buf, ranges, m_ghostStateIDs, stateSize, ghostStatesData);
    }  
  }
  
  createStatesAll(localStatesData, ghostStatesData, states);
  
  if (isWithSolution)  { 
    MPI_Barrier(m_comm);
    MPI_Offset endPos = startListOffset + m_totNbStates*stateSize*sizeof(CFreal);
    MPI_File_seek(*fh, endPos, MPI_SEEK_SET);
  }
  
  CFLogDebugMin("m_localStateIDs.size() = " << m_localStateIDs.size() << "\n");
  CFLogDebugMin("m_ghostStateIDs.size() = " << m_ghostStateIDs.size() << "\n");
  
  CFLogDebugMin( "ParCFmeshBinaryFileReader::readStateList() end\n");
}
      
//////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
