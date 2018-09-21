// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include "Common/PE.hh"
#include "Common/MPI/MPIIOFunctions.hh"
#include "Common/CFMap.hh"
#include "Common/OSystem.hh"
#include "Common/BadValueException.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Framework/MapGeoEnt.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataHandleOutput.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/WriteListMap.hh"
#include "Framework/VarSetTransformer.hh"

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/ParWriteSolution.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Common;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

 MethodCommandProvider<ParWriteSolution, TecWriterData, TecplotWriterModule>
 parWriteSolutionProvider("ParWriteSolution");

 //////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("OnlyNodal", "This flag forces output to be all nodal.");
  options.addConfigOption< std::string>("FileFormat","Format to write Tecplot file."); 
  options.addConfigOption< CFuint >("NbWriters", "Number of writers (and MPI groups)");
  options.addConfigOption< CFuint >("NbWritersPerNode", "Number of writers per node");
  options.addConfigOption< int >("MaxBuffSize", "Maximum buffer size for MPI I/O"); 
}
    
//////////////////////////////////////////////////////////////////////////////

 ParWriteSolution::ParWriteSolution(const std::string& name) : 
   TecWriterCom(name),
   ParFileWriter(),
   socket_nodes("nodes"),
   socket_nstatesProxy("nstatesProxy"),
   _mapTrsName2TecplotData(),
   _mapGlobal2LocalNodeID(),
   m_nodalvars(),
   m_ccvars(),
   _intWordFormatSize(),
   _isNewBFile(false)
 {
   addConfigOptionsTo(this);
   
   m_onlyNodal = true;
   setParameter("OnlyNodal",&m_onlyNodal);
   
   _fileFormatStr = "ASCII";
   setParameter("FileFormat",&_fileFormatStr);
   
   _nbWriters = 1;
   setParameter("NbWriters",&_nbWriters);
  
   _nbWritersPerNode = 0;
   setParameter("NbWritersPerNode",&_nbWritersPerNode);
 
   _maxBuffSize = 2147479200; // (CFuint) std::numeric_limits<int>::max();
   setParameter("MaxBuffSize",&_maxBuffSize);
 }
    
 //////////////////////////////////////////////////////////////////////////////

 ParWriteSolution::~ParWriteSolution() 
 {
   CFAUTOTRACE;
   
   for (CFuint i = 0; i < _mapTrsName2TecplotData.size(); ++i) {
    delete _mapTrsName2TecplotData[i];
   }
 }

 //////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::execute()
{
  CFAUTOTRACE;
 
  CFLog(VERBOSE, "ParWriteSolution::execute() => start\n");
  
  const TecplotTRSType& tt = *_mapTrsName2TecplotData.find("InnerCells");
  if (hasChangedMesh(tt)) {
    buildNodeIDMapping();
  }
    
  if (_fileFormatStr == "ASCII") {    
    const boost::filesystem::path cfgpath = getMethodData().getFilename();
    if (!getMethodData().onlySurface()) {
      // write inner domain data
      writeData(cfgpath, _isNewFile, string("Unstructured grid data"), &ParWriteSolution::writeInnerData);
    }
    
    if (!getMethodData().getSurfaceTRSsToWrite().empty()) {
      // write boundary surface data
      boost::filesystem::path bpath = cfgpath.branch_path() / ( basename(cfgpath) + ".surf" + extension(cfgpath) );
      writeData(bpath, _isNewBFile, string("Boundary data"), &ParWriteSolution::writeBoundaryData);
    }
  }
  
  CFLog(VERBOSE, "ParWriteSolution::execute() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeData(const boost::filesystem::path& filepath, 
				 bool& flag, 
				 const std::string title, 
				 WriterFun fun)
{
  CFLog(VERBOSE, "ParWriteSolution::writeData() [" << title << "] => start\n");
  
  // reset to 0 the new file flag
  flag = false;
  ofstream* file = CFNULL;
  Common::SelfRegistPtr<Environment::FileHandlerOutput>* fhandle =
    Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();
  
  if (_isWriterRank) { 
    // if the file has already been processed once, open in I/O mode
    if (_fileList.count(filepath) > 0) {
      file = &(*fhandle)->open(filepath, ios_base::in | ios_base::out);
    }
    else {
      file = &(*fhandle)->open(filepath);
      
      // if the file is a new one add it to the file list
      _fileList.insert(filepath);
      flag = true;
    }
  }
  
  int isNewFile = (int)flag;
  MPI_Bcast(&isNewFile, 1, MPIStructDef::getMPIType(&isNewFile), _ioRank, _comm);
  
  flag = (isNewFile == 0) ? false : true;
    
  (this->*fun)(filepath, flag, title, file); 
  
  if (_isWriterRank) {
    const string writerName = getMethodData().getNamespace() + "_Writers";
    Group& wg = PE::GetPE().getGroup(writerName);
    MPI_Barrier(wg.comm);
    (*fhandle)->close();
  }

  delete fhandle;

  CFLog(VERBOSE, "ParWriteSolution::writeData() [" << title << "] => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

const std::string ParWriteSolution::getWriterName() const
{
  return getName();
}

//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeToBinaryFile()
{
 CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////
      
void ParWriteSolution::writeInnerData
(const boost::filesystem::path& filepath,
 const bool isNewFile,
 const std::string title,
 std::ofstream* fout)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::writeInnerData() [" << title << "] => start\n");
  CFLog(INFO, "Writing solution to " << filepath.string() << "\n");
  
  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
  datahandle_output->getDataHandles();
  
  // write the TECPLOT file header
  if (isNewFile) {
    writeHeader(fout, title, datahandle_output);
  }
  
  std::vector<SafePtr<TopologicalRegionSet> > trsList =
    MeshDataStack::getActive()->getTrsList();
  
  for(CFuint iTrs= 0; iTrs < trsList.size(); ++iTrs) {
    SafePtr<TopologicalRegionSet> trs = trsList[iTrs];
    
    if ((trs->hasTag("inner")) && (trs->hasTag("cell"))) {
      SafePtr<vector<ElementTypeData> > elementType =
	MeshDataStack::getActive()->getElementTypeData(trs->getName());
      TecplotTRSType& tt = *_mapTrsName2TecplotData.find(trs->getName());
      
      // loop over the element types
      // and create a zone in the tecplot file
      // for each element type
      for (CFuint iType = 0; iType < elementType->size(); ++iType) {
	ElementTypeData& eType = (*elementType)[iType];
	
	const bool meshChanged = hasChangedMesh(iType, tt);
	if (isNewFile || meshChanged) {
	  vector<MPI_Offset>& headerOffset = tt.headerOffset[iType];
	  
	  // print zone header: one zone per element type
	  if (_myRank == _ioRank) {
	    headerOffset[0] = fout->tellp();
	    // write the header for the inner domain
	    writeInnerZoneHeader(fout, iType, eType, trs);
	    headerOffset[1] = fout->tellp();
	  }
	  
	  // communicate the current file position to all processor
	  MPI_Bcast(&headerOffset[0], 2, MPIStructDef::getMPIOffsetType(), _ioRank, _comm);
	  
	  // backup the total counts for this element type
	  tt.oldNbNodesElemsInType[iType].first  = tt.totalNbNodesInType[iType];
	  tt.oldNbNodesElemsInType[iType].second = eType.getNbTotalElems();
	}
	
	// print nodal coordinates and stored node variables
	writeNodeList(fout, iType, trs, false);
	
	if (isNewFile || meshChanged) {
	  // write element-node connectivity
	  writeElementList(fout, iType, eType.getNbNodes(),
			   eType.getNbTotalElems(), eType.getNbElems(),  
			   eType.getStartIdx(), eType.getGeoOrder(), trs);
	}
      }
    } //end if inner cells
  } //end loop over trs
  
  CFLog(VERBOSE, "ParWriteSolution::writeInnerData() [" << title << "] => end\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::setup()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::setup() => start\n");
  
  TecWriterCom::setup();
  ParFileWriter::setWriterGroup();
  
  // store the names of additional variables 
  m_ccvars.clear();
  m_nodalvars.clear();
  
  SafePtr<DataHandleOutput> dh = getMethodData().getDataHOutput();
  dh->getDataHandles();
  
  if (dh.isNotNull()) {
    vector<string> dh_varnames = dh->getVarNames();
    for (CFuint i = 0; i < dh_varnames.size() ; ++i) {
      m_nodalvars.push_back(dh_varnames[i]);
    }
    
    vector<string> dh_ccvarnames = dh->getCCVarNames();
    for (CFuint i = 0; i < dh_ccvarnames.size() ; ++i) {
      m_ccvars.push_back(dh_ccvarnames[i]);
    }
  }  
  
  SafePtr< vector<ElementTypeData> > me = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbElementTypes = me->size();
  cf_assert(nbElementTypes > 0);
  
  // create TecplotTRSType for the domain interior
  TecplotTRSType* tt = new TecplotTRSType
    (nbElementTypes, MeshDataStack::getActive()->getTrs("InnerCells"));
  _mapTrsName2TecplotData.insert("InnerCells", tt);
  
  // create TecplotTRSType for the domain boundaries
  const vector<string>& surfTRS = getMethodData().getSurfaceTRSsToWrite();
  for(CFuint iTrs = 0; iTrs < surfTRS.size(); ++iTrs) {
    SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(surfTRS[iTrs]);
    TecplotTRSType* tb = new TecplotTRSType(trs->getNbTRs(), trs);
    _mapTrsName2TecplotData.insert(surfTRS[iTrs], tb);
  }
  _mapTrsName2TecplotData.sortKeys();
  
  _mapGlobal2LocalNodeID.reserve(socket_nodes.getDataHandle().size());
    
  CFLog(VERBOSE, "ParWriteSolution::setup() => end\n");
}      
      
//////////////////////////////////////////////////////////////////////////////
 
void ParWriteSolution::unsetup()
{
  CFAUTOTRACE;
  
  // cleanupNodeIDMapping();
  // _mapGlobal2LocalNodeID.clear();
  
  TecWriterCom::unsetup(); 
}
    
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution:: cleanupNodeIDMapping()
{  
  CFAUTOTRACE;
  
  if (_mapTrsName2TecplotData.size() > 0) {
    for (CFuint i = 0; i < _mapTrsName2TecplotData.size(); ++i) {
      _mapTrsName2TecplotData[i]->cleanupMappings();
    }
    //_mapTrsName2TecplotData.clear();
  }
  
  _mapGlobal2LocalNodeID.clear();
  cf_assert(_mapGlobal2LocalNodeID.size() == 0);
}  
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution:: buildNodeIDMapping()
{   
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::buildNodeIDMapping() => start\n");
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  
  const CFuint nbElementTypes = MeshDataStack::getActive()->getElementTypeData()->size();
  cf_assert(nbElementTypes > 0);
  
  cleanupNodeIDMapping();
  
  vector<SafePtr<TopologicalRegionSet> > trsList = 
    MeshDataStack::getActive()->getTrsList();
  
  const CFuint totalNodeCount = MeshDataStack::getActive()->getTotalNodeCount();
  cf_assert(totalNodeCount >= nodes.size());
  
  // compute the integer format size
  const CFuint maxInteger = 100*totalNodeCount;
  _intWordFormatSize = (CFuint)std::log10((CFreal)maxInteger);
  
  CFLog(VERBOSE, "ParWriteSolution::buildNodeIDMapping() => integer format size = " << _intWordFormatSize << "\n");
  
  vector<bool> foundGlobalID(totalNodeCount);
  
  std::vector<std::string> sortedSurfTRS = 
    getMethodData().getSurfaceTRSsToWrite();
  
  sort(sortedSurfTRS.begin(), sortedSurfTRS.end(), less<string>());
  
  for(CFuint iTrs= 0; iTrs < trsList.size(); ++iTrs) {
    SafePtr<TopologicalRegionSet> trs = trsList[iTrs];
    
    if ((trs->hasTag("inner")) && (trs->hasTag("cell"))) {
      SafePtr<vector<ElementTypeData> > elementType =
	MeshDataStack::getActive()->getElementTypeData(trs->getName());
      cf_assert(elementType->size() == nbElementTypes);
      
      CFLog(VERBOSE, "ParWriteSolution::buildNodeIDMapping() for " << trs->getName() << " \n");
      
      const CFuint nbElemTypes = elementType->size();
      for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
	ElementTypeData& eType = (*elementType)[iType];
	
	storeMappings(trs, iType, eType.getNbElems(), eType.getNbNodes(), 
		      eType.getStartIdx(), eType.getEndIdx(), foundGlobalID);
      }
    }
    
    if (binary_search(sortedSurfTRS.begin(), sortedSurfTRS.end(),trs->getName())) {
      CFLog(VERBOSE, "ParWriteSolution::buildNodeIDMapping() for " << trs->getName() << " \n");
      
      const CFuint nbTRs = trs->getNbTRs();
      for (CFuint iType = 0; iType < nbTRs; ++iType) {
	SafePtr<TopologicalRegion> tr = (*trs)[iType];
	const CFuint nbTRGeos = tr->getLocalNbGeoEnts();
	
	// compute maximum number of nodes in TR geometric entities
	CFuint maxNbNodesInTRGeo = 0;
	for (CFuint iGeo = 0; iGeo < nbTRGeos; ++iGeo) {
	  maxNbNodesInTRGeo = max(maxNbNodesInTRGeo, tr->getNbNodesInGeo(iGeo));
	}
	
	CFLog(VERBOSE, "ParWriteSolution::buildNodeIDMapping() => maxNbNodesInTRGeo = " 
	      << maxNbNodesInTRGeo << " \n");
	
	storeMappings(trs, iType, nbTRGeos, maxNbNodesInTRGeo, 
		      0, nbTRGeos, foundGlobalID);
      }
    }
  }
  
  // build the inverse mapping global IDs to local IDs for all nodes 
  for (CFuint i = 0; i < nodes.size(); ++i) {
    _mapGlobal2LocalNodeID.insert(nodes[i]->getGlobalID(), nodes[i]->getLocalID());
  }
  _mapGlobal2LocalNodeID.sortKeys();
  
  CFLog(VERBOSE, "ParWriteSolution::buildNodeIDMapping() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeBoundaryData(const boost::filesystem::path& filepath,
					 const bool isNewFile,
					 const std::string title,
					 std::ofstream* fout)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::writeBoundaryData() => start\n");
  
  CFLog(INFO, "Writing solution to " << filepath.string() << "\n");
    
  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy = socket_nstatesProxy.getDataHandle();
  
  // AL: this is a handle for the nodal states which can be
  // stored as arrays of State*, RealVector* or RealVector
  // but they are always used as arrays of RealVector*
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];
  
  // we assume that the element-node == element-state connectivity
  // @todo tecplot writer should not assume element-to-node
  // connectivity to be the same as element-to-state
  
  /// return if the user has not indicated any TRS
  if (getMethodData().getSurfaceTRSsToWrite().empty()) return;
  
  // write the TECPLOT file header
  if (isNewFile) {
    writeHeader(fout, title, CFNULL);
  }
  
  const std::vector<std::string>& surfTRS = getMethodData().getSurfaceTRSsToWrite();
  
  CFuint countTRToWrite = 0;  
  for(vector<string>::const_iterator itr = surfTRS.begin(); itr != surfTRS.end(); ++itr) {
    Common::SafePtr<TopologicalRegionSet> currTrs = MeshDataStack::getActive()->getTrs(*itr);
    const CFuint nbTRs = currTrs->getNbTRs();
    for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
      SafePtr<TopologicalRegion> tr = currTrs->getTopologicalRegion(iTR);
      countTRToWrite++;
    }
  }
  cf_assert(countTRToWrite > 0);
   
  const vector<vector<CFuint> >&  trsInfo =
    MeshDataStack::getActive()->getTotalTRSInfo();
   
  SafePtr<vector<vector<vector<CFuint> > > > globalGeoIDS =
    MeshDataStack::getActive()->getGlobalTRSGeoIDs();
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  // loop over the TRSs selected by the user
  for(vector<string>::const_iterator itr = surfTRS.begin(); itr != surfTRS.end(); ++itr) {
    SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs(*itr);
    const int iTRS = getGlobalTRSID(trs->getName());
    cf_assert(iTRS >= 0);
    
    TecplotTRSType& tt = *_mapTrsName2TecplotData.find(trs->getName());
    
    const CFuint nbTRs = trs->getNbTRs(); 
    cf_assert(nbTRs > 0);
    cf_assert(nbTRs == (*globalGeoIDS)[iTRS].size());
    
    for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
      SafePtr<TopologicalRegion> tr = trs->getTopologicalRegion(iTR);
      
      // count the maximum number of nodes in a boundary face belonging to this TR
      CFuint maxNbNodesInTRGeo = 0;
      const CFuint nbTRGeos = (*trs)[iTR]->getLocalNbGeoEnts();
      for (CFuint iGeo = 0; iGeo < nbTRGeos; ++iGeo) {
	maxNbNodesInTRGeo = max(maxNbNodesInTRGeo, (*trs)[iTR]->getNbNodesInGeo(iGeo));
      }
      
      CFuint maxNbNodesInGeo = 0;
      MPI_Allreduce(&maxNbNodesInTRGeo, &maxNbNodesInGeo, 1, MPIStructDef::getMPIType
		    (&maxNbNodesInGeo), MPI_MAX, _comm);
      cf_assert(maxNbNodesInGeo > 0);
      
      std::string elemShape = "LINESEG";
      if (dim == 3) {
	// AL: the maximum number of nodes in face is considered: an extra virtual 
	// node is added if TETRA and QUADRILATERAL are both present
	elemShape = (maxNbNodesInGeo == 3) ? "TRIANGLE" : "QUADRILATERAL";
      }
      
      const CFuint totalNbTRGeos = trsInfo[iTRS][iTR];
      
      // print zone header: one zone per TR
      const bool meshChanged = hasChangedMesh(iTR, tt);
      if (isNewFile || meshChanged) {
	vector<MPI_Offset>& headerOffset = tt.headerOffset[iTR];
	if (_myRank == _ioRank) {
	  headerOffset[0] = fout->tellp();
	  
	  writeZoneHeader(fout, iTR, trs->getName(), tt.totalNbNodesInType[iTR], 
			  totalNbTRGeos, elemShape, "\n", true);
	  
	  headerOffset[1] = fout->tellp();
	}
	
	// communicate the current file position to all processor
	MPI_Bcast(&headerOffset[0], 2, MPIStructDef::getMPIOffsetType(), _ioRank, _comm);
	
	// backup the total counts for this element type
	tt.oldNbNodesElemsInType[iTR].first  = tt.totalNbNodesInType[iTR];
	tt.oldNbNodesElemsInType[iTR].second = totalNbTRGeos;
      }
      
      // print nodal coordinates and stored node variables
      writeNodeList(fout, iTR, trs, true);
      
      if (isNewFile || meshChanged) {
	// write element-node connectivity
	writeElementList(fout, iTR, maxNbNodesInGeo, 
			 totalNbTRGeos, nbTRGeos, 0, 
			 CFPolyOrder::ORDER1, trs);
      }
    }
  }
  
  CFLog(VERBOSE, "ParWriteSolution::writeBoundaryData() => end\n");
}
    
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ParWriteSolution::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_nstatesProxy);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
      
// void ParWriteSolution::writeHeader(MPI_File* fh)
// {
//   CFAUTOTRACE;
  
//   CFLog(VERBOSE, "ParWriteSolution::writeHeader() start\n");
  
//   if (_myRank == _ioRank) {
    
//     //  Tecplot Header
//     // maximum string size to write includes 30 characters 
//     MPIIOFunctions::writeKeyValue<char>(fh, "TITLE = ");
//     MPIIOFunctions::writeKeyValue<char>(fh, "Unstructured grid data");
//     MPIIOFunctions::writeKeyValue<char>(fh, "\nVARIABLES = ");
    
//     const CFuint dim = PhysicalModelStack::getActive()->getDim();
//     const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
//     for (CFuint i = 0; i < dim; ++i) {
//       MPIIOFunctions::writeKeyValue<char>(fh, " \"x" + StringOps::to_str(i) + '\"');
//     }
    
//     if (!getMethodData().onlyCoordinates()) {
//       SafePtr<ConvectiveVarSet> outputVarSet = getMethodData().getOutputVarSet();
//       const vector<std::string>& varNames = outputVarSet->getVarNames();
//       cf_assert(varNames.size() == nbEqs);
      
//       for (CFuint i = 0 ;  i < nbEqs; ++i)  {
// 	std::string n = varNames[i];
// 	if ( *n.begin()  != '\"' )  n  = '\"' + n;
// 	if ( *n.rbegin() != '\"' )  n += '\"';
// 	MPIIOFunctions::writeKeyValue<char>(fh, " " + n);
//       }
      
//       if (getMethodData().shouldPrintExtraValues()) {
// 	vector<std::string> extraVarNames = outputVarSet->getExtraVarNames();
// 	for (CFuint i = 0 ;  i < extraVarNames.size(); ++i) {
// 	  MPIIOFunctions::writeKeyValue<char>(fh, " " +  extraVarNames[i]);
// 	}
//       }
      
//       SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
//       vector<string> dhVarNames = datahandle_output->getVarNames();
//       for (CFuint i = 0 ;  i < dhVarNames.size(); ++i) {
// 	MPIIOFunctions::writeKeyValue<char>(fh, " " + dhVarNames[i]);
//       }
//     }
//   }
  
//   CFLog(VERBOSE, "ParWriteSolution::writeHeader() end\n");
// }
  
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeHeader(std::ofstream* fout, 
				   const std::string title,
				   SafePtr<DataHandleOutput> dh)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::writeHeader() \"" << title << "\" start\n");
  
  if (_myRank == _ioRank) {
    //  Tecplot Header
    // maximum string size to write includes 30 characters 
    *fout << "TITLE = " << title;
    *fout << "\nVARIABLES = ";
    
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    for (CFuint i = 0; i < dim; ++i) {
      *fout << " \"x" << StringOps::to_str(i) << '\"';
    }
    
    if (!getMethodData().onlyCoordinates()) {
      SafePtr<ConvectiveVarSet> outputVarSet = getMethodData().getOutputVarSet();
      const vector<std::string>& varNames = outputVarSet->getVarNames();
      cf_assert(varNames.size() == nbEqs);
      
      if (getMethodData().withEquations()) {
	for (CFuint i = 0 ;  i < nbEqs; ++i)  {
	  std::string n = varNames[i];
	  if ( *n.begin()  != '\"' )  n  = '\"' + n;
	  if ( *n.rbegin() != '\"' )  n += '\"';
	  *fout << " " << n;
	}
      }
      
      if (getMethodData().shouldPrintExtraValues()) {
	vector<string> extraVarNames = outputVarSet->getExtraVarNames();
	for (CFuint i = 0 ;  i < extraVarNames.size(); ++i) {
	  *fout << " \"" << extraVarNames[i] << "\"";
	}
      }
      
      if (dh.isNotNull()) {
	vector<string> dhVarNames = dh->getVarNames();
	for (CFuint i = 0 ;  i < dhVarNames.size(); ++i) {
	  *fout << " \"" << dhVarNames[i] << "\"";
	}
	
	// discard DataHandle CC if only NODAL
	if (!m_onlyNodal) {
	  vector<string> dhCCVarNames = dh->getCCVarNames();
	  for (CFuint i = 0; i < dhCCVarNames.size() ; ++i) {
	    *fout << " \"" << dhCCVarNames[i] << "\"";
	  }
	}
      }
    }
    *fout << "\n";
  }
  
  if (_isWriterRank) {
    const string writerName = getMethodData().getNamespace() + "_Writers";
    Group& wg = PE::GetPE().getGroup(writerName);
    MPI_Barrier(wg.comm);
  }
  
  CFLog(VERBOSE, "ParWriteSolutionBlock::writeHeader() \"" << title << "\" end\n");
}
 
//////////////////////////////////////////////////////////////////////////////
  
void ParWriteSolution::writeNodeList(ofstream* fout, const CFuint iType,
				     SafePtr<TopologicalRegionSet> elements,
				     const bool isBoundary)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::writeNodeList() => start\n");
  
  const CFuint dim  = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  SafePtr<ConvectiveVarSet> outputVarSet = getMethodData().getOutputVarSet();
  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  SafePtr<VarSetTransformer> updateToOutput = getMethodData().getUpdateToOutputVarSetTransformer();
  vector<string> extraVarNames = outputVarSet->getExtraVarNames();
  const CFuint nbExtraVars = extraVarNames.size();
  
  CFuint nodesStride = dim; 
  if (!getMethodData().onlyCoordinates()) {
    if (getMethodData().withEquations()) {
      nodesStride += nbEqs;
    }
    
    if (getMethodData().shouldPrintExtraValues()) {
      nodesStride += nbExtraVars;
    }
    // nodal data handle variables 
    nodesStride += datahandle_output->getVarNames().size();
    cf_assert(datahandle_output->getVarNames().size() == m_nodalvars.size());    
  }
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  cf_assert(nodes.size() > 0);
  
  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
    socket_nstatesProxy.getDataHandle();

  // this is a sort of handle for the nodal states
  // (which can be stored as arrays of State*, RealVector* or
  // RealVector but they are used as arrays of RealVector*)
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];
  
  TecplotTRSType& tt = *_mapTrsName2TecplotData.find(elements->getName());
  
  const CFuint nSend = _nbWriters;
  const CFuint nbElementTypes = 1; // just nodes
  const CFuint nbLocalElements = tt.nodesInType[iType].size();
  const CFuint totNbNodes = tt.totalNbNodesInType[iType];
  // size should be set in the VarSet but for safety we resize it here too
  RealVector extraValues;
  if (nbExtraVars > 0) extraValues.resize(nbExtraVars);
  
  
  const string writerName = getMethodData().getNamespace() + "_Writers";
  Group& wg = PE::GetPE().getGroup(writerName);
  
  WriteListMap elementList;
  elementList.reserve(nbElementTypes, nSend, nbLocalElements);
  
  // fill in the writer list
  CFuint totalToSend = 0;
  elementList.fill(totNbNodes, nodesStride, totalToSend);
    
  const CFuint wordFormatSize = 22;
  
  vector<MPI_Offset> wOffset(_nbWriters, tt.headerOffset[iType][1]); 
  // set the offsets for the nodes of this element type
  tt.nodesOffset[iType].first  = tt.headerOffset[iType][1];
  tt.nodesOffset[iType].second = tt.nodesOffset[iType].first + totalToSend*wordFormatSize;
  
  CFLog(VERBOSE, "ParWriteSolution::writeNodeList() => offsets = [" 
	<<  tt.nodesOffset[iType].first << ", " << tt.nodesOffset[iType].second << "]\n");
  
  // update the maximum possible element-list size to send
  const CFuint maxElemSendSize = elementList.getMaxElemSize();
  
  // insert in the write list the local IDs of the elements
  // the range ID is automatically determined inside the WriteListMap
  const vector<CFuint>& nodesInType = tt.nodesInType[iType];
  for (CFuint iElem = 0; iElem < nbLocalElements; ++iElem) {
    const CFuint globalTypeID = tt.mapNodeID2NodeIDByEType[iType]->find(nodesInType[iElem]);
    elementList.insertElemLocalID(iElem, globalTypeID, 0);
  }
  elementList.endElemInsertion(_myRank);
  
  // buffer data to send
  vector<CFreal> sendElements(maxElemSendSize, 0);
  vector<CFreal> elementToPrint(maxElemSendSize, 0);
  
  CFint wRank = -1; 
  CFuint wSendSize = 0;
  CFuint rangeID = 0;
  CFuint countElem = 0;
  for (CFuint is = 0; is < nSend; ++is, ++rangeID) {
    bool isRangeFound = false;
    WriteListMap::List elist = elementList.find(rangeID, isRangeFound);
    
    if (isRangeFound) {
      CFuint eSize = 0;
      for (WriteListMap::ListIterator it = elist.first; it != elist.second; ++it, ++eSize) {
	const CFuint localElemID = it->second;
	const CFuint globalElemID = tt.mapNodeID2NodeIDByEType[iType]->find(nodesInType[localElemID]);
	const CFuint sendElemID = globalElemID - countElem;
	const CFuint nodeID = _mapGlobal2LocalNodeID.find(nodesInType[localElemID]);
	
	// this fix has to be added EVERYWHERE when writing states in parallel
	if (nodes[nodeID]->isParUpdatable()) {
	  CFuint isend = sendElemID*nodesStride;
	  for (CFuint in = 0; in < dim; ++in, ++isend) {
	    cf_assert(isend < sendElements.size());
	    cf_assert(nodeID < nodes.size());
	    sendElements[isend] = (*nodes[nodeID])[in]*refL;
	  }
	  
	  if (!getMethodData().onlyCoordinates()) {
            // get state for updateVarSet
	    const RealVector& currState = *nodalStates.getState(nodeID);
	    const CFuint stateID = nodalStates.getStateLocalID(nodeID);

      State tempState;

	    tempState.setLocalID(stateID);
	    // the node is set  in the temporary state
	    tempState.setSpaceCoordinates(nodes[nodeID]);
	    for (CFuint ieq = 0; ieq < nbEqs; ++ieq) {
	      tempState[ieq] = currState[ieq];
	    }
      State* updateDimState = new State(RealVector(nbEqs));
      State* outputDimState;

	    if (getMethodData().shouldPrintExtraValues()) {
	      // dimensionalize the solution
              updateVarSet->setDimensionalValuesPlusExtraValues
                (tempState, *updateDimState, extraValues);
              outputDimState = updateToOutput->transform(updateDimState);

	      if (getMethodData().withEquations()) {
                for (CFuint in = 0; in < outputDimState->size(); ++in, ++isend) {
		  cf_assert(isend < sendElements.size());
                  sendElements[isend] = (*outputDimState)[in];
		}
	      }
	      for (CFuint in = 0; in < extraValues.size(); ++in, ++isend) {
		cf_assert(isend < sendElements.size());
		sendElements[isend] = extraValues[in];
	      }
	    }
	    else {
	      if (getMethodData().withEquations()) {
                updateVarSet->setDimensionalValues(tempState, *updateDimState);
                outputDimState = updateToOutput->transform(updateDimState);
                for (CFuint in = 0; in < outputDimState->size(); ++in, ++isend) {
		  cf_assert(isend < sendElements.size());
                  sendElements[isend] = (*outputDimState)[in];
		}
	      }
	    }	    
	    delete updateDimState; 
	    datahandle_output->fillStateData(&sendElements[0], stateID, isend);
	  }
	}
      }
      
      cf_assert(eSize*nodesStride <= elementList.getSendDataSize(rangeID));
    }
    
    CFLogDebugMax(_myRank << CFPrintContainer<vector<CFreal> >
		  (" sendElements  = ", &sendElements, nodesStride) << "\n");
    
    const CFuint sendSize = elementList.getSendDataSize(rangeID);
    cf_assert(sendSize <= sendElements.size());
    cf_assert(sendSize <= elementToPrint.size());
    
    // if the rank corresponds to a writing process, record the size to send for this writer
    if (_isWriterRank && wg.globalRanks[is] == _myRank) {
      wSendSize = sendSize; // this should be the total sendsize in the range
      wRank = is;
    }
    
    MPI_Op myMpiOp;
    MPI_Op_create((MPI_User_function *)cmpAndTakeMaxAbs3, 1, &myMpiOp);
    MPI_Reduce(&sendElements[0], &elementToPrint[0], sendSize,
	       MPIStructDef::getMPIType(&sendElements[0]), myMpiOp, wg.globalRanks[is], _comm);
    
    CFLogDebugMax(_myRank << CFPrintContainer<vector<CFreal> >
		  (" elementToPrint  = ", &elementToPrint, nodesStride) << "\n");
    
    if (is == 0) {
      CFLog(VERBOSE, "[0] => wOffset = " << wOffset[0] << ", sendSize = " << sendSize << "\n");
    }
    
    // the offsets for all writers with send ID > current must be incremented  
    for (CFuint iw = is+1; iw < wOffset.size(); ++iw) {
      cf_assert(sendSize > 0);
      wOffset[iw] += sendSize*wordFormatSize;
      CFLog(VERBOSE, "[" << is << ", " << iw << "] => wOffset = " << wOffset[iw] << ", sendSize = " << sendSize << "\n");
    }
    
    // reset the all sendElement list to 0
    for (CFuint i = 0; i < maxElemSendSize; ++i) {
      sendElements[i] = 0;
    }
    
    // update the count element for the current element type
    countElem += elementList.getSendDataSize(rangeID)/nodesStride;
  }
  
  if (_isWriterRank) { 
    CFLog(DEBUG_MIN, "ParWriteSolution::writeNodeList() => P[" << _myRank << "] => offset = " << wOffset[wRank] << "\n");
    // each writer can now concurrently write all the collected data (related to one element type)
    cf_assert(wRank >= 0);
    
    // point to the corresponding writing location
    fout->seekp(wOffset[wRank]);
    
    CFLog(VERBOSE, "wSendSize = " << wSendSize << ", nodesStride = " << nodesStride << "\n");
    for (CFuint i = 0; i < wSendSize; ++i) {
      // this fix is needed for ensuring consistency in the format and avoid 
      // to have entries with 3 digits in the exponent (e.g. +A.BC..e-XXX)
      if (std::abs(elementToPrint[i]) < 1e-50) elementToPrint[i] = 0.;
      
      // const CFreal coeff = (countN < dim) ? refL : 1.;
      // this format corresponds to a line of 22*nodesStride bytes  
      fout->precision(14);
      fout->setf(ios::scientific, ios::floatfield); 
      fout->setf(ios::showpos);
      if ((i+1)%nodesStride > 0) {
	*fout << elementToPrint[i] << " ";
      }
      else {
	*fout << elementToPrint[i] << "\n";
      }
    }
    
    MPI_Offset lastpos = fout->tellp();
    MPI_Offset maxpos  = 0;
    MPI_Allreduce(&lastpos, &maxpos, 1, MPIStructDef::getMPIOffsetType(), MPI_MAX, wg.comm);
    
    if (tt.nodesOffset[iType].second != maxpos) {
      CFLog(WARN, "tt.nodesOffset[iType].second (" << tt.nodesOffset[iType].second  
	    << ") != maxpos (" << maxpos << ")\n");
      cf_assert(tt.nodesOffset[iType].second == maxpos);
    }
    
    fout->seekp(maxpos);
  }
  
  //reset the all sendElement list to 0
  for (CFuint i = 0; i < maxElemSendSize; ++i) {
    elementToPrint[i] = 0;
  }
  
  CFLog(VERBOSE, "ParWriteSolution::writeNodeList() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeElementList
(std::ofstream* fout,
 const CFuint iType,
 const CFuint nbNodesInType,
 const CFuint nbElementsInType,
 const CFuint nbLocalElements,
 const CFuint startElementID,
 const CFuint geoOrder,
 SafePtr<TopologicalRegionSet> elements)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ParWriteSolution::writeElementList() => start\n");
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  
  const CFuint nSend = _nbWriters;
  const CFuint nbElementTypes = 1; // one element type at a time
  
  TecplotTRSType& tt = *_mapTrsName2TecplotData.find(elements->getName());
  const bool isCell = elements->hasTag("cell");
  
  WriteListMap elementList;
  elementList.reserve(nbElementTypes, nSend, nbLocalElements);
  
  // store global info about the global ID ranges for sending
  // and the size of each send
  // each send will involve one writing process which wil collect all data 
  const CFuint totalSize = nbElementsInType*nbNodesInType;
  
  // fill in the writer list
  CFuint totalToSend = 0;
  elementList.fill(nbElementsInType, nbNodesInType, totalToSend);
  cf_assert(totalToSend == totalSize);
  
  // update the maximum possible element-list size to send
  const CFuint maxElemSendSize = elementList.getMaxElemSize();
  const CFuint wordFormatSize = _intWordFormatSize+1;
  
  // wOffset is initialized with current offset
  MPI_Offset offset = tt.nodesOffset[iType].second;
  vector<MPI_Offset> wOffset(_nbWriters, offset); 
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  // real number of nodes written for this element (can be more than nbNodesInType:
  // for instance pyramids, prisms and bricks will all have 8 nodes)
  const CFuint nbNodesInTypeToWrite = 
    getWriteNbNodesInType(nbNodesInType, geoOrder, dim, isCell);
  cf_assert(nbNodesInTypeToWrite >= nbNodesInType);
    
  const CFuint totalWriteSize = nbElementsInType*nbNodesInTypeToWrite;
  
  // start element list offset (current position)
  tt.elemsOffset[iType].first = offset;
  // end element list offset
  tt.elemsOffset[iType].second = tt.elemsOffset[iType].first + totalWriteSize*wordFormatSize;
  
  CFLog(VERBOSE, "ParWriteSolution::writeElementList() => offsets = [" 
	<<  tt.elemsOffset[iType].first << ", " << tt.elemsOffset[iType].second << "]\n");
  
  SafePtr< vector<CFuint> > globalElementIDs = MeshDataStack::getActive()->getGlobalElementIDs();
  if (isCell == false) {
    const int iTRS = getGlobalTRSID(elements->getName());
    cf_assert(iTRS >= 0);
    globalElementIDs = &(*MeshDataStack::getActive()->getGlobalTRSGeoIDs())[iTRS][iType];
  }
  cf_assert(globalElementIDs->size() >= nbLocalElements);
  
  // insert in the write list the local IDs of the elements
  // the range ID is automatically determined inside the WriteListMap
  CFuint elemID = startElementID;
  const CFuint nbLocalElementsInType = nbLocalElements;
  for (CFuint iElem = 0; iElem < nbLocalElementsInType; ++iElem, ++elemID) {
    elementList.insertElemLocalID(elemID, (*globalElementIDs)[elemID], 0);
  }
  elementList.endElemInsertion(_myRank);
  
  CFLog(VERBOSE, "P[" << _myRank << "] has maxElemSendSize = " << maxElemSendSize << "\n");
  
  // buffer data to send
  vector<CFuint> sendElements(maxElemSendSize, 0);
  vector<CFuint> elementToPrint(maxElemSendSize, 0);
    
  const string writerName = getMethodData().getNamespace() + "_Writers";
  Group& wg = PE::GetPE().getGroup(writerName);
  
  CFint wRank = -1; 
  CFuint rangeID = 0;
  CFuint wSendSize = 0;
  CFuint countElem = 0;
  for (CFuint is = 0; is < nSend; ++is, ++rangeID) {
    bool isRangeFound = false;
    WriteListMap::List elist = elementList.find(rangeID, isRangeFound);
    
    if (isRangeFound) {
      CFuint eSize = 0;
      for (WriteListMap::ListIterator it = elist.first; it != elist.second; ++it, ++eSize) {
	const CFuint localElemID = it->second;
	const CFuint globalElemID = (*globalElementIDs)[localElemID];
	const CFuint sendElemID = globalElemID - countElem;
	const CFuint nbNodes = (isCell) ? elements->getNbNodesInGeo(localElemID) : 
	  (*elements)[iType]->getNbNodesInGeo(localElemID);
	cf_assert(nbNodes <= nbNodesInType);
	
	CFuint isend = sendElemID*nbNodesInType;
	for (CFuint in = 0; in < nbNodesInType; ++in, ++isend) {
	  // fix for degenerated elements (e.g. quads with 2 coincident nodes)
	  const CFuint inID = (in < nbNodes) ? in : in-1;
	  const CFuint localNodeID = (isCell) ? elements->getNodeID(localElemID, inID) : 
	    (*elements)[iType]->getNodeID(localElemID, inID);
	  
	  if (isend >= sendElements.size()) {
	    CFLogInfo(_myRank << " nbNodesInType = " << nbNodesInType
		      << " node isend = " << isend << " , size = " << sendElements.size() << "\n");
	    cf_assert(isend < sendElements.size());
	  }
	  
	  const CFuint globalNodeID = nodes[localNodeID]->getGlobalID();
	  sendElements[isend] = tt.mapNodeID2NodeIDByEType[iType]->find(globalNodeID);
	}
      }
      
      cf_assert(eSize*nbNodesInType <= elementList.getSendDataSize(rangeID));
    }
    
    CFLogDebugMax(_myRank << CFPrintContainer<vector<CFuint> >
		  (" sendElements  = ", &sendElements, nbNodesInType) << "\n");
    
    // size of the data that, in total, those processes will send to the aggregator (writer) with ID=is
    const CFuint sendSize = elementList.getSendDataSize(rangeID);
    cf_assert(sendSize <= sendElements.size());
    cf_assert(sendSize <= elementToPrint.size());
    
    // if the rank corresponds to a writing process, record the size to send for this writer
    if (_isWriterRank && wg.globalRanks[is] == _myRank) {
      wSendSize = sendSize; // this should be the total sendsize in the range
      wRank = is;
    }
    
    // for each send, accumulate data to the corresponding writing process
    // this might allocate too much memory on master processes with OPENMPI
    MPI_Reduce(&sendElements[0], &elementToPrint[0], (int)sendSize,
	       MPIStructDef::getMPIType(&sendElements[0]), MPI_MAX, wg.globalRanks[is], _comm);
    
    // the offsets for all writers with send ID > current must be incremented  
    for (CFuint iw = is+1; iw < wOffset.size(); ++iw) {
      cf_assert(sendSize > 0);
      // offset must be computed taling into account the real number of element nodes to write
      wOffset[iw] += (sendSize/nbNodesInType)*nbNodesInTypeToWrite*wordFormatSize;
      CFLog(DEBUG_MIN, "[" << is << ", " << iw << "] => wOffset = " << wOffset[iw] << ", sendSize = " 
	    << (sendSize/nbNodesInType)*nbNodesInTypeToWrite << "\n");
    }
    
    //reset the all sendElement list to 0
    for (CFuint i = 0; i < maxElemSendSize; ++i) {
      sendElements[i] = 0;
    }
    
    // update the count element for the current element type
    countElem += elementList.getSendDataSize(rangeID)/nbNodesInType;
    
    CFLogDebugMax(_myRank << CFPrintContainer<vector<CFuint> >
		  (" elementToPrint  = ", &elementToPrint, nbNodesInType) << "\n");
    
  } // end sending loop
  
  // need to distinguish between writing on boundary and not
  // need to be able to handle hybrid case ... (with degenerated quads having node[3] = node[2]) 
  
  if (_isWriterRank) { 
    CFLog(VERBOSE, "ParCFmeshBinaryFileWriter::writeElementList() => P[" << _myRank 
	  << "] => offset = " << wOffset[wRank] << "\n");
    
    // each writer can now concurrently write all the collected data (related to one element type)
    cf_assert(wRank >= 0);
    
    // point to the corresponding writing location
    fout->seekp(wOffset[wRank]);
    
    const CFuint nbElementsToWrite = wSendSize/nbNodesInType;
    for (CFuint i = 0; i < nbElementsToWrite; ++i) {
      writeElementConn(*fout, &elementToPrint[i*nbNodesInType],
		       nbNodesInType, geoOrder, dim, isCell);
      *fout << "\n";
    }
    
    MPI_Offset lastpos = fout->tellp();
    MPI_Offset maxpos  = 0;
    MPI_Allreduce(&lastpos, &maxpos, 1, MPIStructDef::getMPIOffsetType(), MPI_MAX, wg.comm);
    
    if (tt.elemsOffset[iType].second != maxpos) {
      CFLog(WARN, "tt.elemsOffset[iType].second != maxpos => " << 
	    tt.elemsOffset[iType].second  <<  " != " << maxpos << "\n");
      cf_assert(tt.elemsOffset[iType].second == maxpos);
    }
    fout->seekp(maxpos);
  }
  
  // reset all the elements to print to 0
  for (CFuint i = 0; i < maxElemSendSize; ++i) {
    elementToPrint[i] = 0;
  }
  
  CFLog(VERBOSE, "ParWriteSolution::writeElementList() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

CFuint ParWriteSolution::getWriteNbNodesInType(const CFuint nbNodes,
					       const CFuint geoOrder,
					       const CFuint dim, 
					       const bool isCell)
{
  /// @todo only 1st order geometry
  cf_assert(geoOrder == CFPolyOrder::ORDER1);
  
  if (isCell) {
    switch(dim) {
      
    case DIM_2D:
      return nbNodes;
      break;
      
    case DIM_3D:
      return (nbNodes > 4) ? 8 : 4;
      break;
      
    default:
      std::string msg = std::string("Wrong dimension. Can only be 2D or 3D: ") +
	Common::StringOps::to_str(dim);
      throw BadValueException(FromHere(),msg);
    }
  }
  else {
    // case for 2D and 3D faces
    cf_assert(!isCell);
    
    switch(dim) {
      
    case DIM_2D:
      cf_assert(nbNodes == 2);
      return 2;
      break;
      
    case DIM_3D:
      switch(nbNodes) {
	
      case 3:  // TRIANGLE
	return 3;
	break;
	
      case 4: // QUADRILATERAL
	return 4;
	break;
	
      default:
	std::string msg = std::string("Wrong number of nodes in 3D for face: ") +
	  Common::StringOps::to_str(nbNodes);
	throw BadValueException(FromHere(),msg);
      }
      break;
      
    default:
      std::string msg = std::string("Wrong dimension. Can only be 2D or 3D: ") +
	Common::StringOps::to_str(dim);
      throw BadValueException(FromHere(),msg);
    }  
  }
  
  // should never get here
  cf_assert(false);
  return nbNodes;
}

//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeElementConn(ofstream& file,
					CFuint* nodeIDs,
					const CFuint nbNodes,
					const CFuint geoOrder,
					const CFuint dim, 
					const bool isCell)
{
  /// @todo only 1st order geometry
  cf_assert(geoOrder == CFPolyOrder::ORDER1);
  
  if (isCell) {
    switch(dim) {
      
    case DIM_2D:
      switch(nbNodes) {
	
      case 3:  // TRIANGLE
	file << std::noshowpos 
	     << setw(_intWordFormatSize) << nodeIDs[0]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[1]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[2]+1;
	break;
	
      case 4: // QUADRILATERAL
	file << std::noshowpos 
	     << setw(_intWordFormatSize) << nodeIDs[0]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[1]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[2]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[3]+1;
	break;
	
      default:
	std::string msg = std::string("Wrong number of nodes in 2D for cell: ") +
	  Common::StringOps::to_str(nbNodes);
	throw BadValueException(FromHere(),msg);
      }
      break;
    
    case DIM_3D:
      switch(nbNodes) {
	
      case 4: // TETRAHEDRON
	file << std::noshowpos 
	     << setw(_intWordFormatSize) << nodeIDs[0]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[1]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[2]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[3]+1;
	break;
	
      case 5: // BRICK with nodes coalesced 5,6,7->4
	file << std::noshowpos 
	     << setw(_intWordFormatSize) << nodeIDs[0]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[1]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[2]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[3]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[4]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[4]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[4]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[4]+1;
	break;
	
      case 6: // BRICK with nodes 2->3 and 6->7 coalesced
	file << std::noshowpos 
	     << setw(_intWordFormatSize) << nodeIDs[0]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[1]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[2]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[2]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[3]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[4]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[5]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[5]+1;
	break;
	
      case 8: // BRICK
	file << std::noshowpos 
	     << setw(_intWordFormatSize) << nodeIDs[0]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[1]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[2]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[3]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[4]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[5]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[6]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[7]+1;
	break;
      default:
	std::string msg = std::string("Wrong number of nodes in 3D for cell: ") +
	  Common::StringOps::to_str(nbNodes);
	throw BadValueException(FromHere(),msg);
      }
      break;
    default:
      std::string msg = std::string("Wrong dimension. Can only be 2D or 3D: ") +
      Common::StringOps::to_str(dim);
      throw BadValueException(FromHere(),msg);
    }
  }
  else {
    // case for 2D and 3D faces
    cf_assert(!isCell);
    
    switch(dim) {
      
    case DIM_2D:
      switch(nbNodes) {
	
      case 2: // LINESEG // boundary faces only
	cf_assert(isCell == false);
	file << std::noshowpos 
	     << setw(_intWordFormatSize) << nodeIDs[0]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[1]+1;
	break;
	
      default:
	std::string msg = std::string("Wrong number of nodes in 2D for face: ") +
	  Common::StringOps::to_str(nbNodes);
	throw BadValueException(FromHere(),msg);
      }
      break;
      
    case DIM_3D:
      switch(nbNodes) {
	
      case 3:  // TRIANGLE
	file << std::noshowpos 
	     << setw(_intWordFormatSize) << nodeIDs[0]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[1]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[2]+1;
	break;
	
      case 4: // QUADRILATERAL
	file << std::noshowpos 
	     << setw(_intWordFormatSize) << nodeIDs[0]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[1]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[2]+1 << " "
	     << setw(_intWordFormatSize) << nodeIDs[3]+1;
	break;
	
      default:
	std::string msg = std::string("Wrong number of nodes in 3D for face: ") +
	  Common::StringOps::to_str(nbNodes);
	throw BadValueException(FromHere(),msg);
      }
      break;
    default:
      std::string msg = std::string("Wrong dimension. Can only be 2D or 3D: ") +
	Common::StringOps::to_str(dim);
      throw BadValueException(FromHere(),msg);
    }  
  }
}

//////////////////////////////////////////////////////////////////////////////
      
bool ParWriteSolution::hasChangedMesh(const TecplotTRSType& tt) const
{
  CFAUTOTRACE;
  
  // use signals monitoring for this in the future
  for (CFuint iType = 0; iType < tt.totalNbNodesInType.size(); ++iType) {
    if (hasChangedMesh(iType, tt)) {
      CFLog(VERBOSE, "ParWriteSolution::hasChangedMesh() => true\n");
      return true;
    }
  }
  
  CFLog(VERBOSE, "ParWriteSolution::hasChangedMesh() => false\n");
  return false;
}
      
//////////////////////////////////////////////////////////////////////////////

bool ParWriteSolution::hasChangedMesh(const CFuint iType, 
				      const TecplotTRSType& tt) const
{
  CFAUTOTRACE;
  
  // use signals monitoring for this in the future
  cf_assert(iType < tt.totalNbNodesInType.size());
  cf_assert(iType < tt.oldNbNodesElemsInType.size());
  if (tt.totalNbNodesInType[iType] != tt.oldNbNodesElemsInType[iType].first) {
    CFLog(VERBOSE, "ParWriteSolution::hasChangedMesh (iType = " << iType << ") => true\n");
    return true;
  }
  
  CFuint elemCount = 0; 
  cf_assert(tt.trs.isNotNull());
  if ((tt.trs->hasTag("inner")) && (tt.trs->hasTag("cell"))) {
    elemCount = (*MeshDataStack::getActive()->getElementTypeData())[iType].getNbTotalElems();
  }
  else {
    const int iTRS = getGlobalTRSID(tt.trs->getName());
    cf_assert(iTRS >= 0);
    elemCount = MeshDataStack::getActive()->getTotalTRSInfo()[iTRS][iType];
  }
  
  if (elemCount != tt.oldNbNodesElemsInType[iType].second) {
    CFLog(VERBOSE, "ParWriteSolution::hasChangedMesh (iType = " << iType << ") => true\n");
    return true;
  }
  
  CFLog(VERBOSE, "ParWriteSolution::hasChangedMesh (iType = " << iType << ") => false\n");
  return false;
}
      
//////////////////////////////////////////////////////////////////////////////

int ParWriteSolution::getGlobalTRSID(const string& name) const
{
  const vector<string>& trsNames = MeshDataStack::getActive()->getTotalTRSNames();
  cf_assert(trsNames.size() > 0);
  for (int i = 0; i < (int)trsNames.size(); ++i) {
    if (trsNames[i] == name) return i;
  }
  CFLog(ERROR, "ParWriteSolution::getGlobalTRSID() => NOT FOUND name = " << name << "\n");
  return -1;
}
      
//////////////////////////////////////////////////////////////////////////////
 
void ParWriteSolution::TecplotTRSType::cleanupMappings() 
{
  for (CFuint i = 0; i < nodesInType.size(); ++i) {
    vector<CFuint>().swap(nodesInType[i]);
  }
  
  for (CFuint i = 0; i < mapNodeID2NodeIDByEType.size(); ++i) {
    delete mapNodeID2NodeIDByEType[i];
  }
}
 
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::storeMappings(SafePtr<TopologicalRegionSet> trs,
				     const CFuint iType,
				     const CFuint nbCellsInType,
				     const CFuint nbNodesInType,
				     const CFuint startIdx, 
				     const CFuint endIdx,
				     vector<bool>& foundGlobalID)
{
  CFLog(VERBOSE, "ParWriteSolution::storeMappings() => start\n");
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  
  TecplotTRSType& tt = *_mapTrsName2TecplotData.find(trs->getName());
  cf_assert(iType < tt.nodesInType.size());
  vector<CFuint>& nodesInType = tt.nodesInType[iType];
  const bool isCell = trs->hasTag("cell");
  SafePtr<TopologicalRegion> tr = (isCell) ? CFNULL : (*trs)[iType];
    
  // find which global nodeIDs are used in the elements of this type
  if (nbCellsInType > 0) {
    cf_assert(nodesInType.size() == 0);
    nodesInType.reserve(nbCellsInType*nbNodesInType); // this array can be oversized
    
    for (CFuint iElem = startIdx; iElem < endIdx; ++iElem) {
      const CFuint nbNodes = (isCell) ? trs->getNbNodesInGeo(iElem) : tr->getNbNodesInGeo(iElem);
      // the following is needed for hybrid meshes for which one TR could contain both 
      // quads and triangles (nbNodesInType=4 but nbNodes can be=3 or =4)
      cf_assert(nbNodes <= nbNodesInType);
      
      for (CFuint in = 0; in < nbNodesInType; ++in) {
	// fix for degenerated elements (e.g. quads with 2 coincident nodes)
	const CFuint inID = (in < nbNodes) ? in : in-1;
	const CFuint localNodeID = (isCell) ? trs->getNodeID(iElem,inID) : tr->getNodeID(iElem,inID);
	nodesInType.push_back(nodes[localNodeID]->getGlobalID());
      }
    }
    
    // sort the vector so we can then remove duplicated nodes
    sort(nodesInType.begin(), nodesInType.end(), std::less<CFuint>());
    
    // remove duplicated nodes
    vector<CFuint>::iterator lastNode = unique(nodesInType.begin(), nodesInType.end());
    nodesInType.erase(lastNode, nodesInType.end());
  }
    
  // now each process has a unique list of global node IDs for the current element type
  // maximum number of unique nodes across all processors for the current element type
  CFuint maxNbNodesInProc = 0;
  CFuint nbLocalNodesInType = nodesInType.size();
  MPIError::getInstance().check
    ("MPI_Allreduce", "ParWriteSolution::setup()",
     MPI_Allreduce(&nbLocalNodesInType, &maxNbNodesInProc, 1, 
		   MPIStructDef::getMPIType(&maxNbNodesInProc), MPI_MAX, _comm));
  
  cf_assert(maxNbNodesInProc >= nbLocalNodesInType);
  
  // temporary buf to store global IDs
  vector<CFuint> buf(maxNbNodesInProc);
  
  // reset the flags for all nodes
  foundGlobalID.assign(foundGlobalID.size(), false);
  
  // use collective broadcast to know all global IDs referencing the current element type
  for (CFuint root = 0; root < _nbProc; ++root) {
    CFuint bufSize = 0;
    if (root == _myRank) {
      for (CFuint i = 0; i < nbLocalNodesInType; ++i) {
	buf[i] = nodesInType[i];
      }
      bufSize = nbLocalNodesInType;
    } 
    
    // communicate the buffer size
    MPIError::getInstance().check
      ("MPI_Bcast", "ParWriteSolution::setup()", 
       MPI_Bcast(&bufSize, 1, MPIStructDef::getMPIType(&bufSize), root, _comm));
    
    // the following check is needed for ensuring consistent collective behaviour
    if (bufSize > 0) {
      cf_assert(bufSize <= foundGlobalID.size());
      
      MPIError::getInstance().check
	("MPI_Bcast", "ParWriteSolution::setup()", 
	 MPI_Bcast(&buf[0], bufSize, MPIStructDef::getMPIType(&buf[0]), root, _comm));
      
      // flag all global IDs that are found
      for (CFuint i = 0; i < bufSize; ++i) {
	foundGlobalID[buf[i]] = true;
      }
    }
  }
  
  cf_assert(iType < tt.mapNodeID2NodeIDByEType.size());
  tt.mapNodeID2NodeIDByEType[iType] = new CFMap<CFuint, CFuint>(nodesInType.size());
  
  // create mapping between global IDs belonging to the current processor
  // and global IDs reordered by element type (starting from 0) 
  long long int etypeNodeID = -1;
  CFuint countMatching = 0;
  for (CFuint globalNodeID = 0; globalNodeID < foundGlobalID.size(); ++globalNodeID) {
    if (foundGlobalID[globalNodeID]) {
      etypeNodeID++; 
      // in TECPLOT IDs start from 1: this is enforced in ParWriteSolution::writeElementConn()
    }
    if (std::binary_search(nodesInType.begin(), nodesInType.end(), globalNodeID)) {
      tt.mapNodeID2NodeIDByEType[iType]->insert(globalNodeID, etypeNodeID);
      countMatching++;
    }
  }
  
  cf_assert(iType < tt.totalNbNodesInType.size());
  tt.totalNbNodesInType[iType] = etypeNodeID+1;
  
  cf_assert(countMatching == nodesInType.size());
  cf_assert(iType < tt.mapNodeID2NodeIDByEType.size());
  tt.mapNodeID2NodeIDByEType[iType]->sortKeys(); 
  
  CFLog(VERBOSE, "ParWriteSolution::storeMappings() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

ParWriteSolution::TecplotTRSType::TecplotTRSType
(const CFuint nbElementTypes, SafePtr<TopologicalRegionSet> in) 
{
  trs = in;
  
  headerOffset.resize(nbElementTypes);
  for (CFuint i = 0; i < nbElementTypes; ++i) {
    headerOffset[i].resize(2, 0);
  }
  
  elemsOffset.resize(nbElementTypes);
  for (CFuint i = 0; i < nbElementTypes; ++i) {
    elemsOffset[i].first = elemsOffset[i].second = 0;
  }
  
  nodesOffset.resize(nbElementTypes);
  for (CFuint i = 0; i < nbElementTypes; ++i) {
    nodesOffset[i].first = nodesOffset[i].second = 0;
  }
  
  oldNbNodesElemsInType.resize(nbElementTypes);
  for (CFuint i = 0; i < nbElementTypes; ++i) {
    oldNbNodesElemsInType[i].first = 0;
    oldNbNodesElemsInType[i].second = 0;
  }
  totalNbNodesInType.resize(nbElementTypes, 0);
  
  nodesInType.resize(nbElementTypes);
  mapNodeID2NodeIDByEType.resize(nbElementTypes);
}
      
//////////////////////////////////////////////////////////////////////////////

ParWriteSolution::TecplotTRSType::~TecplotTRSType() 
{
  cleanupMappings();
}
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeInnerZoneHeader(std::ofstream* fout, 
					    const CFuint iType,
					    ElementTypeData& eType,
					    SafePtr<TopologicalRegionSet> trs) 
{
  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  TecplotTRSType& tt = *_mapTrsName2TecplotData.find(trs->getName());
  
  const string geoType = MapGeoEnt::identifyGeoEntTecplot
    (eType.getNbNodes(), eType.getGeoOrder(), PhysicalModelStack::getActive()->getDim()); 
  
  writeZoneHeader(fout, iType, eType.getShape(), tt.totalNbNodesInType[iType], 
		  eType.getNbTotalElems(), geoType, "\t", false);
  
  if (getMethodData().getAppendAuxData()) {
    *fout << ", AUXDATA TRS=\"" << trs->getName() << "\""
	  << ", AUXDATA Filename=\"" << getMethodData().getFilename().leaf() << "\""
	  << ", AUXDATA ElementType=\"" << eType.getShape() << "\""
	  << ", AUXDATA Iter=\"" << std::noshowpos << setw(_intWordFormatSize) << subSysStatus->getNbIter() << "\""
	  << ", AUXDATA PhysTime=\"";
    fout->precision(14);  
    fout->setf(ios::scientific,ios::floatfield);
    *fout << subSysStatus->getCurrentTimeDim() << "\"";
  }
  *fout << "\n";
}
    
//////////////////////////////////////////////////////////////////////////////
  
void ParWriteSolution::writeZoneHeader(std::ofstream* fout, 
				       const CFuint iType,
				       const string& geoShape,
				       const CFuint nbNodesInType,
				       const CFuint nbElemsInType,
				       const string& geoType,
				       const string& end,
				       const bool isBoundary) 
{
  *fout << "ZONE "
	<< "  T= \"ZONE" << iType << " " << geoShape <<"\""
	<< ", N=" << std::noshowpos << nbNodesInType
	<< ", E=" << std::noshowpos << nbElemsInType
	<< ", F=FEPOINT" << ", ET=" << geoType 
	<< ", SOLUTIONTIME=";
  fout->precision(14);  
  fout->setf(ios::scientific,ios::floatfield);
  const CFreal timeDim = SubSystemStatusStack::getActive()->getCurrentTimeDim();
  const CFreal nbIter  = (CFreal)SubSystemStatusStack::getActive()->getNbIter();
  *fout << ((timeDim > 0.) ? timeDim : nbIter) << end;
}
    
//////////////////////////////////////////////////////////////////////////////
  
    } // namespace TecplotWriter
  
} // namespace COOLFluiD
