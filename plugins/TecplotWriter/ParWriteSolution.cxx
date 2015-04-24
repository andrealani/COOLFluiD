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

#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Framework/MapGeoEnt.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/DataHandleOutput.hh"
#include "Framework/SubSystemStatus.hh"

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
  options.addConfigOption< std::string>("FileFormat","Format to write Tecplot file."); 
  options.addConfigOption< CFuint >("NbWriters", "Number of writers (and MPI groups)");
  options.addConfigOption< int >("MaxBuffSize", "Maximum buffer size for MPI I/O");
}
      
//////////////////////////////////////////////////////////////////////////////

ParWriteSolution::ParWriteSolution(const std::string& name) : 
  TecWriterCom(name),
  ParFileWriter(),
  socket_nodes("nodes"),
  socket_nstatesProxy("nstatesProxy")
{
  addConfigOptionsTo(this);

  _fileFormatStr = "ASCII";
  setParameter("FileFormat",&_fileFormatStr);
  
  _nbWriters = 1;
  setParameter("NbWriters",&_nbWriters);
  
  _maxBuffSize = 2147479200; // (CFuint) std::numeric_limits<int>::max();
  setParameter("MaxBuffSize",&_maxBuffSize);
}

//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::execute()
{
  CFLog(INFO, "Writing solution to " << getMethodData().getFilename().string() << "\n");
  
  if(_fileFormatStr == "ASCII"){    
    writeToFileStream(getMethodData().getFilename());
  }
  else {
    cf_assert(_fileFormatStr == "BINARY");
    
    ///@todo change this to use the tecplot library
    ///this is slow and NOT portable but at least, it takes less space
    writeToFile("tmp");
    std::string transformFile = "$TECHOME/bin/preplot tmp " + getMethodData().getFilename().string();
    CFout << transformFile << "\n";
    
    OSystem::getInstance().executeCommand(transformFile);
    
    //     writeToBinaryFile();
  }
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
      
void ParWriteSolution::writeToFileStream(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;
  
  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  
  if (!getMethodData().onlySurface()) {
    DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
    DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy = socket_nstatesProxy.getDataHandle();
    
    // this is a handle for the nodal states
    // (which can be stored as arrays of State*, RealVector* or
    // RealVector but they are used as arrays of RealVector*)
    ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];
    
    SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
    datahandle_output->getDataHandles();
    
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
    
    PE::Group& wg = PE::getGroup("Writers");
    char* fileName = const_cast<char*>(filepath.string().c_str()); 
    
    if (_isWriterRank) {
      CFLog(INFO, "Opening file " << filepath.string() << "\n");
      MPI_File_open(wg.comm, fileName, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &_fh); 
    }
    
    writeHeader(&_fh);
    
    if (_isWriterRank) {
      MPI_File_close(&_fh);
    }
    
    exit(1);
    
    // array to store the dimensional states
    RealVector dimState(nbEqs);
    RealVector extraValues; // size will be set in the VarSet
    State tempState;
    
    std::vector<SafePtr<TopologicalRegionSet> > trsList =
      MeshDataStack::getActive()->getTrsList();
    
    
    /*
    for(CFuint iTrs= 0; iTrs < trsList.size(); ++iTrs) {
      SafePtr<TopologicalRegionSet> trs = trsList[iTrs];
      
      if((trs->hasTag("inner")) && (trs->hasTag("cell"))) {
	// we will assume that the number of nodes is the same as
	// the number of states but the connectivity might be different
	//   cf_assert(nodes.size() == nodalStates.getSize());
	
	SafePtr<vector<ElementTypeData> > elementType =
	  MeshDataStack::getActive()->getElementTypeData(trs->getName());
	
	
	// loop over the element types
	// and create a zone in the tecplot file
	// for each element type
	for (CFuint iType = 0; iType < elementType->size(); ++iType) {
	  ElementTypeData& eType = (*elementType)[iType];
	  
	  const CFuint nbCellsInType  = eType.getNbElems();
	  
	  // Tecplot doesn't read zones with 0 elements !!! (this was a problem in parallel)
	  if (nbCellsInType > 0) {
	    const CFuint nbNodesInType  = eType.getNbNodes();
	    std::valarray<CFuint> nodeIDs (nbNodesInType);
	    
	    // find which nodeIDs are used in the elements of this type
	    vector<CFuint> nodesInType;
	    
	    for (CFuint iCell = eType.getStartIdx();
		 iCell < eType.getEndIdx();
		 ++iCell) {
	      cf_assert(nbNodesInType == trs->getNbNodesInGeo(iCell));
	      for (CFuint in = 0; in < nbNodesInType; ++in) {
		nodesInType.push_back(trs->getNodeID(iCell,in));
	      }
	    }
	    
	    // sort the vector so we can then remove duplicated nodes
	    sort(nodesInType.begin(), nodesInType.end(), std::less<CFuint>());
	    // remove duplicated nodes
	    vector<CFuint>::iterator lastNode = unique(nodesInType.begin(),
                                                    nodesInType.end());
	    
	    nodesInType.erase(lastNode,nodesInType.end());
	    
	    // create a map from LocalIDs (in CPU) to IDs per ElementType
	    typedef CFuint LocalID;
	    typedef CFuint IDinType;
	    CFMap<LocalID,IDinType> localToTypeID;
	    localToTypeID.reserve(nodesInType.size());
	    
	    for (CFuint i = 0; i < nodesInType.size(); ++i) {
	      // in the following, + 1 is due Tecplot numbering
	      localToTypeID.insert(nodesInType[i],i + 1);
	    }
	    localToTypeID.sortKeys();
	    
	    // print zone header
	    // one sone per element type
	    fout << "ZONE "
		 << "  T=\"P" << PE::GetPE().GetRank()<< " ZONE" << iType << " " << eType.getShape() <<"\""
		 << ", N=" << nodesInType.size()
		 << ", E=" << nbCellsInType
		 << ", F=FEPOINT"
		 << ", ET=" << MapGeoEnt::identifyGeoEntTecplot
	      (eType.getNbNodes(),
	       eType.getGeoOrder(),
	       PhysicalModelStack::getActive()->getDim()) << flush;
	    fout << ", SOLUTIONTIME=" << subSysStatus->getCurrentTimeDim()  << flush;
	    if (getMethodData().getAppendAuxData())
	      fout << ", AUXDATA CPU=\"" << PE::GetPE().GetRank() << "\""
		   << ", AUXDATA TRS=\"" << trs->getName() << "\""
		   << ", AUXDATA Filename=\"" << getMethodData().getFilename().leaf() << "\""
		   << ", AUXDATA ElementType=\"" << eType.getShape() << "\""
		   << ", AUXDATA Iter=\"" << subSysStatus->getNbIter() << "\""
		   << ", AUXDATA PhysTime=\"" << subSysStatus->getCurrentTimeDim() << "\""
		   << flush;
	    fout  << "\n"  << flush;
	    
	    // print nodal coordinates and stored node variables
	    for (vector<CFuint>::iterator itr = nodesInType.begin();
		 itr != nodesInType.end();
		 ++itr) {
	      
	      // current node
	      const CFuint nodeID = *itr;
	      const Node& currNode = *nodes[nodeID];
	      // node has to be printed with the right length
	      for (CFuint iDim = 0; iDim < dim; ++iDim) {
		fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
		fout << currNode[iDim]*refL << " ";
	      }
	      
	      const RealVector& currState = *nodalStates.getState(nodeID);
	      
	      for (CFuint ieq = 0; ieq < nbEqs; ++ieq) {
		tempState[ieq] = currState[ieq];
	      }
	      
	      const CFuint stateID = nodalStates.getStateLocalID(nodeID);
	      tempState.setLocalID(stateID);
	      // the node is set  in the temporary state
	      tempState.setSpaceCoordinates(nodes[nodeID]);
	      
	      if (getMethodData().shouldPrintExtraValues())
		{
		  // dimensionalize the solution
		  updateVarSet->setDimensionalValuesPlusExtraValues
                (tempState, dimState, extraValues);
		  fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
		  fout << dimState << " " << extraValues << "\n";
		}
	      else
		{
		  // set other useful (dimensional) physical quantities
		  updateVarSet->setDimensionalValues(tempState, dimState);
	      fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
	      fout << dimState << "\n";
		}
	      
	      datahandle_output->printStateData(fout,stateID);
	    }
	    
	    // write  Connectivity
	    for (CFuint iCell = eType.getStartIdx();
		 iCell < eType.getEndIdx();
		 ++iCell) {
	      
	      for(CFuint n = 0; n < nbNodesInType; ++n) {
              nodeIDs[n] = localToTypeID.find(trs->getNodeID(iCell, n));
	      }
	      
	      // write their connectivity
	      MapGeoEnt::writeTecplotGeoEntConn(fout,
						nodeIDs,
						eType.getGeoOrder(),
						PhysicalModelStack::getActive()->getDim());
	      fout << "\n";
	    }
	  }
	}
      } //end if inner cells
    } //end loop over trs
    fout.close();
    */
    
  } // if only surface
  
  // write boundary surface data
  // writeBoundarySurface();
}
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::setup()
{
  CFAUTOTRACE;

  TecWriterCom::setup();
  ParFileWriter::setWriterGroup();
}
      
//////////////////////////////////////////////////////////////////////////////

void ParWriteSolution::writeBoundarySurface()
{
  CFAUTOTRACE;

  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
    socket_nstatesProxy.getDataHandle();

  // AL: this is a handle for the nodal states which can be
  // stored as arrays of State*, RealVector* or RealVector
  // but they are always used as arrays of RealVector*
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];

  // we assume that the element-node connectivity
  // is same as element-state connectivity
  /// @todo tecplot writer should not assume element-to-node
  //connectivity to be the same as element-to-state
  //  cf_assert(nodes.size() == nodalStates.getSize());

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  if (getMethodData().getSurfaceTRSsToWrite().empty())
     return;

  const std::vector<std::string>& surfTRS = getMethodData().getSurfaceTRSsToWrite();
  CFuint countTRToWrite = 0;
  std::vector<std::string>::const_iterator itr = surfTRS.begin();

  for(; itr != surfTRS.end(); ++itr) {
    Common::SafePtr<TopologicalRegionSet> currTrs = MeshDataStack::getActive()->getTrs(*itr);
    const CFuint nbTRs = currTrs->getNbTRs();
    for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
      SafePtr<TopologicalRegion> tr = currTrs->getTopologicalRegion(iTR);
      if (tr->getLocalNbGeoEnts() > 0) {
	countTRToWrite++;
      }
    }
  }

  // AL: the file is written only if there is at least one TR with more than 0 nodes
  // else Tecplot cannot handle it and you have manually to skip the file
  if (countTRToWrite > 0) {
    path cfgpath = getMethodData().getFilename();
    path filepath = cfgpath.branch_path() / ( basename(cfgpath) + "-surf" + extension(cfgpath) );

    Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle =
      Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    ofstream& fout = fhandle->open(filepath);

    //  Tecplot Header
    fout << "TITLE      = Boundary data" << "\n";
    fout << "VARIABLES  = ";
    for (CFuint i = 0; i < dim; ++i) {
      fout << " \"x" << i << '\"';
    }

    if (!getMethodData().onlyCoordinates()) {
      SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
      const vector<std::string>& varNames = updateVarSet->getVarNames();
      cf_assert(varNames.size() == nbEqs);

      for (CFuint i = 0 ;  i < nbEqs; ++i) 
      {
	      std::string n = varNames[i];
	      if ( *n.begin()   != '\"' )  n  = '\"' + n;
	      if ( *n.rbegin()  != '\"' )  n += '\"';
	      fout << " " << n;
      }

      if (getMethodData().shouldPrintExtraValues()) {
	vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
	for (CFuint i = 0 ;  i < extraVarNames.size(); ++i) {
	  fout << " " << extraVarNames[i];
	}
      }
    }

    fout << "\n";

#if 0 // DEBUG
    std::string fileNameTST = getMethodData().getFilename() + ".TEST";
    ofstream tstout(fileNameTST.c_str());
#endif // DEBUG

    // array to store the dimensional states
    RealVector dimState(nbEqs);
    RealVector extraValues; // size will be set in the VarSet
    State tempState;

    std::vector<std::string>::const_iterator itr = surfTRS.begin();
    for(; itr != surfTRS.end(); ++itr) {

      Common::SafePtr<TopologicalRegionSet> currTrs = MeshDataStack::getActive()->getTrs(*itr);

      const CFuint nbTRs = currTrs->getNbTRs();
      for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
	SafePtr<TopologicalRegion> tr = currTrs->getTopologicalRegion(iTR);
	const CFuint nbBFaces = tr->getLocalNbGeoEnts();
	if (nbBFaces > 0) {
	  // get all the unique nodes in the TR
	  vector<CFuint> trNodes;
	  tr->putNodesInTR(trNodes);

	  const CFuint nbAllNodes = trNodes.size();
	  CFMap<CFuint,CFuint> mapNodesID(nbAllNodes);
	  for (CFuint i = 0; i < nbAllNodes; ++i) {
	    mapNodesID.insert(trNodes[i],i+1);
	  }
	  mapNodesID.sortKeys();

	  std::string elemShape;
	  CFuint maxNbNodesInGeo = 2;
	  if (dim == 2) {
	    elemShape = "LINESEG";
	  }
	  else if (dim == 3) {
	    maxNbNodesInGeo = 3;
	    const CFuint nbBFaces = tr->getLocalNbGeoEnts();
	    for (CFuint iGeo = 0; iGeo < nbBFaces; ++iGeo) {
	      maxNbNodesInGeo = max(maxNbNodesInGeo, tr->getNbNodesInGeo(iGeo));
	    }

	    // AL: the maximum number of nodes in face is considered:
	    // an extra virtual node is added if TETRA and QUADRILATERAL
	    // are both present
	    elemShape = (maxNbNodesInGeo == 3) ? "TRIANGLE" : "QUADRILATERAL";
	  }

	  // AL: Tecplot doesn't read zones with 0 elements !!! (this was a problem in parallel)
	  if (tr->getLocalNbGeoEnts() > 0) {
	    // print zone header
	    // one zone per TR
	    fout << "ZONE N=" << nbAllNodes
		 << ", T=\"" << currTrs->getName() << ", TR " << iTR << "\""
		 << ", E=" << tr->getLocalNbGeoEnts()
		 << ", F=FEPOINT"
		 << ", ET=" << elemShape
     << ", SOLUTIONTIME=" << subSysStatus->getCurrentTimeDim()
		 << "\n";

	    vector<CFuint>::const_iterator itr;
	    // print  nodal coordinates and stored nodal variables
	    for (itr = trNodes.begin(); itr != trNodes.end(); ++itr) {
	      // node has to be printed with the right length
	      const CFuint nodeID = *itr;
	      const Node& currNode = *nodes[nodeID];
	      for (CFuint iDim = 0; iDim < dim; ++iDim) {
		fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
		fout << currNode[iDim]*refL << " ";
	      }
	      
	      const RealVector& currState = *nodalStates.getState(nodeID);
	      for (CFuint ieq = 0; ieq < nbEqs; ++ieq) {
		tempState[ieq] = currState[ieq];
	      }

	      if (!getMethodData().onlyCoordinates()) {
		const CFuint stateID = nodalStates.getStateLocalID(nodeID);
		tempState.setLocalID(stateID);
		SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();

		if (getMethodData().shouldPrintExtraValues()) {
		  // dimensionalize the solution
		  updateVarSet->setDimensionalValuesPlusExtraValues
		    (tempState, dimState, extraValues);
		  fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
		  fout << dimState << " " << extraValues << "\n";
		}
		else {
		  // set other useful (dimensional) physical quantities
		  updateVarSet->setDimensionalValues(tempState, dimState);
		  fout.precision(14); fout.setf(ios::scientific,ios::floatfield);
		  fout << dimState << "\n";
		}
	      }
	      else {
		fout << "\n";
	      }
	    }

	    // #if 0 // DEBUG
	    //      tstout << "FACES"  << std::endl;
	    //         // write  Connectivity
	    //         GeomEntList::const_iterator itg;
	    //         for (itg = bfaces->begin(); itg != bfaces->end(); ++itg) {
	    //           const vector<Node*>& nodesInFace = *(*itg)->getNodes();
	    //           vector<Node*>::const_iterator itr;
	    //           for (itr = nodesInFace.begin(); itr != nodesInFace.end();
	    //            ++itr) {
	    //             tstout << (*itr)->getGlobalID() << "-" << (*itr)->getLocalID()
	    //                 << " : ";
	    //           }
	    //           tstout << std::endl;
	    //         }
	    //      tstout << "CONNECTIVITY"  << std::endl;
	    // #endif

	    // write  Connectivity
	    for (CFuint iGeo = 0; iGeo < nbBFaces; ++iGeo) {
	      const CFuint nbNodesInGeo = tr->getNbNodesInGeo(iGeo);
	      for (CFuint in = 0; in < nbNodesInGeo; ++in) {
#if 0 // DEBUG
		tstout << tr->getNodeID(iGeo,in); tstout.flush();
		tstout << " : " << mapNodesID.find(tr->getNodeID(iGeo,in)) << endl;
#endif
		fout << mapNodesID.find(tr->getNodeID(iGeo,in)) << " ";
	      }
	      // if the number of face nodes is less than the
	      // maximum number of nodes in a face of this TR
	      // add an extra dummy node equal to the last "real" one
	      if (nbNodesInGeo < maxNbNodesInGeo) {
		cf_assert(maxNbNodesInGeo == nbNodesInGeo + 1);
		fout << mapNodesID.find(tr->getNodeID(iGeo, (nbNodesInGeo-1))) << " ";
	      }
	      fout << "\n";
	      fout.flush();
	    }
	  }
	}
      }
    }

    fout.close();
  }
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
      
void ParWriteSolution::writeHeader(MPI_File* fh)
{
  CFLog(VERBOSE, "ParWriteSolution::writeHeader() start\n");
  
  if (_myRank == _ioRank) {
    //  Tecplot Header
    // maximum string size to write includes 30 characters 
    MPIIOFunctions::writeKeyValue<char>(fh, "TITLE = ");
    MPIIOFunctions::writeKeyValue<char>(fh, "Unstructured grid data");
    MPIIOFunctions::writeKeyValue<char>(fh, "\nVARIABLES = ");
    
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    for (CFuint i = 0; i < dim; ++i) {
      MPIIOFunctions::writeKeyValue<char>(fh, " \"x" + StringOps::to_str(i) + '\"');
    }
    
    SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
    const vector<std::string>& varNames = updateVarSet->getVarNames();
    cf_assert(varNames.size() == nbEqs);
    
    for (CFuint i = 0 ;  i < nbEqs; ++i)  {
      std::string n = varNames[i];
      if ( *n.begin()  != '\"' )  n  = '\"' + n;
      if ( *n.rbegin() != '\"' )  n += '\"';
      MPIIOFunctions::writeKeyValue<char>(fh, " " + n);
    }
    
    if (getMethodData().shouldPrintExtraValues()) {
      vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
      for (CFuint i = 0 ;  i < extraVarNames.size(); ++i) {
	MPIIOFunctions::writeKeyValue<char>(fh, " " +  extraVarNames[i]);
      }
    }
    
    SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
    vector<string> dhVarNames = datahandle_output->getVarNames();
    for (CFuint i = 0 ;  i < dhVarNames.size(); ++i) {
      MPIIOFunctions::writeKeyValue<char>(fh, " " + dhVarNames[i]);
    }
    MPIIOFunctions::writeKeyValue<char>(fh, "\n");
  }
  
  CFLog(VERBOSE, "ParWriteSolution::writeHeader() end\n");
}
  
//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter
  
} // namespace COOLFluiD
