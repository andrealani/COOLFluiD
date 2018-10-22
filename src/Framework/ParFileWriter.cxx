// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ParFileWriter.hh"
#include "Framework/MeshData.hh"
#include "Common/PE.hh"
#include "Common/CFPrintContainer.hh"
#include "Environment/CFEnvVars.hh"
#include "Environment/CFEnv.hh"
#include "MathTools/CFMat.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

ParFileWriter::ParFileWriter() : 
  FileWriter(),
  _comm(0),
  _fh(),
  _status(),
  _myRank(0),
  _nbProc(0),
  _ioRank(0),
  _myGroupID(0),
  _offset(),
  _isNewFile(false),
  _isWriterRank(false),
  _fileList(),
  _mapFileToStartNodeList()
{
  _nbWriters = 1;
  _nbWritersPerNode = 0;
  _maxBuffSize = 2147479200;
  _firstWithoutSolution = false;
}
    
//////////////////////////////////////////////////////////////////////////////

ParFileWriter::~ParFileWriter()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParFileWriter::setWriterGroup()
{
  CFAUTOTRACE;
 
  CFLog(VERBOSE, "ParFileWriter::setWriterGroup() => start\n");

  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  
  // AL: should I put "Default" here??
  _comm   = PE::GetPE().GetCommunicator(nsp);
  _myRank = PE::GetPE().GetRank(nsp);
  _nbProc = PE::GetPE().GetProcessorCount(nsp);
  cf_assert(_nbProc > 0);

  vector<int> writerRanks;

  (_nbWritersPerNode == 0) ? 
    setDefaultWriters(writerRanks) : setNodeWriters(writerRanks);
  
  // create the writers group
  const string writerName = nsp + "_Writers";
  PE::GetPE().createGroup(nsp, writerName, writerRanks, false);
  
  CFLog(VERBOSE, "ParFileWriter::setWriterGroup() => " << 
	CFPrintContainer<vector<int> >(" writerRanks  = ",  &writerRanks) << "\n");  
  
  CFLog(VERBOSE, "ParFileWriter::setWriterGroup() => end\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void ParFileWriter::setDefaultWriters(vector<int>& writerRanks)
{ 
  CFAUTOTRACE;
  
  cf_assert(_nbWritersPerNode == 0);
  
  // the number of writers can be either customized by each derived 
  // algorithm or uniquely defined fixed for all algorithms 
  // by setting CFEnv.NbWriters
  if (_nbWriters == 1) {
    _nbWriters = CFEnv::getInstance().getVars()->NbWriters;
  }
  
  if (_nbWriters > _nbProc) {
    CFLog(WARN, "ParFileWriter::setDefaultWriters() => _nbWriters > _nbProc, => set equal\n");
    _nbWriters = _nbProc;
  }
  
  CFLog(VERBOSE, "ParFileWriter::setDefaultWriters() => _nbWriters = " << _nbWriters << "\n");
  
  writerRanks.resize(_nbWriters);
  vector<int> ranks;
  CFint count = 0;
  bool first = true;
  for (CFuint i = 0; i < _nbWriters; ++i) {    
    const CFuint nbProcPerWriter = _nbProc/_nbWriters;
    const CFuint nranks =  (i < _nbWriters-1) ? 
      nbProcPerWriter : nbProcPerWriter + _nbProc%_nbWriters;
    ranks.resize(nranks);
    for (CFuint r = 0; r < nranks; ++r, ++count) { 
      ranks[r] = count;
      if (ranks[r] == _myRank) {
	// store the ID of the corresponding group
	_myGroupID = i;   // not used
	cf_assert(first); // sanity check to be sure that each rank is associated to one group
	first = false;
      }
    }
    
    // the first rank in each group will be the master process
    writerRanks[i] = ranks[0];
    if (ranks[0] == _myRank) {
      _isWriterRank = true;
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

void ParFileWriter::setNodeWriters(vector<int>& writerRanks)
{ 
  CFAUTOTRACE;
  
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  
  cf_assert(_nbWritersPerNode > 0);
  
  if (_nbWritersPerNode > _nbProc) {
    CFLog(WARN, "ParFileWriter::setNodeWriters() => _nbWritersPerNode > _nbProc => set equal\n");
    _nbWritersPerNode = 1;
  }
  
  vector<int> nodeIDs;
  vector<int> coreIDs;
  vector<int> uniqueNodeIDs;
  PE::GetPE().createNodeToCoreMapping(nsp, nodeIDs, coreIDs, uniqueNodeIDs);
  
  const CFuint nbNodes= uniqueNodeIDs.size();
  _nbWriters = nbNodes*_nbWritersPerNode;
  
  if (_nbWriters > _nbProc) {
    CFLog(WARN, "ParFileWriter::setNodeWriters() => _nbWriters > _nbProc, => set equal\n");
    _nbWriters = _nbProc;
    setDefaultWriters(writerRanks);
  }
  else {
    CFLog(VERBOSE, "ParFileWriter::setNodeWriters() => _nbWriters = " << _nbWriters << "\n");
    
    // choose _nbWritersPerNode per node
    writerRanks.resize(_nbWriters);
    
    // build reverse global to local nodeID
    CFMap<int, int> global2LocalNodeID(nbNodes);
    for (CFuint i = 0; i < uniqueNodeIDs.size(); ++i) {
      global2LocalNodeID.insert(uniqueNodeIDs[i], i);
    }
    global2LocalNodeID.sortKeys();
    
    // nsp nbNodes can be < MPI_COMM_WORLD number of nodes
    CFMat<int> nodeRanks(nbNodes, _nbWritersPerNode, -1);
    cf_assert(nodeRanks.size() == writerRanks.size());
    vector<CFuint> countNodeRanks(nbNodes, (CFuint)0);
    
    cf_assert(nodeIDs.size() == PE::GetPE().GetProcessorCount(nsp));
    
    CFuint countRanks = 0;
    for (CFuint i = 0; i < nodeIDs.size(); ++i) {
      const CFuint nodeID = nodeIDs[i];
      const CFuint localNodeID = global2LocalNodeID.find(nodeID);
      cf_assert(localNodeID < countNodeRanks.size());
      const CFuint countr = countNodeRanks[localNodeID];
      if (countr < _nbWritersPerNode)  {
	nodeRanks(localNodeID, countr) = i;
	countNodeRanks[localNodeID]++;
	countRanks++;
      }
    }
    
    CFLog(INFO, "ParFileWriter::setNodeWriters() => (node, ranks):\n" << nodeRanks);
    
    cf_assert(countRanks == nodeRanks.size());
    for (CFuint i = 0; i < writerRanks.size(); ++i) {
      writerRanks[i] = nodeRanks[i];
      if (nodeRanks[i] == _myRank) {
	_isWriterRank = true;
      }
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

