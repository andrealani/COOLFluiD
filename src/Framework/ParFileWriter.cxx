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

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

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
  _maxBuffSize = 2147479200;
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
  
  CFLog(VERBOSE, "1 ParFileWriter::setWriterGroup() => start\n");

  // the number of writers can be either customized by each derived 
  // algorithm or uniquely defined fixed for all algorithms 
  // by setting CFEnv.NbWriters
  if (_nbWriters == 1) {
    _nbWriters =  CFEnv::getInstance().getVars()->NbWriters;
  }
  
  CFLog(VERBOSE, "2 ParFileWriter::setWriterGroup() => start\n");

  if (_nbWriters > _nbProc) {
    CFLog(WARN, "ParFileWriter::setWriterGroup() => _nbWriters > _nbProc, therefore they are set equal!\n");
    _nbWriters = _nbProc;
  }
  
  CFLog(VERBOSE, "ParFileWriter::setWriterGroup() => _nbWriters = " << _nbWriters << "\n");
  
  // in reality ranks will be user-defined 
  vector<int> ranks;
  vector<int> writerRanks(_nbWriters);
  CFint count = 0;
  bool first = true;
  for (CFuint i = 0; i < _nbWriters; ++i) {    
    const CFuint nbProcPerWriter = _nbProc/_nbWriters;
    const CFuint nranks =  (i < _nbWriters-1) ? nbProcPerWriter : nbProcPerWriter + _nbProc%_nbWriters;
    ranks.resize(nranks);
    for (CFuint r = 0; r < nranks; ++r, ++count) { 
      ranks[r] = count;
      if (ranks[r] == _myRank) {
	// store the ID of the corresponding group
	_myGroupID = i;
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
  
  // create the writers group
  const string writerName = nsp + "_Writers";
  PE::GetPE().createGroup(writerName, writerRanks, false);
  
  CFLog(VERBOSE, "ParFileWriter::setWriterGroup() => " << 
	CFPrintContainer<vector<int> >(" writerRanks  = ",  &writerRanks) << "\n");  
  
  // if (_myRank == 0) {
  //   for (CFuint i = 0; i < _nbProc; ++i) {
  //     cout << "rank [" << i << "] in group " << PE::getGroupName(i) << endl;
  //   }
  // }
  
  // create an info object
  // MPI_Info_create(&_info);
  
  CFLog(VERBOSE, "ParFileWriter::setWriterGroup() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

