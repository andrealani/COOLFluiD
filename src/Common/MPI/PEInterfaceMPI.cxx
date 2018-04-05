// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/MPI/PEInterfaceMPI.hh"
#include "Common/MPI/MPIInitObject.hh"
#include "Common/MPI/MPIError.hh"
#include "Common/MPI/MPIStructDef.hh"
#include "Common/CFPrintContainer.hh"

#include <stdlib.h>

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace Common {

//////////////////////////////////////////////////////////////////////////////

static void ThrowMPI (MPI_Comm * Comm, int * Error, ...)
{
  DoMPIError (*Error);
}

//////////////////////////////////////////////////////////////////////////////

PEInterface<PM_MPI>::PEInterface (int * argc, char *** args)
  : InitOK_(false), StopCalled_(false)
{
  CheckMPIStatus(MPI_Init (argc, args));
  
  InitOK_ = true;
  
  CallInitFunctions ();
}
      
//////////////////////////////////////////////////////////////////////////////

PEInterface<PM_MPI>::~PEInterface ()
{
  CallDoneFunctions ();
  
  clearGroups();
  
  CheckMPIStatus(MPI_Finalize());
  
  InitOK_=false;
}

//////////////////////////////////////////////////////////////////////////////

void PEInterface<PM_MPI>::CallInitFunctions ()
{
  cf_assert (!StopCalled_);
  InitContainerType::iterator Cur = InitList_.begin();
  while (Cur != InitList_.end()) {
    CallInitFunction (*Cur);
    ++Cur;
  }
  CFLogDebugMin( "PEInterface<PM_MPI>::CallInitFunctions()\n");
}

//////////////////////////////////////////////////////////////////////////////

void PEInterface<PM_MPI>::CallDoneFunctions ()
{
  InitContainerType::iterator Cur = InitList_.begin();
  while (Cur != InitList_.end()) {
    CallDoneFunction (*Cur);
    ++Cur;
  }
  StopCalled_ = true;
  // Clear the list
  InitList_.clear ();
  CFLogDebugMin( "PEInterface<PM_MPI>::CallDoneFunctions ()\n");
}

//////////////////////////////////////////////////////////////////////////////

void PEInterface<PM_MPI>::CallInitFunction (MPIInitObject * NewObj) const
{
  // AL: check this, I'm not sure if "Default" should be used here 
  NewObj->MPI_Init (GetCommunicator(std::string("Default")));
  CFLogDebugMin( "Calling MPI_Init on " << NewObj << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void PEInterface<PM_MPI>::CallDoneFunction (MPIInitObject * NewObj) const
{
  NewObj->MPI_Done ();
  CFLogDebugMin( "Calling MPI_Done on " << NewObj << "\n"); 
}
//////////////////////////////////////////////////////////////////////////////
      
void PEInterface<PM_MPI>::RegisterInitObject (MPIInitObject * NewObj)
{
  cf_assert (!StopCalled_);
  
  InitList_.push_back (NewObj);
  
  if (InitOK_) {
    /// Initialization is already started, so we can call init immediately
    CallInitFunction (InitList_.back());
  }
}

//////////////////////////////////////////////////////////////////////////////

void PEInterface<PM_MPI>::createGroup(const std::string nsp, 
				      const std::string name, 
				      const std::vector<int>& ranks, 
				      const bool mapRank2Group) 
{
  CFLog(VERBOSE, "PEInterface<PM_MPI>::createGroup [" << name << "] => start\n");
  
  if (m_groups.count(name) == 0) {
    const CFuint nranks = ranks.size();
    cf_assert(nranks > 0);
    
    Group* g = new Group();
    g->globalRanks.resize(nranks);
    
    CFLog(VERBOSE, "PEInterface<PM_MPI>::createGroup() => inserting ranks \t");
    for (CFuint i = 0; i < nranks; ++i) {
      g->globalRanks[i] = ranks[i];
      CFLog(VERBOSE, " " << ranks[i] << " ");
      if (mapRank2Group) {
	m_rank2Group.insert(std::make_pair(ranks[i], name));
      }
    }
    CFLog(VERBOSE, "\n");
    
    g->groupRanks.resize(nranks);
    
    MPI_Group allGroup; 
    MPIError::getInstance().check
      ("MPI_Comm_group", "PEInterface<PM_MPI>::createGroup()", MPI_Comm_group(GetCommunicator(nsp), &allGroup)); 
    
    int rank = 0;
    MPIError::getInstance().check
      ("MPI_Comm_rank", "PEInterface<PM_MPI>::createGroup()", MPI_Comm_rank(GetCommunicator(nsp), &rank)); 
    
    MPIError::getInstance().check
      ("MPI_Group_incl", "PEInterface<PM_MPI>::createGroup()", MPI_Group_incl(allGroup, nranks, &g->globalRanks[0], &g->group));
    
    MPIError::getInstance().check
      ("MPI_Comm_create", "PEInterface<PM_MPI>::createGroup()", MPI_Comm_create(GetCommunicator(nsp), g->group, &g->comm));
    
    // assign the group ranks corresponding to the given global ranks 
    MPIError::getInstance().check
      ("MPI_Group_translate_ranks", "PEInterface<PM_MPI>::createGroup()", 
       MPI_Group_translate_ranks(allGroup, nranks, &g->globalRanks[0], g->group, &g->groupRanks[0])); 
    
    m_groups.insert(std::make_pair(name, g));
    
    // sort the ranks if they are not already sorted (AL: why was this needed???)
    std::sort(g->globalRanks.begin(), g->globalRanks.end()); // AL: recheck this: very tricky
    std::sort(g->groupRanks.begin(), g->groupRanks.end());   // AL: recheck this: very tricky
  }
  else {
    CFLog(WARN, "WARNING: PEInterface<PM_MPI>::createGroup() => group " << name << " already created!\n");
  } 
  
  CFLog(VERBOSE, "PEInterface<PM_MPI>::createGroup() [" << name << "] => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////
      
void PEInterface<PM_MPI>::clearGroups()
{
  /// @pre cannot be called from ~PE() because static PE object is destroyed 
  ///      after MPI_finalize(): it would try to double delete MPI_group's otherwise 
  for (map<string,Group*>::iterator itr = m_groups.begin(); itr != m_groups.end(); ++itr) {
    delete itr->second;
  } 
}
      
//////////////////////////////////////////////////////////////////////////////
  
void PEInterface<PM_MPI>::createNodeToCoreMapping(const std::string nsp, 
						  std::vector<int>& nodeIDs,
						  std::vector<int>& coreIDs, 
						  std::vector<int>& uniqueNodeIDs)
{
  MPI_Comm comm = GetCommunicator(nsp);
  
  int world_size = 0;
  MPI_Comm_size(comm, &world_size);
  
  int world_rank = 0;
  MPI_Comm_rank(comm, &world_rank);
  
  int nodeID = 0;
  int coreID = 0;
  
#ifdef CF_HAVE_IBMSTATIC
  getIDNew(comm, world_size, world_rank, nodeID, coreID);
#else
  getID(comm, world_size, world_rank, nodeID, coreID);
#endif
  
  nodeIDs.resize(world_size, -1);
  coreIDs.resize(world_size, -1);
  vector<int> nodeIDsIn(world_size, -1);
  vector<int> coreIDsIn(world_size, -1);
  
  // set local nodeID and coreID
  nodeIDsIn[world_rank] = nodeID;
  coreIDsIn[world_rank] = coreID;
  
  MPI_Allreduce(&nodeIDsIn[0], &nodeIDs[0], world_size, MPI_INT, MPI_MAX, comm);
  MPI_Allreduce(&coreIDsIn[0], &coreIDs[0], world_size, MPI_INT, MPI_MAX, comm);
  
  set<int> nodeIDsList;
  for (int i = 0; i < world_size; ++i) {
    nodeIDsList.insert(nodeIDs[i]);
  }
  cf_assert(nodeIDsList.size() > 0);
  
  uniqueNodeIDs.reserve(nodeIDsList.size());
  for (set<int>::const_iterator it = nodeIDsList.begin(); it != nodeIDsList.end(); ++it) {
    uniqueNodeIDs.push_back(*it);
  }
  cf_assert(uniqueNodeIDs.size() > 0);
  
  if (world_rank == 0) {
    for (int i = 0; i < world_size; ++i) {
      CFLog(VERBOSE, "rank["<<i<<"] => nodeID["<<nodeIDs[i]<<"], coreID["<<coreIDs[i]<<"]\n");
    } 
  }
  MPI_Barrier(comm);
}
      
//////////////////////////////////////////////////////////////////////////////
      
void PEInterface<PM_MPI>::getID(MPI_Comm comm, const int size, const int rank, 
				int& nodeID, int& coreID)
{
  const int world_size=size;
  const int world_rank=rank;
  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len = 0;
  MPI_Get_processor_name(processor_name, &name_len);
  
  int i = 0;
  int j = 0;
  int N = 0;
  int n = 0;
  int ret = 0;
  const int MSIZE = MPI_MAX_PROCESSOR_NAME+1;
  char str[world_size][MSIZE];
  
  MPI_Gather(&processor_name, MSIZE, MPI_CHAR,str[0], MSIZE, MPI_CHAR, 0, comm);
  
  char DiffNodeName[world_size][MSIZE];
  
  N=1;
  
  if(world_rank==0) {
    strncpy(DiffNodeName[0],str[0],MSIZE);
    for(i=1;i<world_size;i++) {
      n=1;
      for(j=0;j<N;j++) {    
	ret= strncmp(str[i],DiffNodeName[j],MSIZE);     
	n=n*ret;
      }
      if(n!=0) {
	strncpy(DiffNodeName[N],str[i],MSIZE);
	N++;
      }
    }
  }
  
  // MPI_Barrier(comm);  
  MPI_Bcast(&N, 1, MPI_INT, 0, comm );
  // MPI_Barrier(comm);  
  MPI_Bcast(&DiffNodeName, N*MSIZE, MPI_CHAR, 0, comm);
  // MPI_Barrier(comm);  
  
  int NodeNumber = -1;
  for(i=0;i<N;i++) { 
    ret = strncmp(DiffNodeName[i],processor_name,MSIZE); 
    if(ret==0) {
      NodeNumber=i;
    }
  }
  MPI_Barrier(comm);  
  
  int coreIDlocal[world_size]; 
  
  if(world_rank==0) {
    for(j=0;j<N;j++) {
      int n=0;
      for(i=0;i<world_size;i++) {
	ret= strncmp(str[i],DiffNodeName[j],MSIZE);     
	
	if(ret==0) {
	  coreIDlocal[i]=n;
	  n++;
	}
      }
    }
  }
  
  MPI_Bcast(&coreIDlocal[0], world_size, MPI_INT, 0, comm);
  MPI_Barrier(comm);  
  
  nodeID = NodeNumber;
  coreID = coreIDlocal[world_rank];
}
      
//////////////////////////////////////////////////////////////////////////////
  
void PEInterface<PM_MPI>::getIDNew(MPI_Comm comm, const int size, const int rank, 
				   int& nodeID, int& coreID)
{
  const int world_size=size;
  const int world_rank=rank;
  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len = 0;
  MPI_Get_processor_name(processor_name, &name_len);
  
  const int MSIZE = MPI_MAX_PROCESSOR_NAME+1;
  
  char *node_str=(char *)malloc(MSIZE*sizeof(char)); 
  char *last;
  char *token;
  node_str=strndup(processor_name,MSIZE);
  
  while((token=strsep(&node_str," ")))  last=token;
  
  strncpy(processor_name,last,MSIZE);
  
  int i = 0;
  int j = 0;
  int N = 0;
  int n = 0;
  int ret = 0;
  
  char *str=(char *)malloc(world_size*MSIZE*sizeof(char));       
  
  MPI_Gather(&processor_name,MSIZE, MPI_CHAR,str, MSIZE, MPI_CHAR,0,comm);
  char  *DiffNodeName= (char*)malloc(world_size*MSIZE*sizeof(char));
  
  N=1;
  
  if(world_rank==0) {
    strncpy(DiffNodeName,str,MSIZE);
    
    for(i=1;i<world_size;i++) {
      n=1;
      for(j=0;j<N;j++) {    
	ret= strncmp(str+i*MSIZE,DiffNodeName+j*MSIZE,MSIZE);     
	n=n*ret;
      }
      if(n!=0) {
	strncpy(DiffNodeName+N*MSIZE,str+i*MSIZE,MSIZE);
	N++;
      }
    }
  }
  
  MPI_Barrier(comm);  
  MPI_Bcast(&N, 1, MPI_INT, 0, comm );
  MPI_Barrier(comm);  
  
  MPI_Bcast(DiffNodeName, N*MSIZE, MPI_CHAR, 0, comm );
  MPI_Barrier(comm);  
  
  int NodeNumber = 0;
  for(i=0;i<N;i++) { 
    ret= strncmp(DiffNodeName+i*MSIZE,processor_name,MSIZE); 
    if(ret==0) {
      NodeNumber=i;
    }
  }
  MPI_Barrier(comm);  
  
  int *coreID_local=(int *)malloc(world_size*sizeof(int)); 
  
  if(world_rank==0) {
    for(j=0;j<N;j++) {
      int n=0;
      for(i=0;i<world_size;i++) {
	ret= strncmp(str+i*MSIZE,DiffNodeName+j*MSIZE,MSIZE);     
	
	if(ret==0) {
	  coreID_local[i]=n;
	  n++;
	}
      }
    }
  }
  
  MPI_Bcast(coreID_local, world_size, MPI_INT, 0, comm);
  MPI_Barrier(comm);  
  
  nodeID = NodeNumber;
  coreID = coreID_local[world_rank];
  
  free(str);
  free(DiffNodeName);
  free(coreID_local);
  free(node_str);
}
      
//////////////////////////////////////////////////////////////////////////////
  
    }
}
