// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/MPI/PEInterfaceMPI.hh"
#include "Common/MPI/MPIInitObject.hh"
#include "Common/MPI/MPIError.hh"

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
  
    }
}
