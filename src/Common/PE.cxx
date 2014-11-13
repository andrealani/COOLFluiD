// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"

#ifdef CF_HAVE_MPI
#  include "Common/MPI/PEInterfaceMPI.hh"
#else
#  include "Common/SERIAL/PEInterfaceSERIAL.hh"
#endif // CF_HAVE_MPI

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

// initialize the static data
PEInterface<>      * PE::m_curr_PE = CFNULL;
bool                 PE::m_is_init = false;
char               * PE::m_command_workers = CFNULL;
WorkerStatus::Type   PE::m_current_status = WorkerStatus::NOT_RUNNING;

#ifdef CF_HAVE_MPI
std::map<CFint,CFint>    PE::m_rank2Group;
std::vector<PE::Group*>  PE::m_groups; 
#endif

//////////////////////////////////////////////////////////////////////////////

bool PE::IsInitialised ()
{
  return m_is_init;
}

//////////////////////////////////////////////////////////////////////////////

void PE::InitPE (int * argc, char *** args)
{
    cf_assert (m_curr_PE == CFNULL);

    m_command_workers = new char[strlen(*args[0]) + 1];
    strcpy(m_command_workers, *args[0]);

    m_curr_PE = new PEInterface<> (argc, args);
    cf_assert (m_curr_PE != CFNULL);
    m_is_init = true;
}

//////////////////////////////////////////////////////////////////////////////

PEInterface<> & PE::GetPE ()
{
  cf_assert(m_is_init);
  cf_assert (m_curr_PE != CFNULL);
  return *m_curr_PE;
}

//////////////////////////////////////////////////////////////////////////////

void PE::DonePE ()
{
  cf_assert(m_curr_PE != CFNULL);

  // must be first to make sure all destructors dependent of PE acknowledge MPI is down
  m_is_init = false;
  deletePtr(m_curr_PE);
  deletePtrArray(m_command_workers);
}

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
MPI_Comm PE::spawn(unsigned int count, const char * hosts)
{
 MPI_Info info;
 MPI_Comm comm;
 char command[] = "/nobackup/st/gasper/coolfluid/x86_64/optim/plugins/ClientServer/server/app_server";
 int myRank = m_curr_PE->GetRank();

 MPI_Info_create(&info);

 if(hosts != CFNULL && strlen(hosts) != 0)
 {
  char hostKey[] = "host";
  char * hostsNonConst;

  hostsNonConst = new char[strlen(hosts) + 1];
  strcpy(hostsNonConst, hosts);

  CFLogNotice("Spawning " << count << " worker(s) on the following host(s): "
    << hosts << ".\n");
  MPI_Info_set(info, hostKey, hostsNonConst);

  deletePtrArray(hostsNonConst);
 }
 else
 {
  // not giving the "host" value should be sufficient (localhost is the default
  // when no value is provided)
  // this does not work here since we pass a hostfile : MPI uses this file 
  // to determine which hosts to use if no "host" value is provided...
  // setting "host" explicitly to "localhost" prevents this unwanted behaviour
   char hostKey[] = "host";
   char localhostKey[] = "localhost";
   MPI_Info_set(info, hostKey, localhostKey);
  CFLogNotice("Spawning " << count << " worker(s) on local host.\n");
 }

 MPI_Comm_spawn(command,        // command to run
                MPI_ARGV_NULL,  // arguments to the command
                count,          // number of processes
                info,           // infos
                myRank, // manager (root) rank
                m_curr_PE->GetCommunicator(),
                &comm,
                MPI_ERRCODES_IGNORE);

 return comm;
}
#endif

//////////////////////////////////////////////////////////////////////////////

void PE::setCurrentStatus(WorkerStatus::Type status)
{
 cf_assert ( WorkerStatus::Convert::is_valid(status) );
 m_current_status = status;
}

//////////////////////////////////////////////////////////////////////////////

WorkerStatus::Type PE::getCurrentStatus()
{
 return m_current_status;
}

//////////////////////////////////////////////////////////////////////////////

void PE::createGroup(const std::vector<int>& ranks, const bool mapRank2Group) 
{
  const CFuint nranks = ranks.size();
  cf_assert(nranks > 0);
  
  Group* g = new Group();
  g->globalRanks.resize(nranks);
  
  const int groupID = m_groups.size();
  // std::cout << "inserting ranks ";
  for (CFuint i = 0; i < nranks; ++i) {
    g->globalRanks[i] = ranks[i];
    // std::cout << ranks[i] << " ";
    if (mapRank2Group) {
      m_rank2Group.insert(std::make_pair(ranks[i], groupID));
    }
  }
  // std::cout << "\n";
  
  g->groupRanks.resize(nranks);
  
  MPI_Group allGroup; 
  MPI_Comm_group(MPI_COMM_WORLD, &allGroup); 
  
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  MPI_Group_incl(allGroup, nranks, &g->globalRanks[0], &g->group); 
  MPI_Comm_create(MPI_COMM_WORLD, g->group, &g->comm);
  // assign the group ranks corresponding to the given global ranks 
  MPI_Group_translate_ranks(allGroup, nranks, &g->globalRanks[0],
			    g->group, &g->groupRanks[0]); 
  
  m_groups.push_back(g);
}
    
//////////////////////////////////////////////////////////////////////////////

  } // Common

} // COOLFluiD

