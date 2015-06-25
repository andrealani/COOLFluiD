// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <sstream>
#include <fstream>
#include <iostream>
#include "Common/PE.hh"
#include "Common/CFLog.hh"
#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/MPI/MPIStructDef.hh"
#include "Framework/MeshData.hh"
#include "Framework/PartitionerPeriodicTools.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// AL: this will not work with "long long int", needs to be adapted     

void PartitionerPeriodicTools::writePeriodicInfo(const int ndim, 
						 std::string& name0, 
						 std::vector<int> &gidx0, 
						 std::vector<CFreal> &coord0, 
						 std::string& name1, 
						 std::vector<int> &gidx1, 
						 std::vector<CFreal> &coord1)
{
  // get numbers
  
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  const unsigned nproc = Common::PE::GetPE().GetProcessorCount(nsp);
  MPI_Comm comm = Common::PE::GetPE().GetCommunicator(nsp);
  
  std::vector<int> ngidx0(nproc);
  std::vector<int> ngidx1(nproc);
  std::vector<int> displs0(nproc,0);
  std::vector<int> displs1(nproc,0);
  int iprocsend0=gidx0.size();
  int iprocsend1=gidx1.size();
  
  MPI_Gather(&iprocsend0,1,Common::MPIStructDef::getMPIType(&iprocsend0),
	     &ngidx0[0],1,Common::MPIStructDef::getMPIType(&ngidx0[0]),0,comm);
  MPI_Gather(&iprocsend1,1,Common::MPIStructDef::getMPIType(&iprocsend1),
	     &ngidx1[0],1,Common::MPIStructDef::getMPIType(&ngidx1[0]),0,comm);
  for(unsigned int i=1; i<ngidx0.size(); ++i) displs0[i]=displs0[i-1]+ngidx0[i-1];
  for(unsigned int i=1; i<ngidx1.size(); ++i) displs1[i]=displs1[i-1]+ngidx1[i-1];
  const int totgidx0=displs0[displs0.size()-1]+ngidx0[ngidx0.size()-1];
  const int totgidx1=displs1[displs1.size()-1]+ngidx1[ngidx1.size()-1];
  
  // collect gidx to process 0
  std::vector<int> allgidx0;
  std::vector<int> allgidx1;
  if (Common::PE::GetPE().GetRank(nsp)==0) allgidx0.resize(totgidx0);
  if (Common::PE::GetPE().GetRank(nsp)==0) allgidx1.resize(totgidx1);
  MPI_Gatherv(&gidx0[0],gidx0.size(),Common::MPIStructDef::getMPIType(&gidx0[0]),
	      &allgidx0[0],&ngidx0[0],&displs0[0],Common::MPIStructDef::getMPIType(&allgidx0[0]),0,comm);
  MPI_Gatherv(&gidx1[0],gidx1.size(),Common::MPIStructDef::getMPIType(&gidx1[0]),
	      &allgidx1[0],&ngidx1[0],&displs1[0],Common::MPIStructDef::getMPIType(&allgidx1[0]),0,comm);
  
  // collect coord to process 0
  for(unsigned int i=0; i<ngidx0.size(); ++i) { ngidx0[i]*=ndim; displs0[i]*=ndim; }
  for(unsigned int i=0; i<ngidx1.size(); ++i) { ngidx1[i]*=ndim; displs1[i]*=ndim; }
  std::vector<CFreal> allcoord0;
  std::vector<CFreal> allcoord1;
  if (Common::PE::GetPE().GetRank(nsp)==0) allcoord0.resize(totgidx0*ndim);
  if (Common::PE::GetPE().GetRank(nsp)==0) allcoord1.resize(totgidx1*ndim);

  MPI_Datatype MPI_CFREAL = Common::MPIStructDef::getMPIType(&allcoord0[0]);
  MPI_Gatherv(&coord0[0],coord0.size(),MPI_CFREAL,&allcoord0[0],&ngidx0[0],&displs0[0],MPI_CFREAL,0,comm);
  MPI_Gatherv(&coord1[0],coord1.size(),MPI_CFREAL,&allcoord1[0],&ngidx1[0],&displs1[0],MPI_CFREAL,0,comm);
  
  if (Common::PE::GetPE().GetRank(nsp)==0)
  {
    // sort the darn thing and remove duplicates
    std::vector< PeriodicInfoItem > sv0=fillPeriodicInfoItemVector(ndim,allgidx0,allcoord0);
    std::vector< PeriodicInfoItem > sv1=fillPeriodicInfoItemVector(ndim,allgidx1,allcoord1);
    std::sort(sv0.begin(),sv0.end());
    std::sort(sv1.begin(),sv1.end());
    sv0.erase(std::unique(sv0.begin(),sv0.end()),sv0.end());
    sv1.erase(std::unique(sv1.begin(),sv1.end()),sv1.end());
    if (sv0.size()!=sv1.size()) throw Common::BadValueException(FromHere(),"Periodic boundaries: number of states don't match.");

    // open file
    std::ofstream f("periodic.info");
    if (f.fail()) throw Common::BadValueException(FromHere(),"Couldnt open periodic.info for writing.");
    f.precision(15);

    // check constant offset along the coordinates
    std::vector<CFreal> ds(ndim);
    CFout << "Offset of the periodic boundaries:";
    for (int i=0; i<ndim; ++i) { ds[i]=sv1[0].crd[i]-sv0[0].crd[i]; CFout << " " << ds[i]; }
    CFout << "\n";
    for (unsigned int i=0; i<sv0.size(); ++i)
      for (int j=0; j<ndim; ++j)
        if (fabs((sv1[i].crd[j]-sv0[i].crd[j])-ds[j])>1e-10)
          throw Common::BadValueException(FromHere(),"Incompatible periodic TRSs, the distance is not constant.");

    // write
    f << sv0.size() << "\n" << ndim << "\n";
    f << name0 << "\n" << name1 << "\n";
    for (unsigned int i=0; i<sv0.size(); ++i)
    {
      f << sv0[i].idx; for (int j=0; j<ndim; ++j) f << " " << sv0[i].crd[j]; f << "\n";
      f << sv1[i].idx; for (int j=0; j<ndim; ++j) f << " " << sv1[i].crd[j]; f << "\n";
    }
    f.flush();
    f.close();
  }
  
  // others should wait for process 0 to catch up with the writing
  Common::PE::GetPE().setBarrier(nsp);
}

//////////////////////////////////////////////////////////////////////////////

void PartitionerPeriodicTools::meldNodes(const int numtotalnodes,
					 std::vector<int>& which, 
					 std::vector<int>& with, 
					 std::vector<PartitionerData::IndexT>& in, 
					 std::vector<int>& copy_of_original)
{
  if (which.size()!=with.size()) throw Common::BadValueException(FromHere(),"Sizes of 'which' and 'what' do not match.");
  std::vector<int> allidx(numtotalnodes,-1);
  copy_of_original.resize(in.size(),-1);
  copy_of_original.assign(in.begin(),in.end());
  for (int i=0; i<(const int)which.size(); ++i)
    allidx[which[i]]=with[i];
  for (PartitionerData::IndexT i=0; i< (const PartitionerData::IndexT)in.size(); ++i)
    if (allidx[in[i]]!=-1)
      in[i]=allidx[in[i]];
}

//////////////////////////////////////////////////////////////////////////////

// this is a tricky routine: because parmetis's ParMETIS_V3_PartMeshKway partitions elements, not nodes,
// therefore it has to be made sure that the periodic face pairs are connecting elements residing on the same rank
// order in node0 and node1 must match
// this routine changes the output of ParMETIS_V3_PartMeshKway by:
//  - first looking up the involved nodes and stored on which rank would it be
//  - every element is checked and if has a node on the periodic bc then the node is indexed to be on the smallest rank
//  - taking the mpireduce/minimum to globalize
//  - making sure the pairs are on the same nodes
//  - then every element's rank involved with the periodic bcs will be put on the minimum nodes
// this routine is expected to be fragile
void PartitionerPeriodicTools::fixPeriodicEdges(const int numproc, const int numtotalnodes, 
						std::vector<int>& node0, std::vector<int>& node1, 
						std::vector<PartitionerData::IndexT>& part, 
						std::vector<PartitionerData::IndexT>& eptrn, 
						std::vector<PartitionerData::IndexT>& elemNode)
{
  if (node0.size()!=node1.size()) throw Common::BadValueException(FromHere(),"Sizes of node0 and node1 does not match.");
  
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  MPI_Comm comm = Common::PE::GetPE().GetCommunicator(nsp);
  
  // requires TWO! fuckin global vector of ints because allreduce sux reusing the buffer
  std::vector<int> pernodinfo(numtotalnodes,numproc+1);
  std::vector<int> pernodinfo2(numtotalnodes,numproc+1);
  for (int i=0; i<(const int)node0.size(); ++i) pernodinfo[node0[i]]=numproc;
  for (int i=0; i<(const int)node1.size(); ++i) pernodinfo[node1[i]]=numproc;

  // figuring out the ranks of the nodes
  for (int i=0; i<(const int)(eptrn.size()-1); ++i)
    for (int j=eptrn[i]; j<(const int)eptrn[i+1]; ++j)
      if (pernodinfo[elemNode[j]]!=numproc+1)
        pernodinfo[elemNode[j]]=pernodinfo[elemNode[j]]<part[i]?pernodinfo[elemNode[j]]:part[i];

  // making sure the nodes on the perbcs are to be put on the same rank
  pernodinfo2.assign(pernodinfo.begin(),pernodinfo.end());
  MPI_Allreduce( &pernodinfo2[0], &pernodinfo[0], (int)pernodinfo.size(), MPI_INT, MPI_MIN, comm);
  pernodinfo2.clear();
  pernodinfo2.reserve(0);

  for (int i=0; i<(const int)node0.size(); ++i)
  {
    if (pernodinfo[node0[i]]==numproc) throw Common::BadValueException(FromHere(),"Node0 has an invalid entry.");
    if (pernodinfo[node1[i]]==numproc) throw Common::BadValueException(FromHere(),"Node1 has an invalid entry.");
    if (pernodinfo[node0[i]]!=pernodinfo[node1[i]]) {
//std::cout << "***********XXXXXXXXXXXXX******************** FOUND ONE!!!\n" << std::flush;
      int minidx=pernodinfo[node0[i]]<pernodinfo[node1[i]]?pernodinfo[node0[i]]:pernodinfo[node1[i]];
      pernodinfo[node0[i]]=minidx;
      pernodinfo[node1[i]]=minidx;
    }
  }

//std::cout << "004\n" << std::flush;
//sleep(1);

  // check numbers
  int perctr=0;
  for (int i=0; i<(const int)numtotalnodes; ++i) if (pernodinfo[i]!=numproc+1) perctr++;
  if (perctr!=(int)(2u*node0.size())) throw Common::BadValueException(FromHere(),"Nuber of periodic nodes in the global lookup array does not match expected value.");

//std::cout << "005\n" << std::flush;
//sleep(1);

  // and finally change the rank of the required elements
  for (int i=0; i<(const int)(eptrn.size()-1); ++i)
    for (int j=eptrn[i]; j<(const int)eptrn[i+1]; ++j)
      if (pernodinfo[elemNode[j]]!=numproc+1)
        part[i]=part[i]<pernodinfo[elemNode[j]]?part[i]:pernodinfo[elemNode[j]];

//std::cout << "006\n" << std::flush;
//sleep(1);

}

//////////////////////////////////////////////////////////////////////////////

  } // Framework
} // COOLFluiD


