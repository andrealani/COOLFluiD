// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <sstream>

#include "Common/Stopwatch.hh"
#include "Common/SwapEmpty.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"
#include "Framework/ParMetis.hh"
#include "Framework/MeshData.hh"
#include "Framework/PartitionerPeriodicTools.hh"

/////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Framework {

/////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ParMetis, MeshPartitioner, FrameworkLib, 1>
parMetisProvider("ParMetis");

/////////////////////////////////////////////////////////////////////////////

void ParMetis::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< int >("NCommonNodes", "Parmetis parameter for mesh to graph conversion");
  options.addConfigOption< int >("RND","Random seed to use");
  options.addConfigOption< int >("Options","Parmetis options parameter");
}

/////////////////////////////////////////////////////////////////////////////

ParMetis::ParMetis (const std::string & S)
  : MeshPartitioner(S),
    eptr(),
    eidx(),
    elmdist(),
    part()
{
  addConfigOptionsTo(this);

  IN_NCommonNodes_ = 2;
  setParameter("NCommonNodes",&IN_NCommonNodes_);

  IN_Options_ = 0;
  setParameter("Options",&IN_Options_);

  IN_RND_ = 15;
  setParameter("RND",&IN_RND_);
}

//////////////////////////////////////////////////////////////////////////////

ParMetis::~ParMetis ()
{
}

//////////////////////////////////////////////////////////////////////

void ParMetis::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // the default number of common nodes should be 2 in 2D and 3 in 3D
  IN_NCommonNodes_ = 2;

  MeshPartitioner::configure(args);
}

/////////////////////////////////////////////////////////////////////////////

void ParMetis::doPartition(PartitionerData& pData)
{
  CFAUTOTRACE;

  int CommRank;
  int CommSize;

  MPI_Comm_size (Communicator_, &CommSize);
  MPI_Comm_rank (Communicator_, &CommRank);

  // elmdist, part, eidx, eptr was filled by FillInput()
  PartitionerData::IndexT ncon = 1;
  PartitionerData::IndexT options[] = {1,0,15};
  
  options[0] = 0;
  options[1] = IN_Options_;
  options[2] = IN_RND_;
  
  std::vector<PartitionerData::RealT> ubvec(ncon, 1.05);
  std::vector<PartitionerData::RealT> tpwgts (ncon*CommSize, 1.0/(PartitionerData::RealT)(CommSize));

  PartitionerData::IndexT weightflag=0;
  PartitionerData::IndexT numflag = 0;
  PartitionerData::IndexT ncommonnodes = IN_NCommonNodes_;
  PartitionerData::IndexT edgecut = 0;
  PartitionerData::IndexT* idxdummy = NULL;
  
  CFLogDebugMin( "Calling ParMetis::doPartition()\n");
  Common::Stopwatch<Common::WallTime> MetisTimer;

  // resize output array (are we sure that the local node size will not exceed this ?)
  pData.part->resize(pData.elmdist[CommRank+1] - pData.elmdist[CommRank]);

  // melding periodic nodes
  std::string name0("FILE_NOT_EXISTS"),name1("FILE_NOT_EXISTS");
  std::vector<int> idx0(0),idx1(0);
  std::vector<CFreal> coord0(0),coord1(0);
  std::vector<int> copy_orig(0);
  int locmaxidx=0, totmaxidx;
  PartitionerPeriodicTools::readPeriodicInfo(pData.ndim,name0,idx0,coord0,name1,idx1,coord1);
  if (name0!="FILE_NOT_EXISTS") {
    for (PartitionerData::IndexT i=0; i< (const PartitionerData::IndexT)pData.elemNode.size(); ++i) {
      locmaxidx = (pData.elemNode[i] > locmaxidx) ? pData.elemNode[i] : locmaxidx;
    }
    
    const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
    
    MPI_Allreduce(&locmaxidx,&totmaxidx,1,MPI_INT,MPI_SUM,PE::GetPE().GetCommunicator(nsp));
    PartitionerPeriodicTools::meldNodes(totmaxidx+1,
					idx1,
					idx0,
					pData.elemNode,
					copy_orig);
    
    //    pData.part->assign(pData.part->size(),-1);
    //    std::cout << "PARTSIZE: ";
    //    for (int i=0; i<pData.part->size(); ++i) std::cout << " " << (*pData.part)[i];
    //    std::cout << "\n" << std::flush; sleep(1);
  }
    
  CFLogNotice("ParMetis: ncommonnodes = " << ncommonnodes << "\n");
  MetisTimer.start ();
  PartitionerData::IndexT nbPartitions = (PartitionerData::IndexT)CommSize;
  ParMETIS_V3_PartMeshKway (&pData.elmdist[0], // distribution of the elements (= for every cpu)
			    &pData.eptrn[0],  // contains for each element index of the element nodes
			    &pData.elemNode[0],    // element nodes
			    idxdummy,     // weight of the elements // note here a big difference with ParMETIS 3.1
			    &weightflag,  // 0 -> no weights
			    &numflag,     // numbering starts at index 0
			    &ncon,       // number of weights on each vertex
			    &ncommonnodes,// connectivity degree
			    &nbPartitions,  // Number of partitions
			    &tpwgts[0],      // Vertex weight distribution
			    &ubvec[0],      // Imbalance tolerance
			    &options[0],    // Options
			    &edgecut,       // *output* Partition quality
			    &(*pData.part)[0],  // *output* element ranks, parmetis manual is ambivalent
			    &Communicator_);
  
  // putting back original element connectivity
  /// @todo elements around the two sides of the periodic bc should be on the same rank
  // all the elements around
  if (name0!="FILE_NOT_EXISTS")
 {
    pData.elemNode.assign(copy_orig.begin(),copy_orig.end());
    //std::cout << "XXXXXXXXXXTESTNUM: " << pData.elmdist.size() << " x-x-x-x " << pData.eptrn.size() << "\n" << std::flush;
    PartitionerPeriodicTools::fixPeriodicEdges(CommSize,totmaxidx+1,idx0,idx1,*pData.part,pData.eptrn,pData.elemNode);
    //    std::cout << "PARTSIZE p" << CommRank << " n" << pData.part->size() << ": ";
    //    for (int i=0; i<pData.part->size(); ++i) std::cout << " " << (*pData.part)[i];
    //    std::cout << "\n" << std::flush; sleep(1);
    //    std::cout << "PARTSIZE: ";
    //    for (int i=0; i<pData.elemNode.size(); ++i) std::cout << " " << (pData.elemNode)[i];
    //    std::cout << "\n" << std::flush; sleep(1);
    //    std::cout << "PARTSIZE: " << pData.part->size() << " " << pData.elemNode.size() << "\n" << std::flush; sleep(1);
  }

  MetisTimer.stop ();
  CFLog(NOTICE, "ParMetis::doPartition() took " << MetisTimer << "\n");
}

/////////////////////////////////////////////////////////////////////////////

    }
}
