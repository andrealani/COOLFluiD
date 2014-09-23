// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef PARTITIONPERIDOCITOOLS_HH
#define PARTITIONPERIDOCITOOLS_HH

#include <sstream>
#include <fstream>

#include "Common/PE.hh"
#include "Common/CFLog.hh"
#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Framework/TopologicalRegionSet.hh"
#include "Framework/State.hh"
#include "Framework/DataHandle.hh"
#include "Framework/GlobalCommTypes.hh"
#include "Framework/PartitionerData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Class providing access to some functionalities
class PartitionerPeriodicTools {
public:
  
  /// @todo: Multiple periodicities could be set as one big pair of TRSs and if the user could supply the number of constant-offsets
  /// than a fairly simple algorithm could handle it two lookup node-paris safely.
  /// Then multiple node collapse would show up in the sorted list and the first one would always be the same if coordinate ordered.

  // helper class for the periodic processing
  struct PeriodicInfoItem{
    PeriodicInfoItem(): ndim(-1), idx(-1) { crd[0]=0.; crd[1]=0.; crd[2]=0.; }
    PeriodicInfoItem(int _ndim, int _idx, CFreal *_coord): ndim(_ndim), idx(_idx) { for(int i=0; i<(const int)_ndim; ++i) crd[i]=_coord[i]; }
    PeriodicInfoItem(const PeriodicInfoItem& pii): ndim(pii.ndim), idx(pii.idx) { crd[0]=pii.crd[0]; crd[1]=pii.crd[1]; crd[2]=pii.crd[2]; }
    int ndim;
    int idx;
    CFreal crd[3];
    bool operator < (const PeriodicInfoItem& a) const { /*return idx  < a.idx;*/ for (int i=0; i<ndim; i++) if (fabs(crd[i]-a.crd[i])>1e-10) return crd[i] < a.crd[i]; return true; }
    bool operator ==(const PeriodicInfoItem& a) const { /*return idx == a.idx;*/ for (int i=0; i<ndim; i++) if (fabs(crd[i]-a.crd[i])>1e-10) return false; return true; }
  };
  
  static void writePeriodicInfo(const int ndim, std::string& name0, std::vector<int> &gidx0,
			 std::vector<CFreal> &coord0, std::string& name1, std::vector<int> &gidx1, std::vector<CFreal> &coord1);
  
  static void readPeriodicInfo(const int ndim, std::string& name0, std::vector<int> &gidx0, std::vector<CFreal> &coord0,
			       std::string& name1, std::vector<int> &gidx1, std::vector<CFreal> &coord1)
  {
    // open file
  std::ifstream f("periodic.info");
  if (f.fail())
  {
    name0="FILE_NOT_EXISTS";
    name1="FILE_NOT_EXISTS";
    gidx0.resize(0);
    gidx1.resize(0);
    coord0.resize(0);
    coord1.resize(0);
    return;
  }

  // read num and dimension
  int dim,nstate;
  f >> nstate;
  f >> dim;
  f >> name0;
  f >> name1;

  // set arrays
  gidx0.resize(nstate);       gidx0.reserve(nstate);
  gidx1.resize(nstate);       gidx1.reserve(nstate);
  coord0.resize(nstate*ndim); coord0.reserve(nstate*ndim);
  coord1.resize(nstate*ndim); coord1.reserve(nstate*ndim);

  // read finally
  f.precision(15);
  for (int i=0; i<nstate; ++i)
  {
    f >> gidx0[i]; for (int j=0; j<(const int)ndim; ++j) f >> coord0[i*ndim+j];
    f >> gidx1[i]; for (int j=0; j<(const int)ndim; ++j) f >> coord1[i*ndim+j];
  }
  f.close();
  if (ndim!=dim) throw Common::BadValueException(FromHere(),"Periodic BC: dimension mismatch between computation and periodic.info file.");
  }
  
  static void meldNodes(const int numtotalnodes, std::vector<int>& which, std::vector<int>& with, std::vector<PartitionerData::IndexT>& in, std::vector<int>& copy_of_original);
  
  static std::vector<PeriodicInfoItem> fillPeriodicInfoItemVector(const int ndim, std::vector<int>& idx, std::vector<CFreal>& coord)
  {
    std::vector<PartitionerPeriodicTools::PeriodicInfoItem> r(0);
    r.reserve(idx.size());
    for (int i=0; i<(const int)idx.size(); ++i) r.push_back(PeriodicInfoItem(ndim,idx[i],&coord[ndim*i]));
    return r;
  }
  
  // meaning of pullglobalidx: true==getLocalID or false==getGlobalID
  // it only includes updatable nodes
  static std::vector<PeriodicInfoItem> fillPeriodicInfoItemVector(Common::SafePtr<TopologicalRegionSet> trs, Framework::State** states, bool pullglobalidx=true)
  {
    std::vector<PeriodicInfoItem> r(0);
    Common::SafePtr< std::vector<CFuint> > const stateid = trs->getStatesInTrs();
    const int nstate=stateid->size();
    r.reserve(nstate);
    for (int is=0; is<nstate; ++is) {
      const CFuint sid = (*stateid)[is];
      State* s = states[sid];
      Node& scoord = s->getCoordinates();
      
      if (s->isParUpdatable()) {
      if (pullglobalidx) { r.push_back(PeriodicInfoItem(scoord.size(),s->getGlobalID(),scoord.ptr())); }
      else               { r.push_back(PeriodicInfoItem(scoord.size(),s->getLocalID(), scoord.ptr())); }
      }
    }
    return r;
  }
  
  static void fixPeriodicEdges(const int numproc, const int numtotalnodes, std::vector<int>& node0, std::vector<int>& node1, 
			std::vector<PartitionerData::IndexT>& part, std::vector<PartitionerData::IndexT>& eptrn, std::vector<PartitionerData::IndexT>& elemNode);
  
  // pushes items from 'a' to 'r'!!!
  template <typename T>
  static std::vector< T > findCommonNodes(std::vector< T >& a, std::vector< T >& b)
  {
    std::vector< T > r(0);
    if (a.size()==0) return r;
    if (b.size()==0) return r;
    
    // sort for intersection
    std::sort(a.begin(),a.end());
    std::sort(b.begin(),b.end());
    
    // darn set_intersection does not work correctly, so bruteforcing
    // could be done more smartly with linear tracking algorithm
    for (int i=0; i<(const int)a.size(); ++i)
      {
	for (int j=0; j<(const int)b.size(); ++j)
	  if (a[i]==b[j])
	    {
	      r.push_back(a[i]);
	      break;
	    }
      }
    if (r.size()!=0) std::sort(r.begin(),r.end());
    return r;
  }
  
};

//////////////////////////////////////////////////////////////////////////////

  } // Framework
} // COOLFluiD

#endif // PARTITIONPERIDOCITOOLS_HH
