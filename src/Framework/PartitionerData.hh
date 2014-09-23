// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_PartitionerData_hh
#define COOLFluiD_Framework_PartitionerData_hh

//////////////////////////////////////////////////////////////////////////////

#include <parmetis.h>

#include "Common/NonCopyable.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class groups some data to be exchanged between the @see MeshPartitioner
/// and a mesh reader
/// @author Andrea Lani
struct Framework_API PartitionerData : public Common::NonCopyable<PartitionerData> {
  
  // AL: the following is needed for ensuring backward compatibility with ParMETIS 3.1
#if PARMETIS_MAJOR_VERSION==4
  typedef idx_t IndexT;
  typedef real_t RealT;
#else
  typedef CFint IndexT;
  typedef float RealT;
#endif
  
  /// constructor
  PartitionerData() {}

  /// destructor
  ~PartitionerData() {}
  
  /// dimensionality of the mesh
  CFint ndim;

  /// array to store some element distribution info
  std::vector<IndexT> elmdist;

  /// array to store the size of element-node connectivity
  std::vector<IndexT> sizeElemNodeVec;
  
  /// array to store the size of element-state connectivity
  std::vector<IndexT> sizeElemStateVec;
  
  /// array to store the element-node connectivity
  std::vector<IndexT> elemNode;
  
  /// array to store the element-state connectivity
  std::vector<IndexT> elemState;
  
  /// array to store the element node pointers
  std::vector<IndexT> eptrn;
  
  /// array to store the element state pointers
  std::vector<IndexT> eptrs;
  
  /// array to store the processor IDs of the locally stored
  /// nodes after the call to the MeshPartitioner
  std::vector<IndexT>* part;
};
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_PartitionerData_hh
