// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_ParVector_hh
#define COOLFluiD_Common_ParVector_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Common/MPI/MPICommPattern.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaVector.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace Common {

//////////////////////////////////////////////////////////////////////////////

template < typename T >
class ParVector
{  
public: // functions
  
  /// The communication pattern type
  typedef typename COOLFluiD::Common::ArrayAllocator<T>::IndexType IndexType;

#ifdef CF_HAVE_CUDA
  /// array type
  typedef typename CudaEnv::CudaVector<T> ARRAY;
#else
  /// array type
  typedef typename COOLFluiD::Common::ArrayAllocator<T> ARRAY;
#endif
  
  /// communication pattern
  typedef MPICommPattern<ARRAY> CPATTERN;
  
  /// Constructor.
  /// WARNING: Size parameter is IGNORED!
  ParVector (const std::string& nspaceName, 
	     const T& Init, size_t Size, size_t ESize = 0) : 
    m_data(Init, Size, ESize),
    m_namespace(nspaceName),
    m_init(Init), m_size(Size), m_esize(ESize), m_pattern(CFNULL)
  {     
  }
    
  /// Destructor
  /// Before destructing the class, DoneMPI should be called
  ~ParVector() 
  {
    deletePtr(m_pattern);
  }
  
  /// Return the local size: This is the number of
  /// locally owned points incremented by the number of ghost points
  /// Local operation.
  CFuint size() const {return m_pattern->size();}
  
  /// Returns the total local size (same as size())
  /// Local operation.
  CFuint GetTotalSize () const {return size();}
  
  /// This function returns the global (cross-processes) size of
  /// the underlying parallel array
  /// @return the global size of the parallel array
  CFuint GetGlobalSize() const {return  m_pattern->GetGlobalSize();}
  
  /// This function returns the global (cross-processes) size of
  /// the underlying parallel array
  /// @return the global size of the parallel array
  CFuint GetLocalSize() const {return  m_pattern->GetLocalSize();}
  
  /// Element access (preferred method): 
  /// @param Index   index corresponding to stateID(or nodeID)*stride
  /// Never take the memory location of an element as it
  /// CAN and WILL change on resizing)
  /// [We do it because of the state mess]
  const T & operator () (size_t Index) const {return m_data(Index);}
  
  /// Element access (preferred method)
  /// @param Index   index corresponding to stateID(or nodeID)*stride
  /// Never take the memory location of an element as it
  /// CAN and WILL change on resizing)
  /// [We do it because of the state mess]
  T& operator () (size_t Index) {return m_data(Index);}
  
  /// Make sure we can have up to capacity elements
  /// before needing to allocate
  /// (and possible invalidate pointers & references)
  void reserve (size_t capacity, size_t elementSize, 
		const std::string& nspaceName) 
  {
    if (m_pattern != CFNULL) {deletePtr(m_pattern);}
    m_pattern = new CPATTERN(m_namespace, &m_data, m_init, m_size, m_esize); 
    m_pattern->reserve(capacity, elementSize, nspaceName);
  }
  
  /// begin the synchronization
  void BeginSync() {m_pattern->BeginSync();}
  
  /// end the synchronization
  void EndSync() { m_pattern->EndSync();}
  
  /// Build Sync table
  void BuildGhostMap () { m_pattern->BuildGhostMap();}
  
  /// Returns the list of ghost nodes (by processor rank) to be sent to
  /// another processor
  const std::vector< std::vector< CFuint > >& GetGhostSendList() const
  {return m_pattern->GetGhostSendList();}

  /// Returns the list of ghost nodes (by processor rank) to be received
  /// from another processor
  const std::vector< std::vector< CFuint > >& GetGhostReceiveList() const
  {return m_pattern->GetGhostReceiveList();}
  
  /// Build a continuous global mapping.
  /// Should be called after all points are added,
  /// and after buildMap is called
  /// This function uses indexes if available
  void BuildCGlobal() {m_pattern->BuildCGlobal();}
  
  /// free the continuos global mapping.
  void DestroyIndex() {m_pattern->DestroyIndex();}
  
  /// Lookup the global continuous ID of a local element
  CFuint LocalToCGlobal (CFuint localID) const {return m_pattern->LocalToCGlobal(localID);}
  
  /// Insert a new ghost point
  /// Local operation.
  /// For now, NO Add operations are allowed after
  /// BuildGhostMap is called.
  IndexType AddGhostPoint (IndexType GlobalIndex) {return m_pattern->AddGhostPoint(GlobalIndex);}
  
  /// Insert new local point
  /// Local operation.
  /// For now, NO add operations are allowed after
  /// BuildGhostMap is called.
  IndexType AddLocalPoint (IndexType GlobalIndex) {return m_pattern->AddLocalPoint(GlobalIndex);}
  
  /// Get array 
  Common::SafePtr<ARRAY> getPtr() {return &m_data;}
  
private: // functions
  
  /// This stores the actual element data
  ARRAY m_data;
  
  /// namespace
  std::string m_namespace;
  
  /// initial value
  T m_init; 

  /// size
  size_t m_size;
  
  /// element size
  size_t m_esize;
    
  /// communication pattern
  CPATTERN* m_pattern;
};
      
//////////////////////////////////////////////////////////////////////////////

    } // Common

} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_ParVector_hh
