// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataHandleMPI_hh
#define COOLFluiD_Framework_DataHandleMPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/PE.hh"
#include "Common/ParallelException.hh"
#include "Common/Stopwatch.hh"
#include "Common/SharedPtr.hh"
#include "Common/CFMultiMap.hh"

#include "Framework/GlobalCommTypes.hh"
#include "Framework/GlobalTypeTrait.hh"
#include "Framework/DataHandle.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Implementation for DataHandle for MPI
/// @todo Add special operations (adding point / ghost point / sync)
/// @author Dries Kimpe
/// @author Andrea Lani
/// @author Tiago Quintino
template < class TYPE >
class DataHandle < TYPE , MPI > : public  DataHandleInternal < TYPE , LOCAL >
{

public:

  typedef DataHandleInternal < TYPE , LOCAL > BaseClass;
  typedef typename BaseClass::StorageType   StorageType;
  typedef Common::ParVector<typename Framework::GlobalTypeTrait<TYPE>::GTYPE> GlobalVectorType;
  typedef typename BaseClass::ElemType ElemType;
  
public:

  /// Constructor taking the storage type
  explicit DataHandle (StorageType* const ptr) : BaseClass(ptr), _globalPtr(CFNULL) {}
  
  /// Constructor taking the actual storage types
  DataHandle (StorageType* const ptr, GlobalVectorType *const globalPtr) :
    BaseClass(ptr), _globalPtr(globalPtr) {}
  
  /// Constructor taking void*
  DataHandle (void* ptr, void* globalPtr) :
      BaseClass(static_cast<StorageType*>(ptr)),
      _globalPtr(static_cast<GlobalVectorType*>(globalPtr)) {}

  /// Copy constructor
  DataHandle (const DataHandle<TYPE, MPI> & t) : BaseClass (t), _globalPtr(t._globalPtr) {}

  /// Assignment operator
  const DataHandle <TYPE, MPI>&  operator= (DataHandle <TYPE, MPI> t)
  {
    BaseClass::operator=(t);
    _globalPtr = t._globalPtr;
    return *this;
  }

  /// Assignment operator
  const DataHandle <TYPE, MPI>&  operator= (const TYPE t)
  {
    for (CFuint i = 0; i < _globalPtr->size(); ++i) {
      (*_globalPtr)[i] = t;
    }
    return *this;
  }

  /// Get a pointer to the global data
  typename Framework::GlobalTypeTrait<TYPE>::GTYPE* getGlobalData(CFuint localID) const
  {
    return &(*_globalPtr)(localID);
  }
  
  /// Add a local point
  CFuint addLocalPoint (CFuint GlobalIndex)
  {
    return _globalPtr->AddLocalPoint (GlobalIndex);;
  }

  /// Add a ghost point
  CFuint addGhostPoint (CFuint GlobalIndex)
  {
    if (!Common::PE::GetPE().IsParallel()) {
      throw Common::ParallelException
	(FromHere(),"DataHandle<MPI>:addGhostPoint() => using addGhostPoint in serial run\n");
    }
    CFLogDebugMax("DataHandle<MPI>::addGhostPoint() for global ID = " << GlobalIndex << "\n");
    return _globalPtr->AddGhostPoint (GlobalIndex);
  }
  
  /// Returns the list of ghost nodes (by rank) to be sent to another processor
  const std::vector<std::vector<CFuint> >& getGhostSendList() const
  {
    return _globalPtr->GetGhostSendList();
  }
  
  /// Returns the list of ghost nodes (by rank) to be received from another processor
  const std::vector<std::vector<CFuint > >& getGhostReceiveList() const
  {
    return _globalPtr->GetGhostReceiveList();
  }
  
  /// Set the mapping from global ghost IDs to donor ranks
  void setMapGhost2DonorRanks
    (Common::SharedPtr<Common::CFMultiMap<CFuint, CFuint> >& mapGhost2Donor) 
  {_globalPtr->setMapGhost2DonorRanks(mapGhost2Donor);}
  
  /// Build a continuous global mapping.
  /// Should be called after all points are added,
  /// and after buildMap is called
  /// This function uses indexes if available
  void buildContiguosGlobal()
  {
    _globalPtr->BuildCGlobal();
  }

  /// free the continuos global mapping.
  void freeContiguosGlobal()
  {
    _globalPtr->DestroyIndex();
  }

  /// Lookup the global continuous ID of a local element
  CFuint getContiguosID (CFuint localID) const
  {
    return _globalPtr->LocalToCGlobal(localID);
  }

  /// Build Sync table
  void buildMap(const std::string& algo)
  {
    Common::Stopwatch<Common::WallTime> timer;
    timer.start ();
    _globalPtr->BuildGhostMap(algo);
    timer.stop();
    CFLog(INFO, "DataHandle<MPI>::buildMap(" << algo << ") took " << timer.read() << "s\n");
  }
  
  /// This function returns the global (cross-processes) size of
  /// the underlying parallel array
  /// @return the global size of the parallel array
  CFuint getGlobalSize() const
  {
    return (_globalPtr != CFNULL) ? _globalPtr->GetGlobalSize() : 0;
  }

  /// This function returns the local size of the underlying parallel array
  /// @return the local size of the parallel array
  CFuint getLocalSize() const
  {
    return  (_globalPtr != CFNULL) ? _globalPtr->GetLocalSize() : 0;
  }
  
  /// begin the synchronization
  void beginSync()
  {
    cf_assert(_globalPtr != NULL);
    _globalPtr->BeginSync ();
  }
  
  /// end the synchronization
  void endSync()
  {
    cf_assert(_globalPtr != NULL);
    _globalPtr->EndSync ();
  }
    
  /// execute the synchronization
  void synchronize()
  {
    cf_assert(_globalPtr != NULL);
    _globalPtr->synchronize();
  }

  /// allocate memory dynamically before insertion 
  void reserve (const CFuint Size, 
		const CFuint elementSize, 
		const std::string& nspaceName)
  {
    CFLogDebugMin( "DataHandle<MPI>::reserve() with size=" << Size << "\n");
    _globalPtr->reserve (Size, elementSize, nspaceName);
  }

  /// dump the contents of this array
  void DumpContents ()
  {
    _globalPtr->DumpContents ();
  }
  
  /// dump useful info about the global arrays
  void DumpInfo ()
  {
    CFLog(VERBOSE, "DataHandle<MPI>::DumpInfo() => ptr = " << _globalPtr << "\n");
    CFLog(VERBOSE, "DataHandle<MPI>::DumpInfo() => localSize = " << _globalPtr->GetLocalSize () << "\n");
    CFLog(VERBOSE, "DataHandle<MPI>::DumpInfo() => ghostSize = " << _globalPtr->GetGhostSize () << "\n");
    // _globalPtr->DumpInternalData ();
  }
  
  /// get a pointer to the global array
  Common::SafePtr<typename GlobalVectorType::ARRAY> getGlobalArray() const {return _globalPtr->getPtr();}
  
private:

  /// pointer to array storing data to communicate
  GlobalVectorType* _globalPtr;
};

//////////////////////////////////////////////////////////////////////////////

    }
}

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataHandleMPI_hh
