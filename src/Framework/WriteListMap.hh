// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_WriteListMap_hh
#define COOLFluiD_Framework_WriteListMap_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/CFMultiMap.hh"
#include "Common/CFPrintContainer.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a map for handling a data set to write
/// @author Andrea Lani
class Framework_API WriteListMap {
public:

  /// Constructor.
  WriteListMap();

  /// Destructor.
  ~WriteListMap();

  /// Reserve the map
  void reserve(CFuint nbTypes, CFuint nbSendsPerType, CFuint nbLocalEntries);
  
  /// Fill the map
  void fill(CFuint total, CFuint stride, CFuint& totalToSend);
  
  /// Get the maximum element size to send
  /// @pre this function must be called after fill()
  CFuint getMaxElemSize() const {return _maxElemSize;}
  
  /// Get the minimum element size to send
  /// @pre this function must be called after fill()
  CFuint getMinElemSize() const {return _minElemSize;}
  
  /// Get the number of sends per type
  CFuint getNbSendsPerType() const {return _nbSendsPerType;}
  
  /// Add a range (minumum and maximum global element IDs)
  void addRange(CFuint minID, CFuint maxID) 
  {_minID.push_back(minID); _maxID.push_back(maxID);}
  
  /// Add the number of entries to send
  void addSendDataSize(CFuint dataSize) { _dataSizePerSends.push_back(dataSize); }
  
  /// Get the data to send size
  CFuint getSendDataSize(CFuint rangeID) const
  {
    cf_assert(rangeID < _dataSizePerSends.size());
    return _dataSizePerSends[rangeID];
  }

  /// Insert localID of the element to send in the corresponding range
  void insertElemLocalID(CFuint localElemID,  CFuint globalElemID,  CFuint iType);

  /// End insertion of keys
  void endElemInsertion(CFuint rank);

  /// Get the list of element local IDs to write in the given range
  typedef Common::CFMultiMap<CFuint, CFuint>::MapIterator ListIterator;
  typedef std::pair<ListIterator, ListIterator> List;

  WriteListMap::List find(const CFuint& localElemID, bool& isFound)
  {
    return _rangeToElemLocalID.find(localElemID, isFound);
  }

private: // data

  /// number of types in the storage
  CFuint _nbTypes;

  /// number of sends per type in the storage
  CFuint _nbSendsPerType;

  /// number of local entries in the storage
  CFuint _nbLocalEntries;
  
  /// minimum element size to send
  CFuint _minElemSize;
  
  /// maximum element size to send
  CFuint _maxElemSize;

  /// list of the minimum IDs in send ranges
  std::vector<CFuint> _minID;

  /// list of the maximum IDs in send ranges
  std::vector<CFuint> _maxID;

  /// list of the data sizes of ranges
  std::vector<CFuint> _dataSizePerSends;

  /// mapping between range IDs and local element IDs
  Common::CFMultiMap<CFuint, CFuint> _rangeToElemLocalID;

}; // class WriteListMap

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_WriteListMap_hh
