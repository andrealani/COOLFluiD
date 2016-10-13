// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_LSSIdxMapping_hh
#define COOLFluiD_Framework_LSSIdxMapping_hh

//////////////////////////////////////////////////////////////////////////////

#include <valarray>

#include "Common/COOLFluiD.hh"
#include "Common/NonCopyable.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class provides an index mapping from local IDs to global LSS IDs
/// @author Andrea Lani
class Framework_API LSSIdxMapping : public Common::NonCopyable<LSSIdxMapping> {
public:

  /// Constructor
  LSSIdxMapping() :
    _localToLSSIDs(),
    _isNonLocalRow()
  {
  }

  /// Destructor.
  ~LSSIdxMapping()
  {
  }

  /// Create a mapping from local IDs to global LSS IDS
  /// (sequential case)
  /// @param globalIDs array storing the globall Ids of all the points
  ///                  (all non-ghost) in the only process
  /// @pre globalIDs == localIDs, since this is the sequential case
  void createStdSequential(const std::vector<CFuint>& globalIDs)
  {
    const CFuint nbPoints = globalIDs.size();
    _localToLSSIDs.resize(nbPoints);
    for (CFuint i = 0; i < nbPoints; ++i) {
      _localToLSSIDs[i] = globalIDs[i];
    }

    // in the sequential case all rows are locally owned
    _isNonLocalRow.resize(nbPoints);
    _isNonLocalRow = false;
  }


  void createMapping(const std::valarray<CFuint>& globalIDs,
      const std::valarray<bool>& isNonLocalRow)
  {
    
    const CFuint nbPoints = globalIDs.size();
    _localToLSSIDs.resize(nbPoints);
    _localToLSSIDs = globalIDs;
 //   CFLog(NOTICE,"createMapping size " << nbPoints << "\n");
    // in the sequential case all rows are locally owned
    _isNonLocalRow.resize(nbPoints);
    _isNonLocalRow = isNonLocalRow;
  }

  /// Create a mapping from local IDs to global LSS IDS
  /// @post for the moment a mapping global IDs to global LSS IDs
  ///       is created but could be replaced with a mappping localID to
  ///       global LSS IDs (much more efficient in speed and memory)
  template <class MAPPING_COMPUTER>
  void createMapping(MAPPING_COMPUTER& mappingComputer)
  {
    mappingComputer.computeMapping(_localToLSSIDs,
      _isNonLocalRow);
  }

  /// Get the LSS global ID corresponding to the given local ID
  CFuint getColID(CFuint localID) const
  {
 //   CFLog(NOTICE,"GetCol: localID " << localID << " localToLSSIDs.size() " << _localToLSSIDs.size());
    cf_assert(localID < _localToLSSIDs.size());
 //   CFLog(NOTICE," LSSID " << _localToLSSIDs[localID] << " \n");
    return _localToLSSIDs[localID];
  }

  /// Get the LSS global ID corresponding to the given local ID
  /// and number of equations
  CFuint getColID(CFuint localID, CFuint nbEqs) const
  {
  //  CFLog(NOTICE,"GetCol: localID " << localID << " localToLSSIDs.size() " << _localToLSSIDs.size());
    cf_assert(localID < _localToLSSIDs.size());
  //  CFLog(NOTICE," LSSID " << _localToLSSIDs[localID] << " \n");
    return _localToLSSIDs[localID]*nbEqs;
  }

  /// Get the LSS global ID corresponding to the given local ID
  CFint getRowID(CFuint localID) const
  {
  //  CFLog(NOTICE,"GetRow1: localID " << localID << " localToLSSIDs.size() " << _localToLSSIDs.size());
    cf_assert(localID < _localToLSSIDs.size());
    if (!_isNonLocalRow[localID]) {
  //    CFLog(NOTICE," LSSID " << _localToLSSIDs[localID] << " \n");
      return _localToLSSIDs[localID];
    }
    return -1;
  }

  /// Get the LSS global ID corresponding to the given local ID
  /// and number of equations
  CFint getRowID(CFuint localID, CFuint nbEqs) const
  {
  //  CFLog(NOTICE,"GetRow2: localID " << localID << " localToLSSIDs.size() " << _localToLSSIDs.size());
    cf_assert(localID < _localToLSSIDs.size());
    if (!_isNonLocalRow[localID]) {
      return _localToLSSIDs[localID]*nbEqs;
  //    CFLog(NOTICE," LSSID " << _localToLSSIDs[localID] << " \n");
    }
    return -1;
  }

  CFuint getMappingSize() const
  {
    return _localToLSSIDs.size();
  }

private:

  /// local to global LSS IDs mapping
  std::valarray<CFuint> _localToLSSIDs;

  /// array sstoring flags telling if the current is corresponds to a
  /// non local row entry
  std::valarray<bool> _isNonLocalRow;

}; // end of class LSSIdxMapping

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_LSSIdxMapping_hh
