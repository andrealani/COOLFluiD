// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"
#include "Framework/WriteListMap.hh"

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

WriteListMap::WriteListMap() :
    _nbTypes(0),
    _nbSendsPerType(0),
    _nbLocalEntries(0),
    _minElemSize(0),
    _maxElemSize(0),
    _minID(),
    _maxID(),
    _dataSizePerSends(),
    _rangeToElemLocalID()
{
}

//////////////////////////////////////////////////////////////////////////////

WriteListMap::~WriteListMap() {}

//////////////////////////////////////////////////////////////////////////////

void WriteListMap::reserve (CFuint nbTypes, CFuint nbSendsPerType, CFuint nbLocalEntries)
{
  _nbTypes = nbTypes;
  _nbSendsPerType = nbSendsPerType;
  _nbLocalEntries = nbLocalEntries;
  _minID.reserve(nbSendsPerType*nbTypes);
  _maxID.reserve(nbSendsPerType*nbTypes);
  _dataSizePerSends.reserve(nbSendsPerType*nbTypes);
  _rangeToElemLocalID.reserve(nbLocalEntries);
}

//////////////////////////////////////////////////////////////////////////////

void WriteListMap::insertElemLocalID(CFuint localElemID,  CFuint globalElemID,  CFuint iType)
{
  const CFuint ns = _minID.size();
  const CFuint startType = iType*_nbSendsPerType;
  bool found = false;
  for (CFuint i = startType; i < ns; ++i) {
    if (globalElemID >= _minID[i] && globalElemID < _maxID[i]) {
      _rangeToElemLocalID.insert(i, localElemID);
      found = true;
      break;
    }
  }
  cf_assert(found);
}

//////////////////////////////////////////////////////////////////////////////

void WriteListMap::endElemInsertion(CFuint rank)
{
  cf_assert(_rangeToElemLocalID.size() == _nbLocalEntries);

  _rangeToElemLocalID.sortKeys();
  //    _rangeToElemLocalID.print();
  
  CFLogDebugMin("in proc " << rank << Common::CFPrintContainer<std::vector<CFuint> >  (" global minID  = ", &_minID) << "\n");
  CFLogDebugMin("in proc " << rank << Common::CFPrintContainer<std::vector<CFuint> >  (" global maxID  = ", &_maxID) << "\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void WriteListMap::fill(CFuint total, CFuint stride, CFuint& totalToSend)
{
  const CFuint nSend = _nbSendsPerType;
  _maxElemSize = (total/nSend + total%nSend)*stride;
  _minElemSize = (total/nSend)*stride;
  
  CFuint minElemID = 0;
  for (CFuint is = 0; is < nSend; ++is) {
    // first send per type has maximum size
    // all the other sends (is != 0) per iType have minimum size
    const CFuint currSendSize = (is == 0) ? _maxElemSize : _minElemSize;
    
    // add the global size of the element data (each one of size=nodesPlusStates)
    // to communicate during this send
    addSendDataSize(currSendSize);
    
    const CFuint maxElemID = minElemID + currSendSize/stride;
    
    // insert the range
    addRange(minElemID, maxElemID);
    
    // set the new minimum as the old maximum
    minElemID = maxElemID;
    totalToSend += currSendSize;
  }
}
    
//////////////////////////////////////////////////////////////////////////////
    
  } // namespace Framework
} // namespace COOLFluiD
