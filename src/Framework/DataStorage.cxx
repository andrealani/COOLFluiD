// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/DataStorage.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

//////////////////////////////////////////////////////////////////////////////

DataStorage::DataStorage()
{
}

//////////////////////////////////////////////////////////////////////////////

DataStorage::~DataStorage()
{
  if(!m_dataStorage.empty()) {
    vector<std::string> eraseList;

    MapType::iterator itr = m_dataStorage.begin();
    for(; itr != m_dataStorage.end(); ++itr) {
      eraseList.push_back(itr->first);
    }

    vector<std::string>::iterator jtr = eraseList.begin();
    for(; jtr != eraseList.end(); ++jtr) {
      removeDataPtr(*jtr);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::string DataStorage::dump() const
{
CF_DEBUG_POINT;

  ostringstream res;
  res << "DataStorage dump ...\n";
  for( MapType::const_iterator itr = m_dataStorage.begin(); itr != m_dataStorage.end(); ++itr)
  {
     res << "\t[" << itr->first << "] [" << itr->second << "]\n";
  }
  return res.str();
}

//////////////////////////////////////////////////////////////////////////////

}  //  namespace Framework
}  //  namespace COOLFluiD
