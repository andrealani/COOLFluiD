// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iterator>

#include "Framework/InterpolatorRegister.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

InterpolatorProperties InterpolatorRegister::getInterpolatorProperties(const InterpolatorID& id) const
  throw (Common::NoSuchValueException)
{
  if(id >= _database.size()) {
    throw Common::NoSuchValueException (FromHere(),"No such InterpolatorID present.");
  }
  return InterpolatorProperties(_database[id].first,
                                _database[id].second.first,
                                _database[id].second.second,
                                _database[id].second.third);
}

//////////////////////////////////////////////////////////////////////////////

InterpolatorID InterpolatorRegister::getInterpolatorID(const std::string& name) const
  throw (Common::NoSuchValueException)
{
  database_type::const_iterator begin = _database.begin();
  database_type::const_iterator end   = _database.end();
  database_type::const_iterator itr =
    find_if(begin,end,EqualName(name));

  if(itr != end) {
   return distance(begin,itr);
  }
  else {
    throw Common::NoSuchValueException (FromHere(),"Interpolator: " + name + " not found");
  }
}

//////////////////////////////////////////////////////////////////////////////

InterpolatorID InterpolatorRegister::getInterpolatorID(
                                 const CFPolyForm::Type&  interpolType,
                                 const CFPolyOrder::Type& interpolOrder,
                                 const CFGeoShape::Type&  shape) const
                                 throw (Common::NoSuchValueException)
{
  database_type::const_iterator begin = _database.begin();
  database_type::const_iterator end   = _database.end();
  database_type::const_iterator itr =
    find_if(begin,end,EqualProperties(interpolType,
                                      interpolOrder,
                                      shape));

  if(itr != end) {
   return distance(begin,itr);
  }
  else {
    throw Common::NoSuchValueException (FromHere(),"Interpolator with given properties not found.");
  }
}

//////////////////////////////////////////////////////////////////////////////

InterpolatorID InterpolatorRegister::registInterpolator(
                                 const std::string&    name,
                                 const CFPolyForm::Type&  interpolType,
                                 const CFPolyOrder::Type& interpolOrder,
                                 const CFGeoShape::Type&  shape)
{
  InterpolatorID id = 0;

  database_type::iterator begin = _database.begin();
  database_type::iterator end   = _database.end();
  database_type::iterator itr =
    find_if(begin,end,EqualName(name));

  if(itr != end) {
   id = distance(begin,itr);
  }
  else {
    entry_type newEntry;
    newEntry.first = name;
    newEntry.second.first  = interpolType;
    newEntry.second.second = interpolOrder;
    newEntry.second.third  = shape;

    _database.push_back(newEntry);
    id = _database.size() - 1;
  }
  return id;
}

//////////////////////////////////////////////////////////////////////////////


  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
