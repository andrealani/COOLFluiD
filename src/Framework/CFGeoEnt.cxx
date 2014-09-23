// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/assign/list_of.hpp> // for map_list_of

#include "Framework/CFGeoEnt.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

CFGeoEnt::Convert::FwdMap_t CFGeoEnt::Convert::all_fwd = boost::assign::map_list_of
    ( CFGeoEnt::INVALID, "INVALID" )
    ( CFGeoEnt::EDGE,    "Edge")
    ( CFGeoEnt::FACE,    "Face")
    ( CFGeoEnt::CELL,    "Cell");

CFGeoEnt::Convert::BwdMap_t CFGeoEnt::Convert::all_rev = boost::assign::map_list_of
    ( "INVALID", CFGeoEnt::INVALID )
    ( "Edge",    CFGeoEnt::EDGE)
    ( "Face",    CFGeoEnt::FACE)
    ( "Cell",    CFGeoEnt::CELL);

//////////////////////////////////////////////////////////////////////////////

std::ostream& operator<< ( std::ostream& os, const CFGeoEnt::Type& in )
{
  os << CFGeoEnt::Convert::to_str(in);
  return os;
}

//////////////////////////////////////////////////////////////////////////////

std::istream& operator>> (std::istream& is, CFGeoEnt::Type& in )
{
  std::string tmp;
  is >> tmp;
  in = CFGeoEnt::Convert::to_enum(tmp);
  return is;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD
