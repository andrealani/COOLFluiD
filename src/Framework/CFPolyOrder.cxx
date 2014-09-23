// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/assign/list_of.hpp> // for map_list_of

#include "Framework/CFPolyOrder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

CFPolyOrder::Convert::FwdMap_t CFPolyOrder::Convert::all_fwd = boost::assign::map_list_of
    ( CFPolyOrder::INVALID,  "INVALID" )
    ( CFPolyOrder::ORDER0,   "P0" )
    ( CFPolyOrder::ORDER1,   "P1" )
    ( CFPolyOrder::ORDER2,   "P2" )
    ( CFPolyOrder::ORDER3,   "P3" )
    ( CFPolyOrder::ORDER4,   "P4" )
    ( CFPolyOrder::ORDER5,   "P5" )
    ( CFPolyOrder::ORDER6,   "P6" )
    ( CFPolyOrder::ORDER7,   "P7" )
    ( CFPolyOrder::ORDER8,   "P8" )
    ( CFPolyOrder::ORDER9,   "P9" )
    ( CFPolyOrder::ORDER10,  "P10" )
    ( CFPolyOrder::MAXORDER, "CFPolyOrder::MAXORDER" );


CFPolyOrder::Convert::BwdMap_t CFPolyOrder::Convert::all_rev = boost::assign::map_list_of
    ("INVALID" , CFPolyOrder::INVALID  )
    ("P0",  CFPolyOrder::ORDER0)
    ("P1",  CFPolyOrder::ORDER1)
    ("P2",  CFPolyOrder::ORDER2)
    ("P3",  CFPolyOrder::ORDER3)
    ("P4",  CFPolyOrder::ORDER4)
    ("P5",  CFPolyOrder::ORDER5)
    ("P6",  CFPolyOrder::ORDER6)
    ("P7",  CFPolyOrder::ORDER7)
    ("P8",  CFPolyOrder::ORDER8)
    ("P9",  CFPolyOrder::ORDER9)
    ("P10", CFPolyOrder::ORDER10)
    ("CFPolyOrder::MAXORDER", CFPolyOrder::MAXORDER);

//////////////////////////////////////////////////////////////////////////////

std::ostream& operator<< ( std::ostream& os, const CFPolyOrder::Type& in )
{
  os << CFPolyOrder::Convert::to_str(in);
  return os;
}

std::istream& operator>> (std::istream& is, CFPolyOrder::Type& in )
{
  std::string tmp;
  is >> tmp;
  in = CFPolyOrder::Convert::to_enum(tmp);
  return is;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
