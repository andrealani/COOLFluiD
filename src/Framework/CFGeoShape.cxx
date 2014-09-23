// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/assign/list_of.hpp> // for map_list_of

#include "Framework/CFGeoShape.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

CFGeoShape::Convert::FwdMap_t CFGeoShape::Convert::all_fwd = boost::assign::map_list_of
    ( CFGeoShape::INVALID, "INVALID" )
    ( CFGeoShape::POINT,   "Point")
    ( CFGeoShape::LINE,    "Line")
    ( CFGeoShape::TRIAG,   "Triag")
    ( CFGeoShape::QUAD,    "Quad")
    ( CFGeoShape::TETRA,   "Tetra")
    ( CFGeoShape::PYRAM,   "Pyram")
    ( CFGeoShape::PRISM,   "Prism")
    ( CFGeoShape::HEXA,    "Hexa");

CFGeoShape::Convert::BwdMap_t CFGeoShape::Convert::all_rev = boost::assign::map_list_of
    ("INVALID",  CFGeoShape::INVALID)
    ("Point",    CFGeoShape::POINT)
    ("Line",     CFGeoShape::LINE)
    ("Triag",    CFGeoShape::TRIAG)
    ("Quad",     CFGeoShape::QUAD)
    ("Tetra",    CFGeoShape::TETRA)
    ("Pyram",    CFGeoShape::PYRAM)
    ("Prism",    CFGeoShape::PRISM)
    ("Hexa",     CFGeoShape::HEXA);

//////////////////////////////////////////////////////////////////////////////

std::ostream& operator<< ( std::ostream& os, const CFGeoShape::Type& in )
{
  os << CFGeoShape::Convert::to_str(in);
  return os;
}

std::istream& operator>> (std::istream& is, CFGeoShape::Type& in )
{
  std::string tmp;
  is >> tmp;
  in = CFGeoShape::Convert::to_enum(tmp);
  return is;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
