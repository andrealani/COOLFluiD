// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/assign/list_of.hpp> // for map_list_of

#include "Framework/CFIntegration.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

CFIntegration::Convert::FwdMap_t CFIntegration::Convert::all_fwd = boost::assign::map_list_of
    ( CFIntegration::INVALID, "INVALID")
    ( CFIntegration::CONTOUR, "CONTOUR")
    ( CFIntegration::VOLUME,  "VOLUME");

CFIntegration::Convert::BwdMap_t CFIntegration::Convert::all_rev = boost::assign::map_list_of
    ( "INVALID", CFIntegration::INVALID )
    ( "CONTOUR", CFIntegration::CONTOUR )
    ( "VOLUME",  CFIntegration::VOLUME );

//////////////////////////////////////////////////////////////////////////////

std::ostream& operator<< ( std::ostream& os, const CFIntegration::Type& in )
{
  os << CFIntegration::Convert::to_str(in);
  return os;
}

std::istream& operator>> (std::istream& is, CFIntegration::Type& in )
{
  std::string tmp;
  is >> tmp;
  in = CFIntegration::Convert::to_enum(tmp);
  return is;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD
