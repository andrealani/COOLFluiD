// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/assign/list_of.hpp> // for map_list_of

#include "Framework/CFQuadrature.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

CFQuadrature::Convert::FwdMap_t CFQuadrature::Convert::all_fwd = boost::assign::map_list_of
    ( CFQuadrature::INVALID,         "INVALID")
    ( CFQuadrature::GAUSSLEGENDRE,   "GaussLegendre")
    ( CFQuadrature::DGGAUSSLEGENDRE, "DGGaussLegendre")
    ( CFQuadrature::GAUSSLOBATO,     "GaussLobato")
    ( CFQuadrature::NEWTONCOTES,     "NewtonCotes")
    ( CFQuadrature::SIMPSON,         "Simpson")
    ( CFQuadrature::TRAPEZIUM,       "Trapezium");

CFQuadrature::Convert::BwdMap_t CFQuadrature::Convert::all_rev = boost::assign::map_list_of
    ( "INVALID", CFQuadrature::INVALID )
    ( "GaussLegendre", CFQuadrature::GAUSSLEGENDRE )
    ( "DGGaussLegendre", CFQuadrature::DGGAUSSLEGENDRE )
    ( "GaussLobato", CFQuadrature::GAUSSLOBATO )
    ( "NewtonCotes", CFQuadrature::NEWTONCOTES )
    ( "Simpson", CFQuadrature::SIMPSON )
    ( "Trapezium", CFQuadrature::TRAPEZIUM );

//////////////////////////////////////////////////////////////////////////////

std::ostream& operator<< ( std::ostream& os, const CFQuadrature::Type& in )
{
  os << CFQuadrature::Convert::to_str(in);
  return os;
}

std::istream& operator>> (std::istream& is, CFQuadrature::Type& in )
{
  std::string tmp;
  is >> tmp;
  in = CFQuadrature::Convert::to_enum(tmp);
  return is;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD
