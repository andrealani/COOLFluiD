// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/assign/list_of.hpp> // for map_list_of

#include "Framework/CFPolyForm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

CFPolyForm::Convert::FwdMap_t CFPolyForm::Convert::all_fwd = boost::assign::map_list_of
    ( CFPolyForm::INVALID,            "INVALID" )
    ( CFPolyForm::LAGRANGE,           "Lagrange" )
    ( CFPolyForm::HERMITE,            "Hermite" )
    ( CFPolyForm::SERENDIPITY,        "Serendipity" )
    ( CFPolyForm::SPECTRALFV,         "SpectralFV" )
    ( CFPolyForm::LEGENDRE,           "Legendre" )
    ( CFPolyForm::SPECTRALFD,         "SpectralFD" )
    ( CFPolyForm::FLUXRECONSTRUCTION, "FluxReconstruction");


CFPolyForm::Convert::BwdMap_t CFPolyForm::Convert::all_rev = boost::assign::map_list_of
    ( "INVALID",            CFPolyForm::INVALID )
    ( "Lagrange",           CFPolyForm::LAGRANGE )
    ( "Hermite",            CFPolyForm::HERMITE )
    ( "Serendipity",        CFPolyForm::SERENDIPITY )
    ( "SpectralFV",         CFPolyForm::SPECTRALFV )
    ( "Legendre",           CFPolyForm::LEGENDRE )
    ( "SpectralFD",         CFPolyForm::SPECTRALFD )
    ( "FluxReconstruction", CFPolyForm::FLUXRECONSTRUCTION);

//////////////////////////////////////////////////////////////////////////////

std::ostream& operator<< ( std::ostream& os, const CFPolyForm::Type& in )
{
  os << CFPolyForm::Convert::to_str(in);
  return os;
}

std::istream& operator>> (std::istream& is, CFPolyForm::Type& in )
{
  std::string tmp;
  is >> tmp;
  in = CFPolyForm::Convert::to_enum(tmp);
  return is;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
