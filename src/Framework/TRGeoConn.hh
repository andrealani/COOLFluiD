// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_TRGeoConn_hh
#define COOLFluiD_Framework_TRGeoConn_hh

//////////////////////////////////////////////////////////////////////////////

#include <valarray>
#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

typedef std::valarray<CFuint> GeoConnElementPart;

typedef std::pair<GeoConnElementPart,
            GeoConnElementPart> GeoConnElement;

/// GeometricEntity's Connectivity in a given TR
typedef std::vector<GeoConnElement> GeoConn;

/// GeometricEntity's Connectivity in a given TRS
typedef std::vector< GeoConn > TRGeoConn;

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_TRGeoConn_hh
