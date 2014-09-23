// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS_CGNSData_hh
#define COOLFluiD_CGNS_CGNSData_hh

//////////////////////////////////////////////////////////////////////////////

#include "CGNS2CFmesh/BaseData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CGNS {

//////////////////////////////////////////////////////////////////////////////

/// Storage of CGNS field data
class CGNSData {

public: // functions

  /// constructor
  CGNSData() {}

  /// destructor
  ~CGNSData() {}

public: // data

  /// vector of bases
  std::vector< BaseData > m_bases;

}; // end class CGNSData

//////////////////////////////////////////////////////////////////////////////

  } // namespace CGNS
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CGNS_CGNSData_hh
