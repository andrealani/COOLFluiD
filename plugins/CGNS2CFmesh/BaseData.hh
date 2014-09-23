// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS_BaseData_hh
#define COOLFluiD_CGNS_BaseData_hh

//////////////////////////////////////////////////////////////////////////////

#include "CGNS2CFmesh/ZoneData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CGNS {

//////////////////////////////////////////////////////////////////////////////

/// Storage of CGNS field data
class BaseData {

public: // functions

  /// constructor
  BaseData() {}

  /// destructor
  ~BaseData() {}

public: // data

  /// name of the base
  char basename[CGNS_CHAR_MAX];
  /// dimensionality of the cells
  /// could be smaller than the physical dimensions
  int m_cell_dim;
  /// dimensionality of the physical domain
  /// is the number  of coordinates necessary to define a vector in the domain
  int m_phys_dim;

  /// total number of zones in file
  int m_nbzones;
  /// vector of zones
  std::vector< ZoneData > m_zones;

}; // end class BaseData

//////////////////////////////////////////////////////////////////////////////

  } // namespace CGNS
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CGNS_BaseData_hh
