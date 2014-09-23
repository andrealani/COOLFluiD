// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS_ZoneData_hh
#define COOLFluiD_CGNS_ZoneData_hh

//////////////////////////////////////////////////////////////////////////////

#include "CGNS2CFmesh/GridData.hh"
#include "CGNS2CFmesh/SolutionData.hh"
#include "CGNS2CFmesh/BocoData.hh"
#include "CGNS2CFmesh/SectionData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CGNS {

//////////////////////////////////////////////////////////////////////////////

/// Storage of CGNS field data
class ZoneData {

public: // functions

  /// constructor
  ZoneData() {}

  /// destructor
  ~ZoneData() {}

  /// @return the number of vertices in this zone
  int getNbVertices() { return isize[CGNS_VERT_IDX][0]; };

  /// @return the number of elements in this zone
  int getNbElems() { return isize[CGNS_CELL_IDX][0]; };

public: // data

  /// name of the zone
  char zonename[CGNS_CHAR_MAX];
  /// there is only one index dimension, and for that one we store
  /// the number of vertices, the number of cells and the number of boundary vertices
  int isize[3][1];
  /// lower range index
  int irmin;
  /// upper range index
  int irmax;
  /// the type of the Zone
  CGNSLIB::ZoneType_t ztype;
  /// number of grids in zone
  int m_nbgrids;
  /// vector of grids
  std::vector< GridData > m_grids;
  /// number of solutions per zone
  int m_nsols;
  /// vector of solutions
  std::vector< SolutionData > m_solutions;
  /// number of section
  int m_nsections;
  /// vector of sections
  std::vector< SectionData > m_sections;
  /// number of boundary consitions
  int nbocos;
  /// vector of boundary conditions data
  std::vector< BocoData > m_bocos;

}; // end class ZoneData

//////////////////////////////////////////////////////////////////////////////

  } // namespace CGNS
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CGNS_ZoneData_hh
