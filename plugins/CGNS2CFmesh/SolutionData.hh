// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS_SolutionData_hh
#define COOLFluiD_CGNS_SolutionData_hh

//////////////////////////////////////////////////////////////////////////////

#include "CGNS2CFmesh/FieldData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CGNS {

//////////////////////////////////////////////////////////////////////////////

/// Storage of CGNS field data
class SolutionData {

public: // functions

  /// constructor
  SolutionData() {}

  /// destructor
  ~SolutionData() {}

public: // data

    /// name of the solution
    char solname[CGNS_CHAR_MAX];
    /// grid location of the solution
    CGNSLIB::GridLocation_t loc;
    /// number of fields
    int m_nfields;
    /// data storage for the fields
    std::vector< FieldData > m_fields;

}; // end class SolutionData

//////////////////////////////////////////////////////////////////////////////

  } // namespace CGNS
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CGNS_SolutionData_hh
