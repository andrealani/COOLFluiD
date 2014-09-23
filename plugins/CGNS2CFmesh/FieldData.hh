// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS_FieldData_hh
#define COOLFluiD_CGNS_FieldData_hh

//////////////////////////////////////////////////////////////////////////////


#include "CGNS2CFmesh/CGNSDefinitions.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CGNS {

//////////////////////////////////////////////////////////////////////////////

/// Storage of CGNS field data
class FieldData
{
  public: // functions
    /// constructor
    FieldData() : field(CFNULL) {}
    /// destructor
    ~FieldData() { deallocate(); }
    /// helper function to allocate the solution data
    void allocate( int in_size)
    {
      cf_assert(field == CFNULL);
      size = in_size;
      field = new double[size];
    }
    /// helper function to deallocate the solution data
    void deallocate()
    {
      deletePtr(field);
      size = 0;
    }
    public: // data
      /// name of the field
      char fieldname[CGNS_CHAR_MAX];
      /// data type of this field
      CGNSLIB::DataType_t fielddatatype;
      /// storage of this solution field
      double * field;
      /// size of the solution storage
      unsigned int size;
}; // end class FieldData

//////////////////////////////////////////////////////////////////////////////

  } // namespace CGNS
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CGNS_FieldData_hh
