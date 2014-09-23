// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS_GridData_hh
#define COOLFluiD_CGNS_GridData_hh

//////////////////////////////////////////////////////////////////////////////


#include "CGNS2CFmesh/CGNSDefinitions.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CGNS {

//////////////////////////////////////////////////////////////////////////////

/// Storage of CGNS field data
class GridData
{
  public: // functions
    /// constructor
    GridData() : x(NULL), y(NULL), z(NULL) {}
    /// destructor
    ~GridData() { deallocate(); }
    /// helper function to allocate the solution data
    void allocate_x( int in_size )
    {
      cf_assert(x == NULL);
      size = in_size;
      x = new double[size];
    }
    /// helper function to allocate the solution data
    void allocate_y( int in_size )
    {
      cf_assert(y == NULL);
      size = in_size;
      y = new double[size];
    }
    /// helper function to allocate the solution data
    void allocate_z( int in_size )
    {
      cf_assert(z == NULL);
      size = in_size;
      z = new double[size];
    }
    /// helper function to deallocate the solution data
    void deallocate()
    {
      deletePtrArray(x);
      deletePtrArray(y);
      deletePtrArray(z);
      size = 0;
    }
    public: // data
    /// cordinates of nodes in X
    double * x;
    /// cordinates of nodes in Y
    double * y;
    /// cordinates of nodes in Z
    double * z;
    /// size of the coordinates vector
    unsigned int size;

}; // end class GridData

//////////////////////////////////////////////////////////////////////////////

  } // namespace CGNS
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CGNS_GridData_hh
