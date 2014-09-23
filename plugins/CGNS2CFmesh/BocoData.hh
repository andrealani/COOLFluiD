// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS_BocoData_hh
#define COOLFluiD_CGNS_BocoData_hh

//////////////////////////////////////////////////////////////////////////////

#include "CGNS2CFmesh/CGNSDefinitions.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CGNS {

//////////////////////////////////////////////////////////////////////////////

/// Storage of CGNS section data
class BocoData {
public: // functions

  /// constructor
  BocoData() : ibcpnts(CFNULL) {}
  /// destructor
  ~BocoData() { deallocate(); }
  /// helper function to allocate the solution data
  /// @pre must set first nbcpts
  void allocate()
  {
    cf_assert(nbcpts > 0);
    cf_assert(ibcpnts == CFNULL);
    ibcpnts = new int[nbcpts];
  }

  /// helper function to deallocate the solution data
  void deallocate()
  {
    deletePtr(ibcpnts);
    nbcpts = 0;
  }

public: // data
  /// name for boundary condition
  char boconame[CGNS_CHAR_MAX];
  /// type of the boundary condition
  CGNSLIB::BCType_t ibocotype;
  /// type of the point set
  CGNSLIB::PointSetType_t iptset;
  /// number of points in bc
  int nbcpts;
  /// points on the boundary condition
  int * ibcpnts;
  /// data type used for definition of the normals
  /// either RealSingle or RealDouble
  CGNSLIB::DataType_t normaldatatype;
  ///  data relevant to the normals that we dont use
  int normalindex[3];
  ///  data relevant to the normals that we dont use
  int normallistflag;
  ///  data relevant to the normals that we dont use
  int ndataset;
  ///  data relevant to the normals that we dont use
  int normallist;

}; // end class BocoData

//////////////////////////////////////////////////////////////////////////////

  } // namespace CGNS
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CGNS_BocoData_hh
