// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS_SectionData_hh
#define COOLFluiD_CGNS_SectionData_hh

//////////////////////////////////////////////////////////////////////////////


#include "CGNS2CFmesh/CGNSDefinitions.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CGNS {

//////////////////////////////////////////////////////////////////////////////

/// Storage of CGNS section data
class SectionData
{
  public: // functions
    /// constructor
    SectionData() : ielem(CFNULL) {}
    /// destructor
    ~SectionData() { deallocate(); }
    /// helper function to allocate the solution data
    /// @pre must set first istart, iend and nbnodes_per_elem
    void allocate()
    {
      cf_assert(ielem == CFNULL);
      nbelems = iend-istart+1;
      ielem = new int[nbelems * nbnodes_per_elem];
    }
    /// helper function to deallocate the solution data
    void deallocate()
    {
      deletePtr(ielem);
      nbelems = 0;
    }
    public: // data
    /// name of the section
    char sectionname[CGNS_CHAR_MAX];
    /// element to node connectivity
    /// sized with number of cells in section then number of nodes
    int * ielem;
    /// number of elements in this section
    int nbelems;
    /// number of nodes per element
    int nbnodes_per_elem;
    /// begin index
    int istart;
    /// end index
    int iend;
    /// element type
    CGNSLIB::ElementType_t itype;

}; // end class SectionData

//////////////////////////////////////////////////////////////////////////////

  } // namespace CGNS
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CGNS_SectionData_hh
