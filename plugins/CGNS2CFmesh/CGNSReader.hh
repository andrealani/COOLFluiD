// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CGNS_CGNSReader_hh
#define COOLFluiD_CGNS_CGNSReader_hh

//////////////////////////////////////////////////////////////////////////////




#include "CGNS2CFmesh/CGNSData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CGNS {

//////////////////////////////////////////////////////////////////////////////

/// Reader of CGNS files
class CGNSReader
{
  public: // functions

    /// constructor
    CGNSReader();
    /// destructor
    ~CGNSReader();
    /// reads a CGNS file
    void read(const std::string& filename);
    /// sets the storage where the reader should put the data
    /// that will be read from the file
    /// @pre pointer is not CFNULL
    void setCGNSData(CGNSData * data);

  protected: // helper functions

    /// handles the return error code from a CGNS call
    /// @throw CGNSException if an error occured
    void handle_CGNS_error(int ierr);

    /// reads data from a base in the CGNS file
    /// @throw CGNSException if an error occured
    void readBase();

    /// reads data from a zone in the CGNS file
    /// @throw CGNSException if an error occured
    void readZone();

    /// reads data from a grid in the CGNS file
    /// @throw CGNSException if an error occured
    void readGrid();

    /// reads data from a solution in the CGNS file
    /// @throw CGNSException if an error occured
    void readSolution();

    /// reads data from a section in the CGNS file
    /// @throw CGNSException if an error occured
    void readSection();

    /// reads data from a field in the CGNS file
    /// @throw CGNSException if an error occured
    void readField();

    /// reads data from a boundary confition in the CGNS file
    /// @throw CGNSException if an error occured
    void readBoco();

  private: // data

    /// pointer to the storage of the CGNSData
    /// where this reader will place the data read
    CGNSData * m_data;
    /// index of the file that is open
    int index_file;
    /// number of bases in file
    int m_nbbases;

    /// iterator to the current base
    std::vector< BaseData >::iterator bitr;
    /// index of the current base
    int index_base;
    /// iterator to the current base
    std::vector< ZoneData >::iterator zitr;
    /// index of the zone
    int index_zone;
    /// iterator to the current grid
    std::vector< GridData >::iterator gitr;

    /// iterator to the current solution
    std::vector< SolutionData >::iterator sitr;
    /// index of the solution
    int index_sol;
    /// iterator to the current field
    std::vector< FieldData >::iterator fitr;
    /// iterator to the current section
    std::vector< SectionData >::iterator secitr;
    /// index of the section
    int index_sect;
    /// iterator to the current boundary condition
    std::vector< BocoData >::iterator bcitr;
    /// index of the boundary condition
    int index_boco;

}; // end class CGNSReader

//////////////////////////////////////////////////////////////////////////////

  } // namespace CGNS
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CGNS_CGNSReader_hh
