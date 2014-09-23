// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iostream>

#include "Common/CFLog.hh"
#include "CGNS2CFmesh/CGNSException.hh"
#include "CGNS2CFmesh/CGNSReader.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace CGNS {

using namespace CGNSLIB;

//////////////////////////////////////////////////////////////////////////////

#define CALL_CGNS(x) handle_CGNS_error ( x )

//////////////////////////////////////////////////////////////////////////////

CGNSReader::CGNSReader()
{
}

//////////////////////////////////////////////////////////////////////////////

CGNSReader::~CGNSReader()
{
}

//////////////////////////////////////////////////////////////////////////////

void CGNSReader::setCGNSData(CGNSData * data)
{
  cf_assert(data != CFNULL);
  m_data = data;
}

//////////////////////////////////////////////////////////////////////////////

void CGNSReader::read(const std::string& filename)
{
  CFAUTOTRACE;

  // open file in read mode
  CALL_CGNS ( cg_open(filename.c_str(),MODE_READ,&index_file) );

  // check how many bases we have
  CALL_CGNS ( cg_nbases(index_file,&m_nbbases) );
  m_data->m_bases.resize(m_nbbases);

  // for now handle a single base
  if (m_nbbases != 1) throw CGNSException (FromHere(),"More than one base detected in file. Dont know how to handle multi-bases.");

  bitr = m_data->m_bases.begin();
  index_base = 1;
  for (; bitr != m_data->m_bases.end(); ++bitr, ++index_base)
  {
    readBase();
  }

  // close the CGNS file
  CALL_CGNS ( cg_close(index_file) );
}

//////////////////////////////////////////////////////////////////////////////

void CGNSReader::readBase()
{
  CFAUTOTRACE;

  // get the dimensions
  CALL_CGNS ( cg_base_read(index_file,index_base,bitr->basename,&bitr->m_cell_dim,&bitr->m_phys_dim) );
  if (bitr->m_cell_dim != bitr->m_phys_dim) throw CGNSException (FromHere(),"Cannot handle cell dimensions different than physical dimensions.");

  CFLog (INFO,"base cell dim : " << bitr->m_cell_dim  << "\n");
  CFLog (INFO,"base phys dim : " << bitr->m_phys_dim  << "\n");

  // check how many zones we have
  CALL_CGNS ( cg_nzones(index_file,index_base,&bitr->m_nbzones) );
  bitr->m_zones.resize(bitr->m_nbzones);

  // handle only single zone files for now
  if (bitr->m_nbzones != 1) throw CGNSException (FromHere(),"More than one zone detected in file. Dont know how to handle multi-zones.");

  zitr = bitr->m_zones.begin();
  index_zone = 1;
  for (; zitr != bitr->m_zones.end(); ++zitr, ++index_zone)
  {
    readZone();
  }
}

//////////////////////////////////////////////////////////////////////////////

void CGNSReader::readZone()
{
  CFAUTOTRACE;

  // get zone type
  CALL_CGNS ( cg_zone_type(index_file,index_base,index_zone,&zitr->ztype) );
  if (zitr->ztype == Structured) throw CGNSException (FromHere(),"Cannot handle structured meshes.");

  // get zone size and name
  CALL_CGNS ( cg_zone_read(index_file,index_base,index_zone,zitr->zonename,zitr->isize[0]) );
  // get the number of grids
  CALL_CGNS ( cg_ngrids(index_file,index_base,index_zone,&zitr->m_nbgrids) );
  // find out number of solutions
  CALL_CGNS ( cg_nsols(index_file,index_base,index_zone,&zitr->m_nsols) );
  // find out how many sections
  CALL_CGNS ( cg_nsections(index_file,index_base,index_zone,&zitr->m_nsections) );
  // find out number of BCs that exist under this zone
  CALL_CGNS ( cg_nbocos(index_file,index_base,index_zone,&zitr->nbocos) );

  CFLog (INFO, "zone name   : " << zitr->zonename << "\n");
  CFLog (INFO, "nb nodes    : " << zitr->isize[CGNS_VERT_IDX][0]    << "\n");
  CFLog (INFO, "nb cells    : " << zitr->isize[CGNS_CELL_IDX][0]    << "\n");
  CFLog (INFO, "nb bnodes   : " << zitr->isize[CGNS_BVRT_IDX][0]    << "\n");
  CFLog (INFO, "nb grids    : " << zitr->m_nbgrids   << "\n");
  CFLog (INFO, "nb sols     : " << zitr->m_nsols     << "\n");
  CFLog (INFO, "nb sections : " << zitr->m_nsections << "\n");
  CFLog (INFO, "nb bcs      : " << zitr->nbocos << "\n");

  // only handle one grid per cgns file
  if (zitr->m_nbgrids != 1) throw CGNSException (FromHere(),"Cannot handle nore than one grid per zone.");
  // only handle one solution per cgns file
  if(zitr->m_nsols > 1) throw CGNSException (FromHere(),"Cannot handle more than one solution per zone");
  // there must be more at least one section in the file, containing the element connectivity
  if(zitr->m_nsections < 1) throw CGNSException (FromHere(),"Could not find any section in file.");

  // lower range index
  zitr->irmin=1;
  // upper range index of vertices
  zitr->irmax=zitr->isize[CGNS_VERT_IDX][0];

  zitr->m_grids.resize(bitr->m_nbzones);

  gitr = zitr->m_grids.begin();
  for (; gitr != zitr->m_grids.end(); ++gitr)
  {
    /// @warning the readGrid() function only supports reading the primary grid
    readGrid();
  }

  zitr->m_solutions.resize(zitr->m_nsols);
  sitr = zitr->m_solutions.begin();
  index_sol = 1;
  for (; sitr != zitr->m_solutions.end(); ++sitr, ++index_sol)
  {
    readSolution();
  }

  // read element connectivity
  zitr->m_sections.resize(zitr->m_nsections);
  secitr = zitr->m_sections.begin();
  index_sect = 1;
  for (; secitr != zitr->m_sections.end(); ++secitr, ++index_sect)
  {
    readSection();
  }

  // resize storage of bcs
  zitr->m_bocos.resize(zitr->nbocos);
  bcitr = zitr->m_bocos.begin();
  index_boco = 1;
  for (; bcitr != zitr->m_bocos.end(); ++bcitr, ++index_boco)
  {
    readBoco();
  }
}

//////////////////////////////////////////////////////////////////////////////

void CGNSReader::readGrid()
{
  CFAUTOTRACE;

  // alsitr->locate the grid storage
  gitr->allocate_x(zitr->irmax);
  gitr->allocate_y(zitr->irmax);
  if (bitr->m_phys_dim == DIM_3D)
    gitr->allocate_z(zitr->irmax);

  // read grid coordinates
  CALL_CGNS ( cg_coord_read(index_file,index_base,index_zone,"CoordinateX", RealDouble,&zitr->irmin,&zitr->irmax,gitr->x) );
  CALL_CGNS ( cg_coord_read(index_file,index_base,index_zone,"CoordinateY", RealDouble,&zitr->irmin,&zitr->irmax,gitr->y) );

  if (bitr->m_phys_dim == DIM_3D)
    CALL_CGNS ( cg_coord_read(index_file,index_base,index_zone,"CoordinateZ", RealDouble,&zitr->irmin,&zitr->irmax,gitr->z) );
}

//////////////////////////////////////////////////////////////////////////////

void CGNSReader::readSolution()
{
  CFAUTOTRACE;

    // get information about solution
    CALL_CGNS ( cg_sol_info(index_file,index_base,index_zone,index_sol,sitr->solname,&sitr->loc) );
    CFLog (INFO,"solution name : " << sitr->solname << "\n");

    // we dont handle yet solutions other than vertex based
    if (sitr->loc != Vertex)
      throw CGNSException (FromHere(),"Cannot handle sitr->location of solution different than vertex based");

    CALL_CGNS ( cg_nfields(index_file,index_base,index_zone,index_sol,&sitr->m_nfields) );
    CFLog (INFO,"nb fields : " << sitr->m_nfields << "\n");

    // resize the storage for the solution
    sitr->m_fields.resize(sitr->m_nfields);

    int ifield = 1;
    for (fitr = sitr->m_fields.begin(); fitr != sitr->m_fields.end(); ++fitr, ++ifield )
    {
      // alsitr->locate the storage for the solution of this field
      fitr->allocate(zitr->irmax - zitr->irmin+1);
      // get the info about this field
      CALL_CGNS ( cg_field_info(index_file,index_base,index_zone,index_sol, ifield, &fitr->fielddatatype, fitr->fieldname) );
      CFLog (INFO,"field id : "   << ifield << "\n");
      CFLog (INFO,"field name : " << fitr->fieldname << "\n");
      // read the data of this field
      CALL_CGNS ( cg_field_read(index_file,index_base,index_zone,index_sol, fitr->fieldname, RealDouble,&zitr->irmin,&zitr->irmax,fitr->field) );
    }
}

//////////////////////////////////////////////////////////////////////////////

void CGNSReader::readSection()
{
  CFAUTOTRACE;

  // unused data
  int nbndry,iparent_flag,iparentdata;

  // reading section data
  CALL_CGNS ( cg_section_read(index_file,index_base,index_zone,index_sect,secitr->sectionname, &secitr->itype,&secitr->istart,&secitr->iend,&nbndry,&iparent_flag ) );

  CFLog(INFO,"section name : " << secitr->sectionname << "\n");
  CFLog(INFO,"section type : " << ElementTypeName[secitr->itype] << "\n");
  CFLog(INFO,"secitr->istart,secitr->iend  : " << secitr->istart << " " << secitr->iend << "\n");

  // get the number of nodes per element
  CALL_CGNS ( cg_npe(secitr->itype,&secitr->nbnodes_per_elem) );

  // allocate the connecivity storage
  secitr->allocate();

  // get the connecivity for this element type
  CALL_CGNS ( cg_elements_read(index_file,index_base,index_zone,index_sect,secitr->ielem, &iparentdata) );
}

//////////////////////////////////////////////////////////////////////////////

void CGNSReader::readBoco()
{
  CFAUTOTRACE;

  // get BC info
  CALL_CGNS ( cg_boco_info(index_file,index_base,index_zone,index_boco,bcitr->boconame,&bcitr->ibocotype,&bcitr->iptset,&bcitr->nbcpts,bcitr->normalindex,&bcitr->normallistflag,&bcitr->normaldatatype,&bcitr->ndataset) );

  // check for types of BCs that we can handle
  if (bcitr->iptset != PointList) throw CGNSException (FromHere(),"Cannot handle boundary conditions different than PointList");

  CFLog(INFO, "BC number : " << index_boco << "\n");
  CFLog(INFO, "BC name   : " << bcitr->boconame << "\n");
  CFLog(INFO, "BC type   : " << BCTypeName[bcitr->ibocotype] << "\n");
  CFLog(INFO, "BC nb elms: " << bcitr->nbcpts << "\n");

  // allocate the storage for the point list
  bcitr->allocate();

  // read point list
  CALL_CGNS( cg_boco_read(index_file,index_base,index_zone,index_boco,bcitr->ibcpnts,&bcitr->normallist) );
}

//////////////////////////////////////////////////////////////////////////////

void CGNSReader::handle_CGNS_error(int ierr)
{
  if (ierr)
  {
    const char * error_msg = cg_get_error();
    throw CGNSException (FromHere(),error_msg);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace CGNS
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////


