#include "Environment/ObjectProvider.hh"

#include "UFEM/TetraP1P1Cell/CellProps.hh"
#include "UFEM/TetraP1P1Cell/UFEMTetraP1P1Cell.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TetraP1P1Cell {

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider < CellProps,
                              ElemProps,
                              UFEMTetraP1P1CellPlugin,
                              ElemProps::NARGS >
aTetraP1P1CellProps_Provider ( "TetraP1P1CellProps" );

//////////////////////////////////////////////////////////////////////////////

void CellProps::defineConfigOptions ( Config::OptionList& options )
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

CellProps::CellProps ( const std::string& name ) : ElemProps ( name )
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

CellProps::~CellProps() 
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void CellProps::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  ElemProps::configure ( args );
}

//////////////////////////////////////////////////////////////////////////////

void CellProps::setup ()
{
  CFAUTOTRACE;
  ElemProps::setup ();
}

//////////////////////////////////////////////////////////////////////////////

void CellProps::unsetup ()
{
  CFAUTOTRACE;
  ElemProps::unsetup ();
}

//////////////////////////////////////////////////////////////////////////////

void CellProps::prepare ( const Framework::GeometricEntity& cell )
{
  CFAUTOTRACE;
  // nodal and elem coordinates
  const std::vector<Node*>& nodes = cell.getNodes();
  Node& node0 = *nodes[0];
  Node& node1 = *nodes[1];
  Node& node2 = *nodes[2];
  Node& node3 = *nodes[3];

    // volume of the elem
  CFreal d[3][3];
  
    d[0][0] = node1[XX] - node0[XX];
    d[0][1] = node1[YY] - node0[YY];
    d[0][2] = node1[ZZ] - node0[ZZ];
    d[1][0] = node2[XX] - node0[XX];
    d[1][1] = node2[YY] - node0[YY];
    d[1][2] = node2[ZZ] - node0[ZZ];
    d[2][0] = node3[XX] - node0[XX];
    d[2][1] = node3[YY] - node0[YY];
    d[2][2] = node3[ZZ] - node0[ZZ];

  m_celldata.vol =(  d[0][0]*(d[1][1]*d[2][2]-d[2][1]*d[1][2])
                    +d[1][0]*(d[2][1]*d[0][2]-d[0][1]*d[2][2])
                    +d[2][0]*(d[0][1]*d[1][2]-d[1][1]*d[0][2]) ) / 6.0;

  
  // normals of elem faces
  m_celldata.nx[0] = ((node3[YY]-node1[YY])*(node2[ZZ]-node1[ZZ]) - (node2[YY]-node1[YY])*(node3[ZZ]-node1[ZZ])) / 2.0;
  m_celldata.ny[0] = ((node3[ZZ]-node1[ZZ])*(node2[XX]-node1[XX]) - (node2[ZZ]-node1[ZZ])*(node3[XX]-node1[XX])) / 2.0;
  m_celldata.nz[0] = ((node3[XX]-node1[XX])*(node2[YY]-node1[YY]) - (node2[XX]-node1[XX])*(node3[YY]-node1[YY])) / 2.0;
  m_celldata.nx[1] = ((node2[YY]-node0[YY])*(node3[ZZ]-node0[ZZ]) - (node3[YY]-node0[YY])*(node2[ZZ]-node0[ZZ])) / 2.0;
  m_celldata.ny[1] = ((node2[ZZ]-node0[ZZ])*(node3[XX]-node0[XX]) - (node3[ZZ]-node0[ZZ])*(node2[XX]-node0[XX])) / 2.0;
  m_celldata.nz[1] = ((node2[XX]-node0[XX])*(node3[YY]-node0[YY]) - (node3[XX]-node0[XX])*(node2[YY]-node0[YY])) / 2.0;
  m_celldata.nx[2] = ((node3[YY]-node0[YY])*(node1[ZZ]-node0[ZZ]) - (node1[YY]-node0[YY])*(node3[ZZ]-node0[ZZ])) / 2.0;
  m_celldata.ny[2] = ((node3[ZZ]-node0[ZZ])*(node1[XX]-node0[XX]) - (node1[ZZ]-node0[ZZ])*(node3[XX]-node0[XX])) / 2.0;
  m_celldata.nz[2] = ((node3[XX]-node0[XX])*(node1[YY]-node0[YY]) - (node1[XX]-node0[XX])*(node3[YY]-node0[YY])) / 2.0;
  m_celldata.nx[3] = ((node1[YY]-node0[YY])*(node2[ZZ]-node0[ZZ]) - (node2[YY]-node0[YY])*(node1[ZZ]-node0[ZZ])) / 2.0;
  m_celldata.ny[3] = ((node1[ZZ]-node0[ZZ])*(node2[XX]-node0[XX]) - (node2[ZZ]-node0[ZZ])*(node1[XX]-node0[XX])) / 2.0;
  m_celldata.nz[3] = ((node1[XX]-node0[XX])*(node2[YY]-node0[YY]) - (node2[XX]-node0[XX])*(node1[YY]-node0[YY])) / 2.0;

}

//////////////////////////////////////////////////////////////////////////////

    }  // end namespace TetraP1P1Cell
  }  // end namespace UFEM
}  // end namespace COOLFluiD
