#include "Environment/ObjectProvider.hh"

#include "UFEM/TriagP1P1Cell/CellProps.hh"
#include "UFEM/TriagP1P1Cell/UFEMTriagP1P1Cell.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TriagP1P1Cell {

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider < CellProps,
                              ElemProps,
                              UFEMTriagP1P1CellPlugin,
                              ElemProps::NARGS >
aTriagP1P1CellProps_Provider ( "TriagP1P1CellProps" );

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

  // volume of the elem
  CFreal d[2][2];
  d[0][0] = node1[XX] - node0[XX];
  d[0][1] = node1[YY] - node0[YY];
  d[1][0] = node2[XX] - node0[XX];
  d[1][1] = node2[YY] - node0[YY];
  m_celldata.vol =(  d[0][0]*d[1][1] - d[1][0]*d[0][1] ) * 0.5;

  // normals of elem faces
  m_celldata.nx[0] = node1[YY] - node2[YY];
  m_celldata.ny[0] = node2[XX] - node1[XX];
  m_celldata.nx[1] = node2[YY] - node0[YY];
  m_celldata.ny[1] = node0[XX] - node2[XX];
  m_celldata.nx[2] = node0[YY] - node1[YY];
  m_celldata.ny[2] = node1[XX] - node0[XX];
}

//////////////////////////////////////////////////////////////////////////////

    }  // end namespace TriagP1P1Cell
  }  // end namespace UFEM
}  // end namespace COOLFluiD
