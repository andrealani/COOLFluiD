#include "P3Normal.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

const CFreal P3Normal::m_XI[10] = { 0.0, 1.0, 0.0, 1.0/3.0, 2.0/3.0, 2.0/3.0, 1.0/3.0, 0.0, 0.0, 1.0/3.0 };
const CFreal P3Normal::m_ETA[10] = { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0/3.0, 2.0/3.0, 2.0/3.0, 1.0/3.0, 1.0/3.0 };
const CFreal P3Normal::m_dLdt_base[3][3] = { {0,1,-1}, {-1,0,1}, {-1,1,0} };
const CFuint P3Normal::m_nVert = 10;
const CFuint P3Normal::m_nSubFaces = 18;


/**
   * Default constructor without arguments
   */
  P3Normal::P3Normal()
  {
    for(CFuint i=0;i<m_nVert;++i) m_dndt[i] = 0.0;
  }

  /**
   * Default destructor
   */
  P3Normal::~P3Normal()
  {
  }





///Compute normal perpendicular to a face of P3 triangle. The face has index 0...17
///xi and eta define where the normal is computed in reference space
///The orientation of normals is such that it corresponds to orientation of normals in HOCRD strategy for P3P3 triangles

void P3Normal::ComputeNormal(const std::vector<Framework::Node*>& nodes, const CFuint faceIdx, const CFreal xi, const CFreal eta, RealVector& normal)
{
   cf_assert (nodes.size() == 10);  

    ///initialize variables
  for(CFuint i=0;i<m_nVert;++i) {
    m_dndt[i] = 0.0;
    m_w[i] = 0.0;
  }

  m_L[0] = 1.0-xi-eta;
  m_L[1] = xi;
  m_L[2] = eta;

  const CFuint baseFaceIdx = faceIdx % 3;
  m_dLdt = m_dLdt_base[baseFaceIdx];

//   m_dndt[0] = 0.5*m_dLdt[0]*(3.0*(3.0*m_L[0]-2.0)*m_L[0] + 3.0*(3.0*m_L[0]-1.0)*m_L[0] + (3.0*m_L[0]-1.0)*(3.0*m_L[0]-2.0));
//   m_dndt[1] = 0.5*m_dLdt[1]*(3.0*(3.0*m_L[1]-2.0)*m_L[1] + 3.0*(3.0*m_L[1]-1.0)*m_L[1] + (3.0*m_L[1]-1.0)*(3.0*m_L[1]-2.0));
//   m_dndt[2] = 0.5*m_dLdt[2]*(3.0*(3.0*m_L[2]-2.0)*m_L[2] + 3.0*(3.0*m_L[2]-1.0)*m_L[2] + (3.0*m_L[2]-1.0)*(3.0*m_L[2]-2.0));

  m_dndt[0] = 0.5*m_dLdt[0]*(9.0*(2.0*m_L[0]-1.0)*m_L[0] + (3.0*m_L[0]-1.0)*(3.0*m_L[0]-2.0));
  m_dndt[1] = 0.5*m_dLdt[1]*(9.0*(2.0*m_L[1]-1.0)*m_L[1] + (3.0*m_L[1]-1.0)*(3.0*m_L[1]-2.0));
  m_dndt[2] = 0.5*m_dLdt[2]*(9.0*(2.0*m_L[2]-1.0)*m_L[2] + (3.0*m_L[2]-1.0)*(3.0*m_L[2]-2.0));

  m_dndt[3] = 4.5*(m_dLdt[0]*m_L[1]*(3.0*m_L[0]-1.0) + m_L[0]*m_dLdt[1]*(3.0*m_L[0]-1.0) + 3.0*m_L[0]*m_L[1]*m_dLdt[0]);
  m_dndt[4] = 4.5*(m_dLdt[0]*m_L[1]*(3.0*m_L[1]-1.0) + m_L[0]*m_dLdt[1]*(3.0*m_L[1]-1.0) + 3.0*m_L[0]*m_L[1]*m_dLdt[1]);
  m_dndt[5] = 4.5*(m_dLdt[1]*m_L[2]*(3.0*m_L[1]-1.0) + m_L[1]*m_dLdt[2]*(3.0*m_L[1]-1.0) + 3.0*m_L[1]*m_L[2]*m_dLdt[1]);
  m_dndt[6] = 4.5*(m_dLdt[1]*m_L[2]*(3.0*m_L[2]-1.0) + m_L[1]*m_dLdt[2]*(3.0*m_L[2]-1.0) + 3.0*m_L[1]*m_L[2]*m_dLdt[2]);
  m_dndt[7] = 4.5*(m_dLdt[0]*m_L[2]*(3.0*m_L[2]-1.0) + m_L[0]*m_dLdt[2]*(3.0*m_L[2]-1.0) + 3.0*m_L[0]*m_L[2]*m_dLdt[2]);
  m_dndt[8] = 4.5*(m_dLdt[0]*m_L[2]*(3.0*m_L[0]-1.0) + m_L[0]*m_dLdt[2]*(3.0*m_L[0]-1.0) + 3.0*m_L[0]*m_L[2]*m_dLdt[0]);

  m_dndt[9] = 27.0*(m_dLdt[0]*m_L[1]*m_L[2] + m_L[0]*m_dLdt[1]*m_L[2] + m_L[0]*m_L[1]*m_dLdt[2]);

  normal[XX] = normal[YY] = 0.0;
  for(CFuint j=0;j<m_nVert;++j) {
    normal[XX] += (*nodes[j])[YY]*m_dndt[j];
    normal[YY] -= (*nodes[j])[XX]*m_dndt[j];
  }

  const CFreal orient = faceIdx %3 == 2 ? 1.0 : -1.0;

  normal *= orient;

  ///don't scale the normal!
  ///this will be done outside if necessary!
//   const CFreal invnorm2 = 1.0/std::sqrt(normal[XX]*normal[XX]+normal[YY]*normal[YY]);
//   normal[XX] *= invnorm2*orient;
//   normal[YY] *= invnorm2*orient;

}

void P3Normal::ComputeBNormal(const std::vector<Framework::Node*>& nodes, const CFreal xi, RealVector& normal)
{

  cf_assert (nodes.size() == 4); //P3 boundary face has 4 nodes 

  m_L[0] = 1.0 - xi;
  m_L[1] = xi;

  m_dN0dxi = -0.5*( 9.0*m_L[0]*(2.0*m_L[0]-1.0) + (3.0*m_L[0]-1.0)*(3.0*m_L[0]-2.0) );
  m_dN1dxi =  0.5*( 9.0*m_L[1]*(2.0*m_L[1]-1.0) + (3.0*m_L[1]-1.0)*(3.0*m_L[1]-2.0) );
  m_dN2dxi = 4.5*( (m_L[0]-m_L[1])*(3.0*m_L[0]-1.0) - 3.0*m_L[0]*m_L[1] );
  m_dN3dxi = 4.5*( (m_L[0]-m_L[1])*(3.0*m_L[1]-1.0) + 3.0*m_L[0]*m_L[1] );

  normal[XX] =  m_dN0dxi * (*nodes[0])[YY] + m_dN1dxi * (*nodes[1])[YY] + m_dN2dxi * (*nodes[2])[YY] + m_dN3dxi * (*nodes[3])[YY];
  normal[YY] = -m_dN0dxi * (*nodes[0])[XX] - m_dN1dxi * (*nodes[1])[XX] - m_dN2dxi * (*nodes[2])[XX] - m_dN3dxi * (*nodes[3])[XX];

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
