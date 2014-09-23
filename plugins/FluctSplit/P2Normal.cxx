#include "P2Normal.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////


/**
   * Default constructor without arguments
   */
  P2Normal::P2Normal(): m_weights(3), m_qpPos(3), m_shift(2)
  {

      //Weights at Gauss points
      m_weights[0] = 5./36.;
      m_weights[1] = 8./36.;
      m_weights[2] = 5./36.;

      //Position of quadrature points 
      //Integration domain <-0.25;0.25> considered

      const CFreal s = std::sqrt(0.6);

      m_qpPos[0] = -0.25*s;
      m_qpPos[1] = 0.0;
      m_qpPos[2] = 0.25*s;

      m_shift[0] = 0.25;
      m_shift[1] = 0.75;

  }

  /**
   * Default destructor
   */
  P2Normal::~P2Normal()
  {
  }





///Compute normal perpendicular to a face of P2 triangle. The face has index 0...8
///xi and eta define where the normal is computed in reference space
///The orientation of normals is such that it corresponds to orientation of normals in HOCRD strategy for P1P2 or P2P2 triangles

void P2Normal::ComputeNormal(const std::vector<Framework::Node*>& nodes, const CFuint faceIdx, const CFreal xi, const CFreal eta, RealVector& normal)
{
   cf_assert (nodes.size() == 6);  

   m_workFaceIdx = faceIdx;

   switch (faceIdx) 
   {
    case 5: m_workFaceIdx = 2; //The formula for subfaces 2 and 5 is the same
            break;
    case 6: m_workFaceIdx = 3; //The formula for subfaces 3 and 6 is the same
            break;
    case 7: m_workFaceIdx = 1; //The formula for subfaces 7 and 1 is the same
            break;
    }


   switch (m_workFaceIdx)
   {

   case 0: {
    normal[XX] = -(4*xi-1)*(*nodes[1])[YY]-(4*xi-1)*(*nodes[2])[YY]-2*(*nodes[3])[YY]-\
			(2-8*xi)*(*nodes[4])[YY]+2*(*nodes[5])[YY];
    normal[YY] =  (4*xi-1)*(*nodes[1])[XX]+(4*xi-1)*(*nodes[2])[XX]+2*(*nodes[3])[XX]+\
			(2-8*xi)*(*nodes[4])[XX]-2*(*nodes[5])[XX];
    }
    break;

   case 1: {
    normal[XX] = -(-3+4*eta)*(*nodes[0])[YY]-(4*eta-1)*(*nodes[2])[YY]-(4-8*eta)*(*nodes[5])[YY];
    normal[YY] =  (-3+4*eta)*(*nodes[0])[XX]+(4*eta-1)*(*nodes[2])[XX]+(4-8*eta)*(*nodes[5])[XX];
    }
    break;

   case 2: {
    normal[XX] =  (-3+4*xi)*(*nodes[0])[YY]+(4*xi-1)*(*nodes[1])[YY]+(4-8*xi)*(*nodes[3])[YY];
    normal[YY] = -(-3+4*xi)*(*nodes[0])[XX]-(4*xi-1)*(*nodes[1])[XX]-(4-8*xi)*(*nodes[3])[XX];
    }
    break;

   case 3: {
    normal[XX] = -(4*xi-1)*(*nodes[1])[YY]-(-3+4*xi)*(*nodes[2])[YY]-(4-8*xi)*(*nodes[4])[YY];
    normal[YY] =  (4*xi-1)*(*nodes[1])[XX]+(-3+4*xi)*(*nodes[2])[XX]+(4-8*xi)*(*nodes[4])[XX];
    }
    break;

   case 4: {
    normal[XX] = -(4*eta-1)*(*nodes[0])[YY]-(4*eta-1)*(*nodes[2])[YY]+2*(*nodes[3])[YY]-2*(*nodes[4])[YY]-(2-8*eta)*(*nodes[5])[YY];
    normal[YY] =  (4*eta-1)*(*nodes[0])[XX]+(4*eta-1)*(*nodes[2])[XX]-2*(*nodes[3])[XX]+2*(*nodes[4])[XX]+(2-8*eta)*(*nodes[5])[XX];
    }
    break;

   case 8: {
    normal[XX] =  (4*xi-1)*(*nodes[0])[YY]+(4*xi-1)*(*nodes[1])[YY]+(2-8*xi)*(*nodes[3])[YY]+\
				2*(*nodes[4])[YY]-2*(*nodes[5])[YY];
    normal[YY] = -(4*xi-1)*(*nodes[0])[XX]-(4*xi-1)*(*nodes[1])[XX]-(2-8*xi)*(*nodes[3])[XX]-\
				2*(*nodes[4])[XX]+2*(*nodes[5])[XX];
    }
    break;
   }///switch for m_workFaceIdx 

}




void P2Normal::ComputeBNormal(const std::vector<Framework::Node*>& nodes, const CFreal xi, RealVector& normal)
{

  cf_assert (nodes.size() == 3); //P2 boundary face has 3 nodes 

  normal[XX] =  (4*xi-3)*(*nodes[0])[YY] + (4*xi-1)*(*nodes[1])[YY] + (4-8*xi)*(*nodes[2])[YY];
  normal[YY] = -(4*xi-3)*(*nodes[0])[XX] - (4*xi-1)*(*nodes[1])[XX] - (4-8*xi)*(*nodes[2])[XX];
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
