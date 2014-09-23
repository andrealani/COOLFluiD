


#include "SpectralFV/SimplexGaussIntegrator.hh"

#include "Common/NotImplementedException.hh"
#include "Common/ShouldNotBeHereException.hh"
#include "Framework/BadFormatException.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

SimplexGaussIntegrator::SimplexGaussIntegrator() :
  m_quadPntsMappedCoor(),
  m_quadPntsWheights(),
  m_quadPntsNormals()
{
  CFAUTOTRACE;

  m_dimensionality  = DIM_0D;
  m_integratorOrder = CFPolyOrder::ORDER0;

  resetQuadPntsMappedCoorAndWheights(m_dimensionality, m_integratorOrder);
}

//////////////////////////////////////////////////////////////////////////////

SimplexGaussIntegrator::SimplexGaussIntegrator(const CFDim dimensionality, const CFPolyOrder::Type integratorOrder) :
  m_quadPntsMappedCoor(),
  m_quadPntsWheights()
{
  CFAUTOTRACE;

  m_dimensionality  = dimensionality;
  m_integratorOrder = integratorOrder;

  resetQuadPntsMappedCoorAndWheights(dimensionality, integratorOrder);
}

//////////////////////////////////////////////////////////////////////////////

SimplexGaussIntegrator::~SimplexGaussIntegrator()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void SimplexGaussIntegrator::resetQuadPntsMappedCoorAndWheights(const CFDim dimensionality,
                                                                const CFPolyOrder::Type integratorOrder)
{
  CFAUTOTRACE;

  switch (dimensionality)
  {
    case DIM_0D:
    {

      // Set the number of quadrature points
      m_nbrQuadPnts = 1;

      // resize the variables
      m_quadPntsWheights.resize(m_nbrQuadPnts);
      m_quadPntsMappedCoor.resize(m_nbrQuadPnts);
      m_quadPntsMappedCoor[0].resize(dimensionality);

      // set the coordinates

      // set the wheights
      m_quadPntsWheights[0] = 1.;

    } break;
    case DIM_1D:
    {
      switch (integratorOrder)
      {

        case CFPolyOrder::ORDER0:
        case CFPolyOrder::ORDER1:
        {
          // Set the number of quadrature points
          m_nbrQuadPnts = 1;

          // resize the variables
          m_quadPntsWheights.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor[0].resize(dimensionality);

          // set the coordinates
          m_quadPntsMappedCoor[0][KSI] = 0.5;

          // set the wheights
          m_quadPntsWheights[0] = 1.;

        } break;
        case CFPolyOrder::ORDER2:
        case CFPolyOrder::ORDER3:
        {
          // Set the number of quadrature points
          m_nbrQuadPnts = 2;

          // resize the variables
          m_quadPntsWheights.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor[0].resize(dimensionality);
          m_quadPntsMappedCoor[1].resize(dimensionality);

          // set the coordinates
          m_quadPntsMappedCoor[0][KSI] = 0.5*(1. - 1./sqrt(3.));

          m_quadPntsMappedCoor[1][KSI] = 0.5*(1. + 1./sqrt(3.));

          // set the wheights
          m_quadPntsWheights[0] = 0.5;

          m_quadPntsWheights[1] = 0.5;

        } break;
        case CFPolyOrder::ORDER4:
        case CFPolyOrder::ORDER5:
        {
          // Set the number of quadrature points
          m_nbrQuadPnts = 3;

          // resize the variables
          m_quadPntsWheights.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor.resize(m_nbrQuadPnts);
          for (CFuint iPnt = 0; iPnt < m_nbrQuadPnts; ++iPnt)
            m_quadPntsMappedCoor[iPnt].resize(dimensionality);

          // set the coordinates
          m_quadPntsMappedCoor[0][KSI] = 0.5;

          m_quadPntsMappedCoor[1][KSI] = 0.5*(1. - sqrt(3./5.));

          m_quadPntsMappedCoor[2][KSI] = 0.5*(1. + sqrt(3./5.));

          // set the wheights
          m_quadPntsWheights[0] = 8.0/18.0;

          m_quadPntsWheights[1] = 5.0/18.0;

          m_quadPntsWheights[2] = 5.0/18.0;

        } break;
        default:
        {
          throw Common::NotImplementedException (FromHere(),"SimplexGaussIntegrator only implemented for polynomials up to 5th order.");
        }
      }

    } break;
    case DIM_2D:
    {

      switch (integratorOrder)
      {
        case CFPolyOrder::ORDER0:
        case CFPolyOrder::ORDER1:
        {
          // Set the number of quadrature points
          m_nbrQuadPnts = 1;

          // resize the variables
          m_quadPntsWheights.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor[0].resize(dimensionality);

          // set the coordinates
          m_quadPntsMappedCoor[0][KSI] = 1./3.;
          m_quadPntsMappedCoor[0][ETA] = 1./3.;

          // set the wheights
          m_quadPntsWheights[0] = 0.5;

        } break;
        case CFPolyOrder::ORDER2:
        {
          // Set the number of quadrature points
          m_nbrQuadPnts = 3;

          // resize the variables
          m_quadPntsWheights.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor.resize(m_nbrQuadPnts);
          for (CFuint iPnt = 0; iPnt < m_nbrQuadPnts; ++iPnt)
            m_quadPntsMappedCoor[iPnt].resize(dimensionality);

          // set the coordinates
          m_quadPntsMappedCoor[0][KSI] = 0.5;
          m_quadPntsMappedCoor[0][ETA] = 0.5;

          m_quadPntsMappedCoor[1][KSI] = 0.0;
          m_quadPntsMappedCoor[1][ETA] = 0.5;

          m_quadPntsMappedCoor[2][KSI] = 0.5;
          m_quadPntsMappedCoor[2][ETA] = 0.0;

          // set the wheights
          m_quadPntsWheights[0] = 1./6.;
          m_quadPntsWheights[1] = 1./6.;
          m_quadPntsWheights[2] = 1./6.;

        } break;
        case CFPolyOrder::ORDER3:
        {
          // Set the number of quadrature points
          m_nbrQuadPnts = 4;

          // resize the variables
          m_quadPntsWheights.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor.resize(m_nbrQuadPnts);
          for (CFuint iPnt = 0; iPnt < m_nbrQuadPnts; ++iPnt)
            m_quadPntsMappedCoor[iPnt].resize(dimensionality);

          // set the coordinates
          m_quadPntsMappedCoor[0][KSI] = 1./3.;
          m_quadPntsMappedCoor[0][ETA] = 1./3.;

          m_quadPntsMappedCoor[1][KSI] = 0.2;
          m_quadPntsMappedCoor[1][ETA] = 0.2;

          m_quadPntsMappedCoor[2][KSI] = 0.2;
          m_quadPntsMappedCoor[2][ETA] = 0.6;

          m_quadPntsMappedCoor[3][KSI] = 0.6;
          m_quadPntsMappedCoor[3][ETA] = 0.2;

          // set the wheights
          m_quadPntsWheights[0] = -27./96.;
          m_quadPntsWheights[1] = 25./96.;
          m_quadPntsWheights[2] = 25./96.;
          m_quadPntsWheights[3] = 25./96.;

        } break;
        case CFPolyOrder::ORDER4:
        case CFPolyOrder::ORDER5:
        {
          // Set the number of quadrature points
          m_nbrQuadPnts = 7;

          // resize the variables
          m_quadPntsWheights.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor.resize(m_nbrQuadPnts);
          for (CFuint iPnt = 0; iPnt < m_nbrQuadPnts; ++iPnt)
            m_quadPntsMappedCoor[iPnt].resize(dimensionality);

          // set the coordinates
          m_quadPntsMappedCoor[0][KSI] = 1./3.;
          m_quadPntsMappedCoor[0][ETA] = 1./3.;

          m_quadPntsMappedCoor[1][KSI] = 0.0597158717;
          m_quadPntsMappedCoor[1][ETA] = 0.4701420641;

          m_quadPntsMappedCoor[2][KSI] = 0.4701420641;
          m_quadPntsMappedCoor[2][ETA] = 0.0597158717;

          m_quadPntsMappedCoor[3][KSI] = 0.4701420641;
          m_quadPntsMappedCoor[3][ETA] = 0.4701420641;

          m_quadPntsMappedCoor[4][KSI] = 0.7974269853;
          m_quadPntsMappedCoor[4][ETA] = 0.1012865073;

          m_quadPntsMappedCoor[5][KSI] = 0.1012865073;
          m_quadPntsMappedCoor[5][ETA] = 0.7974269853;

          m_quadPntsMappedCoor[6][KSI] = 0.1012865073;
          m_quadPntsMappedCoor[6][ETA] = 0.1012865073;

          // set the wheights
          m_quadPntsWheights[0] = 0.5 * 0.225;

          m_quadPntsWheights[1] = 0.5 * 0.1323941527;

          m_quadPntsWheights[2] = 0.5 * 0.1323941527;

          m_quadPntsWheights[3] = 0.5 * 0.1323941527;

          m_quadPntsWheights[4] = 0.5 * 0.1259391805;

          m_quadPntsWheights[5] = 0.5 * 0.1259391805;

          m_quadPntsWheights[6] = 0.5 * 0.1259391805;

        } break;
        default:
        {
          throw Common::NotImplementedException (FromHere(),"SimplexGaussIntegrator only implemented for polynomials up to 5th order.");
        }
      }

    } break;
    case DIM_3D:
    {

      switch (integratorOrder)
      {
        case CFPolyOrder::ORDER0:
        case CFPolyOrder::ORDER1:
        {
          // Set the number of quadrature points
          m_nbrQuadPnts = 1;

          // resize the variables
          m_quadPntsWheights.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor[0].resize(dimensionality);

          // set the coordinates
          m_quadPntsMappedCoor[0][KSI] = 0.25;
          m_quadPntsMappedCoor[0][ETA] = 0.25;
          m_quadPntsMappedCoor[0][ZTA] = 0.25;

          // set the wheights
          m_quadPntsWheights[0] = 1.0/6.0;

        } break;
        case CFPolyOrder::ORDER2:
        {
          // Set the number of quadrature points
          m_nbrQuadPnts = 4;

          // resize the variables
          m_quadPntsWheights.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor.resize(m_nbrQuadPnts);
          for (CFuint iPnt = 0; iPnt < m_nbrQuadPnts; ++iPnt)
            m_quadPntsMappedCoor[iPnt].resize(dimensionality);

          // set the coordinates
          const CFreal a = (5. - sqrt(5.0))/20.;
          const CFreal b = (5. + 3.*sqrt(5.0))/20.;

          m_quadPntsMappedCoor[0][KSI] = a;
          m_quadPntsMappedCoor[0][ETA] = a;
          m_quadPntsMappedCoor[0][ZTA] = a;

          m_quadPntsMappedCoor[1][KSI] = a;
          m_quadPntsMappedCoor[1][ETA] = a;
          m_quadPntsMappedCoor[1][ZTA] = b;

          m_quadPntsMappedCoor[2][KSI] = a;
          m_quadPntsMappedCoor[2][ETA] = b;
          m_quadPntsMappedCoor[2][ZTA] = a;

          m_quadPntsMappedCoor[3][KSI] = b;
          m_quadPntsMappedCoor[3][ETA] = a;
          m_quadPntsMappedCoor[3][ZTA] = a;

          // set the wheights
          m_quadPntsWheights[0] = 1.0/24.0;

          m_quadPntsWheights[1] = 1.0/24.0;

          m_quadPntsWheights[2] = 1.0/24.0;

          m_quadPntsWheights[3] = 1.0/24.0;

        } break;
        case CFPolyOrder::ORDER3:
        {
          // Set the number of quadrature points
          m_nbrQuadPnts = 5;

          // resize the variables
          m_quadPntsWheights.resize(m_nbrQuadPnts);
          m_quadPntsMappedCoor.resize(m_nbrQuadPnts);
          for (CFuint iPnt = 0; iPnt < m_nbrQuadPnts; ++iPnt)
            m_quadPntsMappedCoor[iPnt].resize(dimensionality);

          // set the coordinates
          const CFreal a = 0.25 ;
          const CFreal b = 1./6.;
          const CFreal c = 0.5;

          m_quadPntsMappedCoor[0][KSI] = a;
          m_quadPntsMappedCoor[0][ETA] = a;
          m_quadPntsMappedCoor[0][ZTA] = a;

          m_quadPntsMappedCoor[1][KSI] = b;
          m_quadPntsMappedCoor[1][ETA] = b;
          m_quadPntsMappedCoor[1][ZTA] = b;

          m_quadPntsMappedCoor[2][KSI] = b;
          m_quadPntsMappedCoor[2][ETA] = b;
          m_quadPntsMappedCoor[2][ZTA] = c;

          m_quadPntsMappedCoor[3][KSI] = b;
          m_quadPntsMappedCoor[3][ETA] = c;
          m_quadPntsMappedCoor[3][ZTA] = b;

          m_quadPntsMappedCoor[4][KSI] = c;
          m_quadPntsMappedCoor[4][ETA] = b;
          m_quadPntsMappedCoor[4][ZTA] = b;

          // set the wheights
          m_quadPntsWheights[0] = -2.0/15.0;

          m_quadPntsWheights[1] = 3.0/40.0;

          m_quadPntsWheights[2] = 3.0/40.0;

          m_quadPntsWheights[3] = 3.0/40.0;

          m_quadPntsWheights[4] = 3.0/40.0;

        } break;
        default:
        {
          throw Common::NotImplementedException (FromHere(),"SimplexGaussIntegrator not implemented for 3nd order tetra's.");
        }
      }

    } break;
  }
}

//////////////////////////////////////////////////////////////////////////////

void SimplexGaussIntegrator::setDimensionality(const CFDim dimensionality)
{
  m_dimensionality = dimensionality;
  resetQuadPntsMappedCoorAndWheights(m_dimensionality,m_integratorOrder);
  calcNbrPolyTerms4ExIntegr();
}

//////////////////////////////////////////////////////////////////////////////

void SimplexGaussIntegrator::setIntegratorOrder(const CFPolyOrder::Type integratorOrder)
{
  m_integratorOrder = integratorOrder;
  resetQuadPntsMappedCoorAndWheights(m_dimensionality,m_integratorOrder);
  calcNbrPolyTerms4ExIntegr();
}

//////////////////////////////////////////////////////////////////////////////

void SimplexGaussIntegrator::calcNbrPolyTerms4ExIntegr()
{
  m_nbrPolyTerms4ExIntegr = m_integratorOrder + 1;
  for (CFuint iDim = 2; iDim <= m_dimensionality; ++iDim)
    m_nbrPolyTerms4ExIntegr *= (m_integratorOrder + iDim)/iDim;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector > SimplexGaussIntegrator::getQuadPntsCoords(const vector< RealVector >& simplexNodeCoord)
{
  // some checks for correct format
  cf_assert(simplexNodeCoord.size() == m_dimensionality+1);
  cf_assert(m_dimensionality != DIM_0D);
  cf_assert(simplexNodeCoord[0].size() == m_dimensionality);

  // return variable
  vector< RealVector > qNodeCoord(m_nbrQuadPnts);

  // compute the quadrature node global coordinates
  for (CFuint iQNode = 0; iQNode < m_nbrQuadPnts; ++iQNode)
  {
    // resize the variable
    qNodeCoord[iQNode].resize(m_dimensionality);

    // compute global coordinates
    qNodeCoord[iQNode] = simplexNodeCoord[0]*(1.0 - m_quadPntsMappedCoor[iQNode].sum());
    for (CFuint iSNode = 0; iSNode < m_dimensionality; ++iSNode)
    {
      qNodeCoord[iQNode] += simplexNodeCoord[iSNode+1]*m_quadPntsMappedCoor[iQNode][iSNode];
    }
  }

  return qNodeCoord;
}

//////////////////////////////////////////////////////////////////////////////

vector< CFreal > SimplexGaussIntegrator::getQuadPntsWheights(const vector< RealVector >& simplexNodeCoord)
{
  // some checks for correct format
  cf_assert(simplexNodeCoord.size() == m_dimensionality+1);
  cf_assert(simplexNodeCoord[0].size() == m_dimensionality);

  // return variable
  vector< CFreal > qNodeWheights(m_nbrQuadPnts);

  // variable for jacobian determinant
  CFreal jacobDet;

  // compute jacobian determinant
  switch (m_dimensionality)
  {
    case DIM_0D:
    {
      jacobDet = 1.;
    }
    case DIM_1D:
    {
      jacobDet = simplexNodeCoord[1][XX] - simplexNodeCoord[0][XX];
    } break;
    case DIM_2D:
    {
      RealMatrix jacobian(2,2);

      for (CFuint iVec = 0; iVec < 2; ++iVec)
      {
        jacobian.setColumn(RealVector(simplexNodeCoord[iVec+1] - simplexNodeCoord[0]),iVec);
      }

      jacobDet = abs(jacobian.determ2());
    } break;
    case DIM_3D:
    {
      RealMatrix jacobian(3,3);

      for (CFuint iVec = 0; iVec < 3; ++iVec)
      {
        jacobian.setColumn(RealVector(simplexNodeCoord[iVec+1] - simplexNodeCoord[0]),iVec);
      }

      jacobDet = abs(jacobian.determ3());
    } break;
    default:
      throw Common::ShouldNotBeHereException (FromHere(),"Invalid dimensionality");
  }

  // compute the quadrature node wheights
  for (CFuint iQNode = 0; iQNode < m_nbrQuadPnts; ++iQNode)
  {
    qNodeWheights[iQNode] = jacobDet*m_quadPntsWheights[iQNode];
  }

  return qNodeWheights;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector > SimplexGaussIntegrator::getQuadPntsCoordsPlus1D
                                                              (const vector< RealVector >& simplexNodeCoord)
{
  // helper variable
  const CFuint dimP1 = m_dimensionality + 1;

  // some checks for correct format
  cf_assert(m_dimensionality != DIM_3D);
  cf_assert(simplexNodeCoord.size() == dimP1);
  cf_assert(simplexNodeCoord[0].size() == dimP1);

  // return variable
  vector< RealVector > qNodeCoord(m_nbrQuadPnts);

  // compute the quadrature node global coordinates
  if (m_dimensionality != DIM_0D) // not entirely sure if this if() is necessary...
  {
    for (CFuint iQNode = 0; iQNode < m_nbrQuadPnts; ++iQNode)
    {
      // resize the variable
      qNodeCoord[iQNode].resize(dimP1);

      // compute global coordinates
      qNodeCoord[iQNode] = simplexNodeCoord[0]*(1.0 - m_quadPntsMappedCoor[iQNode].sum());
      for (CFuint iSNode = 0; iSNode < m_dimensionality; ++iSNode)
      {
        qNodeCoord[iQNode] += simplexNodeCoord[iSNode+1]*m_quadPntsMappedCoor[iQNode][iSNode];
      }
    }
  }
  else
  {
    qNodeCoord[0] = simplexNodeCoord[0];
  }

  return qNodeCoord;
}

//////////////////////////////////////////////////////////////////////////////

vector< CFreal > SimplexGaussIntegrator::getQuadPntsWheightsPlus1D(const vector< RealVector >& simplexNodeCoord)
{
  // helper variable
  const CFuint dimP1 = m_dimensionality + 1;

  // some checks for correct format
  cf_assert(m_dimensionality != DIM_3D);
  cf_assert(simplexNodeCoord.size() == dimP1);
  cf_assert(simplexNodeCoord[0].size() == dimP1);

  // return variable
  vector< CFreal > qNodeWheights(m_nbrQuadPnts);

  // variable for jacobian determinant
  CFreal jacobDet;

  // compute jacobian determinant
  switch (m_dimensionality)
  {
    case DIM_0D:
    {
      jacobDet = 1.;
    }
    case DIM_1D:
    {
      jacobDet = sqrt(pow(simplexNodeCoord[1][XX] - simplexNodeCoord[0][XX],2.0) +
                      pow(simplexNodeCoord[1][YY] - simplexNodeCoord[0][YY],2.0));
    } break;
    case DIM_2D:
    {
      RealVector vec1 = simplexNodeCoord[1] - simplexNodeCoord[0];
      RealVector vec2 = simplexNodeCoord[2] - simplexNodeCoord[0];

      RealVector normVect(3);
      normVect[XX] = vec1[YY]*vec2[ZZ] - vec1[ZZ]*vec2[YY];
      normVect[YY] = vec1[ZZ]*vec2[XX] - vec1[XX]*vec2[ZZ];
      normVect[ZZ] = vec1[XX]*vec2[YY] - vec1[YY]*vec2[XX];

      jacobDet = sqrt(pow(normVect[XX],2.0) + pow(normVect[YY],2.0) + pow(normVect[ZZ],2.0));
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Invalid dimensionality");
    }
  }

  // compute the quadrature node wheights
  for (CFuint iQNode = 0; iQNode < m_nbrQuadPnts; ++iQNode)
  {
    qNodeWheights[iQNode] = jacobDet*m_quadPntsWheights[iQNode];
  }

  return qNodeWheights;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV
}  // namespace COOLFluiD
