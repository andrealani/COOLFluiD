#include "Framework/BadFormatException.hh"
#include "Common/ShouldNotBeHereException.hh"

#include "MathTools/RealMatrix.hh"

#include "Common/NotImplementedException.hh"

#include "SpectralFD/TensorProductGaussIntegrator.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

namespace COOLFluiD {
  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

TensorProductGaussIntegrator::TensorProductGaussIntegrator() :
  m_quadPntsMappedCoor(),
  m_quadPntsWheights(),
  m_quadPntsMappedCoor1D(),
  m_quadPntsWheights1D()
{
  CFAUTOTRACE;

  m_dimensionality  = DIM_0D;
  m_integratorOrder = CFPolyOrder::ORDER0;

  resetQuadPntsMappedCoorAndWheights1D(m_integratorOrder);
  resetQuadPntsMappedCoorAndWheights  (m_dimensionality );
}

//////////////////////////////////////////////////////////////////////////////

TensorProductGaussIntegrator::TensorProductGaussIntegrator(const CFDim dimensionality, const CFPolyOrder::Type integratorOrder) :
  m_quadPntsMappedCoor(),
  m_quadPntsWheights()
{
  CFAUTOTRACE;

  m_dimensionality  = dimensionality;

  // taking into account the contribution of the Jacobian
  m_integratorOrder = static_cast<CFPolyOrder::Type>( integratorOrder + dimensionality - 1);

  resetQuadPntsMappedCoorAndWheights1D(m_integratorOrder);
  resetQuadPntsMappedCoorAndWheights  (m_dimensionality );
}

//////////////////////////////////////////////////////////////////////////////


TensorProductGaussIntegrator::~TensorProductGaussIntegrator()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void TensorProductGaussIntegrator::resetQuadPntsMappedCoorAndWheights(const CFDim dimensionality)
{
  CFAUTOTRACE;

  // set the number of quadrature points
  m_nbrQuadPnts = 1;
  for (CFuint iDim = 0; iDim < static_cast<CFuint>(m_dimensionality); ++iDim)
  {
    m_nbrQuadPnts *= m_nbrQuadPnts1D;
  }

  // resize the variables
  m_quadPntsWheights  .resize(m_nbrQuadPnts);
  m_quadPntsMappedCoor.resize(m_nbrQuadPnts);
  for (CFuint iQPnt = 0; iQPnt < m_nbrQuadPnts; ++iQPnt)
  {
    m_quadPntsMappedCoor[iQPnt].resize(dimensionality);
  }

  switch (dimensionality)
  {
    case DIM_0D:
    {
      m_quadPntsWheights[0] = 1;
    } break;
    case DIM_1D:
    {
      // set the coordinates and the wheights
      CFuint iQPnt = 0;
      for (CFuint iKsi = 0; iKsi < m_nbrQuadPnts1D; ++iKsi, ++iQPnt)
      {
        m_quadPntsMappedCoor[iQPnt][KSI] = m_quadPntsMappedCoor1D[iKsi];

        m_quadPntsWheights  [iQPnt] = m_quadPntsWheights1D[iKsi];
      }
    } break;
    case DIM_2D:
    {
      // set the coordinates and the wheights
      CFuint iQPnt = 0;
      for (CFuint iKsi = 0; iKsi < m_nbrQuadPnts1D; ++iKsi)
      {
        for (CFuint iEta = 0; iEta < m_nbrQuadPnts1D; ++iEta, ++iQPnt)
        {
          m_quadPntsMappedCoor[iQPnt][KSI] = m_quadPntsMappedCoor1D[iKsi];
          m_quadPntsMappedCoor[iQPnt][ETA] = m_quadPntsMappedCoor1D[iEta];

          m_quadPntsWheights  [iQPnt] = m_quadPntsWheights1D[iKsi]*
                                        m_quadPntsWheights1D[iEta];
        }
      }
    } break;
    case DIM_3D:
    {
      // set the coordinates and the wheights
      CFuint iQPnt = 0;
      for (CFuint iKsi = 0; iKsi < m_nbrQuadPnts1D; ++iKsi)
      {
        for (CFuint iEta = 0; iEta < m_nbrQuadPnts1D; ++iEta)
        {
          for (CFuint iZta = 0; iZta < m_nbrQuadPnts1D; ++iZta, ++iQPnt)
          {
            m_quadPntsMappedCoor[iQPnt][KSI] = m_quadPntsMappedCoor1D[iKsi];
            m_quadPntsMappedCoor[iQPnt][ETA] = m_quadPntsMappedCoor1D[iEta];
            m_quadPntsMappedCoor[iQPnt][ZTA] = m_quadPntsMappedCoor1D[iZta];

            m_quadPntsWheights  [iQPnt] = m_quadPntsWheights1D[iKsi]*
                                          m_quadPntsWheights1D[iEta]*
                                          m_quadPntsWheights1D[iZta];
          }
        }
      }
    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"Unsupported dimensionality");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void TensorProductGaussIntegrator::resetQuadPntsMappedCoorAndWheights1D(const CFPolyOrder::Type integratorOrder)
{
  CFAUTOTRACE;

  switch (integratorOrder)
  {
    case CFPolyOrder::ORDER0:
    case CFPolyOrder::ORDER1:
    {
      // Set the number of quadrature points
      m_nbrQuadPnts1D = 1;

      // resize the variables
      m_quadPntsWheights1D  .resize(m_nbrQuadPnts1D);
      m_quadPntsMappedCoor1D.resize(m_nbrQuadPnts1D);

      // set the coordinates
      m_quadPntsMappedCoor1D[0] = 0.0;

      // set the wheights
      m_quadPntsWheights1D  [0] = 2.;

    } break;
    case CFPolyOrder::ORDER2:
    case CFPolyOrder::ORDER3:
    {
      // Set the number of quadrature points
      m_nbrQuadPnts1D = 2;

      // resize the variables
      m_quadPntsWheights1D  .resize(m_nbrQuadPnts1D);
      m_quadPntsMappedCoor1D.resize(m_nbrQuadPnts1D);

      // set the coordinates
      m_quadPntsMappedCoor1D[0] = -1./sqrt(3.);
      m_quadPntsMappedCoor1D[1] = +1./sqrt(3.);

      // set the wheights
      m_quadPntsWheights1D  [0] = 1.0;
      m_quadPntsWheights1D  [1] = 1.0;

    } break;
    case CFPolyOrder::ORDER4:
    case CFPolyOrder::ORDER5:
    {
      // Set the number of quadrature points
      m_nbrQuadPnts1D = 3;

      // resize the variables
      m_quadPntsWheights1D  .resize(m_nbrQuadPnts1D);
      m_quadPntsMappedCoor1D.resize(m_nbrQuadPnts1D);

      // set the coordinates
      m_quadPntsMappedCoor1D[0] = 0.0;

      m_quadPntsMappedCoor1D[1] = -sqrt(3./5.);

      m_quadPntsMappedCoor1D[2] = +sqrt(3./5.);

      // set the wheights
      m_quadPntsWheights1D  [0] = 8.0/9.0;

      m_quadPntsWheights1D  [1] = 5.0/9.0;

      m_quadPntsWheights1D  [2] = 5.0/9.0;

    } break;
    case CFPolyOrder::ORDER6:
    case CFPolyOrder::ORDER7:
    {
      // Set the number of quadrature points
      m_nbrQuadPnts1D = 4;

      // resize the variables
      m_quadPntsWheights1D  .resize(m_nbrQuadPnts1D);
      m_quadPntsMappedCoor1D.resize(m_nbrQuadPnts1D);

      // set the coordinates
      m_quadPntsMappedCoor1D[0] = -sqrt((3.+2.*sqrt(6./5.))/7.);

      m_quadPntsMappedCoor1D[1] = -sqrt((3.-2.*sqrt(6./5.))/7.);

      m_quadPntsMappedCoor1D[2] = +sqrt((3.-2.*sqrt(6./5.))/7.);

      m_quadPntsMappedCoor1D[3] = +sqrt((3.+2.*sqrt(6./5.))/7.);

      // set the wheights
      m_quadPntsWheights1D  [0] = (18.-sqrt(30.))/36.;

      m_quadPntsWheights1D  [1] = (18.+sqrt(30.))/36.;

      m_quadPntsWheights1D  [2] = (18.+sqrt(30.))/36.;

      m_quadPntsWheights1D  [3] = (18.-sqrt(30.))/36.;

    } break;
    default:
    {
      throw Common::NotImplementedException (FromHere(),"TensorProductGaussIntegrator only implemented for polynomials up to 5th order.");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void TensorProductGaussIntegrator::setDimensionality(const CFDim dimensionality)
{
  m_integratorOrder = static_cast<CFPolyOrder::Type>( m_integratorOrder - m_dimensionality + dimensionality );// taking into account the contribution of the Jacobian
  m_dimensionality = dimensionality;
  resetQuadPntsMappedCoorAndWheights1D(m_integratorOrder);
  resetQuadPntsMappedCoorAndWheights  (m_dimensionality );
  calcNbrPolyTerms4ExIntegr();
}

//////////////////////////////////////////////////////////////////////////////

void TensorProductGaussIntegrator::setIntegratorOrder(const CFPolyOrder::Type integratorOrder)
{
  m_integratorOrder = static_cast<CFPolyOrder::Type>(integratorOrder + m_dimensionality - 1);// taking into account the contribution of the Jacobian
  resetQuadPntsMappedCoorAndWheights1D(m_integratorOrder);
  resetQuadPntsMappedCoorAndWheights  (m_dimensionality );
  calcNbrPolyTerms4ExIntegr();
}

//////////////////////////////////////////////////////////////////////////////

void TensorProductGaussIntegrator::calcNbrPolyTerms4ExIntegr()
{
  m_nbrPolyTerms4ExIntegr = 1;
  for (CFuint iDim = 0; iDim < static_cast<CFuint>(m_dimensionality); ++iDim)
  {
    m_nbrPolyTerms4ExIntegr *= (m_integratorOrder + 1);
  }
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector > TensorProductGaussIntegrator::getQuadPntsCoords(const vector< RealVector >& nodeCoord)
{
  // number of nodes
  const CFuint nbrNodes = nodeCoord.size();

  // some checks for correct format
  cf_assert(nbrNodes == static_cast<CFuint>(pow(static_cast< CFreal >(2),static_cast< CFreal >(m_dimensionality))));
  cf_assert(nodeCoord[0].size() == static_cast<CFuint>(m_dimensionality));

  // return variable
  vector< RealVector > qNodeCoord(m_nbrQuadPnts,RealVector(m_dimensionality));

  // compute the quadrature node global coordinates
  for (CFuint iQPnt = 0; iQPnt < m_nbrQuadPnts; ++iQPnt)
  {
    // evaluate geometric basis functions
    vector< CFreal > geoBasisFuncs = evaluateBasisFunction(m_dimensionality,
                                                           CFPolyOrder::ORDER1,
                                                           m_quadPntsMappedCoor[iQPnt]);
    cf_assert(nbrNodes == geoBasisFuncs.size());

    // compute global coordinates
    for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
    {
      cf_assert(nodeCoord[iNode].size() == static_cast<CFuint>(m_dimensionality));
      qNodeCoord[iQPnt] += nodeCoord[iNode]*geoBasisFuncs[iNode];
//       CF_DEBUG_OBJ(qNodeCoord[iQPnt]);
    }
  }

  return qNodeCoord;
}

//////////////////////////////////////////////////////////////////////////////

vector< CFreal > TensorProductGaussIntegrator::getQuadPntsWheights(const vector< RealVector >& nodeCoord)
{
  // number of nodes
  const CFuint nbrNodes = nodeCoord.size();

  // some checks for correct format
  cf_assert(nbrNodes == static_cast<CFuint>(pow(static_cast< CFreal >(2),static_cast< CFreal >(m_dimensionality))));
  cf_assert(nodeCoord[0].size() == static_cast<CFuint>(m_dimensionality));

  // return variable
  vector< CFreal > qNodeWheights(m_nbrQuadPnts);

  // compute jacobian determinant
  switch (m_dimensionality)
  {
    case DIM_0D:
    {
      qNodeWheights[0] = 1.0;
    } break;
    case DIM_1D:
    {
      // compute Jacobian determinant
      const CFreal jacobDet = 0.5*(nodeCoord[1][XX] - nodeCoord[0][XX]);

      // compute the quadrature node wheights
      for (CFuint iQPnt = 0; iQPnt < m_nbrQuadPnts; ++iQPnt)
      {
        qNodeWheights[iQPnt] = jacobDet*m_quadPntsWheights[iQPnt];
      }
    } break;
    case DIM_2D:
    {
      // compute the quadrature node wheights
      for (CFuint iQPnt = 0; iQPnt < m_nbrQuadPnts; ++iQPnt)
      {
        // evaluate derivatives of geometrical basis functions
        vector< RealVector > geoBasisFuncDerivs = evaluateDerivBasisFunction(m_dimensionality,
                                                                             CFPolyOrder::ORDER1,
                                                                             m_quadPntsMappedCoor[iQPnt]);
        cf_assert(nbrNodes == geoBasisFuncDerivs.size());

        // construct Jacobian matrix
        RealMatrix jacobian(m_dimensionality,m_dimensionality);
        for (CFuint iCoor = 0; iCoor < static_cast<CFuint>(m_dimensionality); ++iCoor)
        {
          RealVector column(m_dimensionality);
          for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
          {
            column += geoBasisFuncDerivs[iNode][iCoor]*nodeCoord[iNode];
          }
          jacobian.setColumn(column,iCoor);
        }

        // compute determinant of jacobian matrix
        const CFreal jacobDet = jacobian.determ2();

        // compute quadrature wheight
        qNodeWheights[iQPnt] = jacobDet*m_quadPntsWheights[iQPnt];
      }
    } break;
    case DIM_3D:
    {
      // compute the quadrature node wheights
      for (CFuint iQPnt = 0; iQPnt < m_nbrQuadPnts; ++iQPnt)
      {
        // evaluate derivatives of geometrical basis functions
        vector< RealVector > geoBasisFuncDerivs = evaluateDerivBasisFunction(m_dimensionality,
                                                                             CFPolyOrder::ORDER1,
                                                                             m_quadPntsMappedCoor[iQPnt]);
        cf_assert(nbrNodes == geoBasisFuncDerivs.size());

        // construct Jacobian matrix
        RealMatrix jacobian(m_dimensionality,m_dimensionality);
        for (CFuint iCoor = 0; iCoor < static_cast<CFuint>(m_dimensionality); ++iCoor)
        {
          RealVector column(m_dimensionality);
          for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
          {
            column += geoBasisFuncDerivs[iNode][iCoor]*nodeCoord[iNode];
          }
          jacobian.setColumn(column,iCoor);
        }

        // compute determinant of jacobian matrix
        const CFreal jacobDet = jacobian.determ3();

        // compute quadrature wheight
        qNodeWheights[iQPnt] = jacobDet*m_quadPntsWheights[iQPnt];
      }
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Invalid dimensionality");
    }
  }

  return qNodeWheights;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector > TensorProductGaussIntegrator::getQuadPntsCoordsPlus1D
                                                              (const vector< RealVector >& nodeCoord)
{
  // helper variable
  const CFuint dimP1 = m_dimensionality + 1;

  // number of nodes
  const CFuint nbrNodes = nodeCoord.size();

  // some checks for correct format
  cf_assert(m_dimensionality != DIM_3D);
  cf_assert(nbrNodes == static_cast<CFuint>(pow(static_cast< CFreal >(2),static_cast< CFreal >(m_dimensionality))));
  cf_assert(nodeCoord[0].size() == dimP1);

  // return variable
  vector< RealVector > qNodeCoord(m_nbrQuadPnts,RealVector(dimP1));

  // compute the quadrature node global coordinates
  if (m_dimensionality != DIM_0D) // not entirely sure if this if() is necessary...
  {
    // compute the quadrature node global coordinates
    for (CFuint iQPnt = 0; iQPnt < m_nbrQuadPnts; ++iQPnt)
    {
      // evaluate geometric basis functions
      vector< CFreal > geoBasisFuncs = evaluateBasisFunction(m_dimensionality,
                                                             CFPolyOrder::ORDER1,
                                                             m_quadPntsMappedCoor[iQPnt]);
      cf_assert(nbrNodes == geoBasisFuncs.size());

      // compute global coordinates
      for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
      {
        cf_assert(nodeCoord[iNode].size() == dimP1);
        qNodeCoord[iQPnt] += nodeCoord[iNode]*geoBasisFuncs[iNode];
      }
    }
  }
  else
  {
    qNodeCoord[0] = nodeCoord[0];
  }

  return qNodeCoord;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector > TensorProductGaussIntegrator::getQuadPntsUnitNormalsPlus1D
                                                                    (const vector< RealVector >& nodeCoord)
{
  // helper variable
  const CFuint dimP1 = m_dimensionality + 1;

  // number of nodes
  const CFuint nbrNodes = nodeCoord.size();

  // some checks for correct format
  cf_assert(m_dimensionality != DIM_3D);
  cf_assert(nbrNodes == static_cast<CFuint>(pow(static_cast< CFreal >(2),static_cast< CFreal >(m_dimensionality))));
  cf_assert(nodeCoord[0].size() == dimP1);

  // return variable
  vector< RealVector > qNodeNormals(m_nbrQuadPnts,RealVector(dimP1));

  // compute jacobian determinant
  switch (m_dimensionality)
  {
    case DIM_0D:
    {
      // compute the quadrature node wheights
      for (CFuint iQNode = 0; iQNode < m_nbrQuadPnts; ++iQNode)
      {
        qNodeNormals[iQNode][XX] = 1.0;
      }
    }
    case DIM_1D:
    {
      // compute the quadrature node wheights
      for (CFuint iQPnt = 0; iQPnt < m_nbrQuadPnts; ++iQPnt)
      {
        // evaluate derivatives of geometrical basis functions
        vector< RealVector > geoBasisFuncDerivs = evaluateDerivBasisFunction(m_dimensionality,
            CFPolyOrder::ORDER1,
            m_quadPntsMappedCoor[iQPnt]);
        cf_assert(nbrNodes == geoBasisFuncDerivs.size());

        // compute local normal vector
        vector< RealVector > vecs(m_dimensionality,RealVector(dimP1));
        for (CFuint iCoor = 0; iCoor < static_cast<CFuint>(m_dimensionality); ++iCoor)
        {
          for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
          {
            vecs[iCoor] += geoBasisFuncDerivs[iNode][iCoor]*nodeCoord[iNode];
          }
        }
        qNodeNormals[iQPnt][XX] =  vecs[KSI][YY];
        qNodeNormals[iQPnt][YY] = -vecs[KSI][XX];

        // divide by length
        qNodeNormals[iQPnt] /= qNodeNormals[iQPnt].norm2();
      }
    } break;
    case DIM_2D:
    {
       // compute the quadrature node wheights
      for (CFuint iQPnt = 0; iQPnt < m_nbrQuadPnts; ++iQPnt)
      {
        // evaluate derivatives of geometrical basis functions
        vector< RealVector > geoBasisFuncDerivs = evaluateDerivBasisFunction(m_dimensionality,
            CFPolyOrder::ORDER1,
            m_quadPntsMappedCoor[iQPnt]);
        cf_assert(nbrNodes == geoBasisFuncDerivs.size());

        // compute local normal vector
        vector< RealVector > vecs(m_dimensionality,RealVector(dimP1));
        for (CFuint iCoor = 0; iCoor < static_cast<CFuint>(m_dimensionality); ++iCoor)
        {
          for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
          {
            vecs[iCoor] += geoBasisFuncDerivs[iNode][iCoor]*nodeCoord[iNode];
          }
        }
        qNodeNormals[iQPnt][XX] = vecs[KSI][YY]*vecs[ETA][ZZ] - vecs[KSI][ZZ]*vecs[ETA][YY];
        qNodeNormals[iQPnt][YY] = vecs[KSI][ZZ]*vecs[ETA][XX] - vecs[KSI][XX]*vecs[ETA][ZZ];
        qNodeNormals[iQPnt][ZZ] = vecs[KSI][XX]*vecs[ETA][YY] - vecs[KSI][YY]*vecs[ETA][XX];

        // divide by length
        qNodeNormals[iQPnt] /= qNodeNormals[iQPnt].norm2();
      }
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Invalid dimensionality");
    }
  }

  return qNodeNormals;
}

//////////////////////////////////////////////////////////////////////////////

vector< CFreal > TensorProductGaussIntegrator::getQuadPntsWheightsPlus1D(const vector< RealVector >& nodeCoord)
{
  // helper variable
  const CFuint dimP1 = m_dimensionality + 1;

  // number of nodes
  const CFuint nbrNodes = nodeCoord.size();

  // some checks for correct format
  cf_assert(m_dimensionality != DIM_3D);
  cf_assert(nbrNodes == static_cast<CFuint>(pow(static_cast< CFreal >(2),static_cast< CFreal >(m_dimensionality))));
  cf_assert(nodeCoord[0].size() == dimP1);

  // return variable
  vector< CFreal > qNodeWheights(m_nbrQuadPnts);

  // compute jacobian determinant
  switch (m_dimensionality)
  {
    case DIM_0D:
    {
      const CFreal jacobDet = 1.;

      // compute the quadrature node wheights
      for (CFuint iQNode = 0; iQNode < m_nbrQuadPnts; ++iQNode)
      {
        qNodeWheights[iQNode] = jacobDet*m_quadPntsWheights[iQNode];
      }
    }
    case DIM_1D:
    {
      // compute the quadrature node wheights
      for (CFuint iQPnt = 0; iQPnt < m_nbrQuadPnts; ++iQPnt)
      {
        // evaluate derivatives of geometrical basis functions
        vector< RealVector > geoBasisFuncDerivs = evaluateDerivBasisFunction(m_dimensionality,
                                                                             CFPolyOrder::ORDER1,
                                                                             m_quadPntsMappedCoor[iQPnt]);
        cf_assert(nbrNodes == geoBasisFuncDerivs.size());

        // compute local normal vector
        vector< RealVector > vecs(m_dimensionality,RealVector(dimP1));
        for (CFuint iCoor = 0; iCoor < static_cast<CFuint>(m_dimensionality); ++iCoor)
        {
          for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
          {
            vecs[iCoor] += geoBasisFuncDerivs[iNode][iCoor]*nodeCoord[iNode];
          }
        }
        RealVector normVect(dimP1);
        normVect[XX] =  vecs[KSI][YY];
        normVect[YY] = -vecs[KSI][XX];

        // normal vector length
        const CFreal normVectLength = sqrt(pow(normVect[XX],2.0) + pow(normVect[YY],2.0));

        // compute quadrature wheight
        qNodeWheights[iQPnt] = normVectLength*m_quadPntsWheights[iQPnt];
      }
    } break;
    case DIM_2D:
    {
       // compute the quadrature node wheights
      for (CFuint iQPnt = 0; iQPnt < m_nbrQuadPnts; ++iQPnt)
      {
        // evaluate derivatives of geometrical basis functions
        vector< RealVector > geoBasisFuncDerivs = evaluateDerivBasisFunction(m_dimensionality,
                                                                             CFPolyOrder::ORDER1,
                                                                             m_quadPntsMappedCoor[iQPnt]);
        cf_assert(nbrNodes == geoBasisFuncDerivs.size());

        // compute local normal vector
        vector< RealVector > vecs(m_dimensionality,RealVector(dimP1));
        for (CFuint iCoor = 0; iCoor < static_cast<CFuint>(m_dimensionality); ++iCoor)
        {
          for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
          {
            vecs[iCoor] += geoBasisFuncDerivs[iNode][iCoor]*nodeCoord[iNode];
          }
        }
        RealVector normVect(dimP1);
        normVect[XX] = vecs[KSI][YY]*vecs[ETA][ZZ] - vecs[KSI][ZZ]*vecs[ETA][YY];
        normVect[YY] = vecs[KSI][ZZ]*vecs[ETA][XX] - vecs[KSI][XX]*vecs[ETA][ZZ];
        normVect[ZZ] = vecs[KSI][XX]*vecs[ETA][YY] - vecs[KSI][YY]*vecs[ETA][XX];

        // normal vector length
        const CFreal normVectLength = sqrt(pow(normVect[XX],2.0) + pow(normVect[YY],2.0) + pow(normVect[ZZ],2.0));

        // compute quadrature wheight
        qNodeWheights[iQPnt] = normVectLength*m_quadPntsWheights[iQPnt];
      }
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Invalid dimensionality");
    }
  }

  return qNodeWheights;
}

//////////////////////////////////////////////////////////////////////////////

vector< CFreal > TensorProductGaussIntegrator::evaluateBasisFunction(const CFDim dimensionality,
                                                                     const CFPolyOrder::Type spatialOrder,
                                                                     const RealVector& localCoord)
{
  switch (spatialOrder)
  {
    case CFPolyOrder::ORDER0:
    case CFPolyOrder::ORDER1:
    {
      switch (dimensionality)
      {
        case DIM_1D:
        {
          cf_assert(localCoord.size() == 1);
          vector< CFreal > evalBaseFunc(2);

          // evaluate the basis functions
          evalBaseFunc[0] = 0.5*(1.0 - localCoord[KSI]);
          evalBaseFunc[1] = 0.5*(1.0 + localCoord[KSI]);

          return evalBaseFunc;
        }
        case DIM_2D:
        {
          cf_assert(localCoord.size() == 2);
          vector< CFreal > evalBaseFunc(4);

          // evaluate the basis functions
          const CFreal ksiNeg = (1.0 - localCoord[KSI]);
          const CFreal ksiPos = (1.0 + localCoord[KSI]);
          const CFreal etaNeg = (1.0 - localCoord[ETA]);
          const CFreal etaPos = (1.0 + localCoord[ETA]);

          evalBaseFunc[0] = 0.25*ksiNeg*etaNeg;
          evalBaseFunc[1] = 0.25*ksiPos*etaNeg;
          evalBaseFunc[2] = 0.25*ksiPos*etaPos;
          evalBaseFunc[3] = 0.25*ksiNeg*etaPos;

          return evalBaseFunc;
        }
        case DIM_3D:
        {
          cf_assert(localCoord.size() == 3);
          vector< CFreal > evalBaseFunc(8);

          // evaluate the basis functions
          const CFreal ksiNeg = (1.0 - localCoord[KSI]);
          const CFreal ksiPos = (1.0 + localCoord[KSI]);
          const CFreal etaNeg = (1.0 - localCoord[ETA]);
          const CFreal etaPos = (1.0 + localCoord[ETA]);
          const CFreal ztaNeg = (1.0 - localCoord[ZTA]);
          const CFreal ztaPos = (1.0 + localCoord[ZTA]);

          evalBaseFunc[0] = 0.125*ksiNeg*etaNeg*ztaNeg;
          evalBaseFunc[1] = 0.125*ksiPos*etaNeg*ztaNeg;
          evalBaseFunc[2] = 0.125*ksiPos*etaPos*ztaNeg;
          evalBaseFunc[3] = 0.125*ksiNeg*etaPos*ztaNeg;
          evalBaseFunc[4] = 0.125*ksiNeg*etaNeg*ztaPos;
          evalBaseFunc[5] = 0.125*ksiPos*etaNeg*ztaPos;
          evalBaseFunc[6] = 0.125*ksiPos*etaPos*ztaPos;
          evalBaseFunc[7] = 0.125*ksiNeg*etaPos*ztaPos;

          return evalBaseFunc;
        }
        default:
        {
          throw Common::ShouldNotBeHereException (FromHere(),"Invalid dimensionality");
        }
      }
    }
    default:
    {
      throw Common::NotImplementedException (FromHere(),"TensorProductGaussIntegrator not implemented for higher spatial order than P1.");
    }
  }

  return vector< CFreal >(0);
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector > TensorProductGaussIntegrator::evaluateDerivBasisFunction(const CFDim dimensionality,
                                                                              const CFPolyOrder::Type spatialOrder,
                                                                              const RealVector& localCoord)
{
  switch (spatialOrder)
  {
    case CFPolyOrder::ORDER0:
    case CFPolyOrder::ORDER1:
    {
      switch (dimensionality)
      {
        case DIM_1D:
        {
          cf_assert(localCoord.size() == 1);
          vector< RealVector > derivBaseFunc(2,RealVector(1));

          // evaluate the basis function derivatives
          derivBaseFunc[0][KSI] = -0.5;
          derivBaseFunc[1][KSI] = +0.5;

          return derivBaseFunc;
        }
        case DIM_2D:
        {
          cf_assert(localCoord.size() == 2);
          vector< RealVector > derivBaseFunc(4,RealVector(2));

          // evaluate the basis function derivatives
          const CFreal ksiNeg = (1.0 - localCoord[KSI]);
          const CFreal ksiPos = (1.0 + localCoord[KSI]);
          const CFreal etaNeg = (1.0 - localCoord[ETA]);
          const CFreal etaPos = (1.0 + localCoord[ETA]);

          derivBaseFunc[0][KSI] = -0.25*etaNeg;
          derivBaseFunc[1][KSI] = +0.25*etaNeg;
          derivBaseFunc[2][KSI] = +0.25*etaPos;
          derivBaseFunc[3][KSI] = -0.25*etaPos;

          derivBaseFunc[0][ETA] = -0.25*ksiNeg;
          derivBaseFunc[1][ETA] = -0.25*ksiPos;
          derivBaseFunc[2][ETA] = +0.25*ksiPos;
          derivBaseFunc[3][ETA] = +0.25*ksiNeg;

          return derivBaseFunc;
        }
        case DIM_3D:
        {
          cf_assert(localCoord.size() == 3);
          vector< RealVector > derivBaseFunc(8,RealVector(3));

          // evaluate the basis function derivatives
          const CFreal ksiNeg = (1.0 - localCoord[KSI]);
          const CFreal ksiPos = (1.0 + localCoord[KSI]);
          const CFreal etaNeg = (1.0 - localCoord[ETA]);
          const CFreal etaPos = (1.0 + localCoord[ETA]);
          const CFreal ztaNeg = (1.0 - localCoord[ZTA]);
          const CFreal ztaPos = (1.0 + localCoord[ZTA]);

          derivBaseFunc[0][KSI] = -0.125*etaNeg*ztaNeg;
          derivBaseFunc[1][KSI] = +0.125*etaNeg*ztaNeg;
          derivBaseFunc[2][KSI] = +0.125*etaPos*ztaNeg;
          derivBaseFunc[3][KSI] = -0.125*etaPos*ztaNeg;
          derivBaseFunc[4][KSI] = -0.125*etaNeg*ztaPos;
          derivBaseFunc[5][KSI] = +0.125*etaNeg*ztaPos;
          derivBaseFunc[6][KSI] = +0.125*etaPos*ztaPos;
          derivBaseFunc[7][KSI] = -0.125*etaPos*ztaPos;

          derivBaseFunc[0][ETA] = -0.125*ksiNeg*ztaNeg;
          derivBaseFunc[1][ETA] = -0.125*ksiPos*ztaNeg;
          derivBaseFunc[2][ETA] = +0.125*ksiPos*ztaNeg;
          derivBaseFunc[3][ETA] = +0.125*ksiNeg*ztaNeg;
          derivBaseFunc[4][ETA] = -0.125*ksiNeg*ztaPos;
          derivBaseFunc[5][ETA] = -0.125*ksiPos*ztaPos;
          derivBaseFunc[6][ETA] = +0.125*ksiPos*ztaPos;
          derivBaseFunc[7][ETA] = +0.125*ksiNeg*ztaPos;

          derivBaseFunc[0][ZTA] = -0.125*ksiNeg*etaNeg;
          derivBaseFunc[1][ZTA] = -0.125*ksiPos*etaNeg;
          derivBaseFunc[2][ZTA] = -0.125*ksiPos*etaPos;
          derivBaseFunc[3][ZTA] = -0.125*ksiNeg*etaPos;
          derivBaseFunc[4][ZTA] = +0.125*ksiNeg*etaNeg;
          derivBaseFunc[5][ZTA] = +0.125*ksiPos*etaNeg;
          derivBaseFunc[6][ZTA] = +0.125*ksiPos*etaPos;
          derivBaseFunc[7][ZTA] = +0.125*ksiNeg*etaPos;

          return derivBaseFunc;
        }
        default:
        {
          throw Common::ShouldNotBeHereException (FromHere(),"Invalid dimensionality");
        }
      }
    }
    default:
    {
      throw Common::NotImplementedException (FromHere(),"TensorProductGaussIntegrator not implemented for higher spatial order than P1.");
    }
  }

  return vector<RealVector>(0);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD
}  // namespace COOLFluiD
