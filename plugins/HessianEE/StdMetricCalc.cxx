#include "HessianEE/HessianEE.hh"
#include "StdMetricCalc.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MathTools/MatrixEigenSolver.hh"
#include "MathTools/MatrixInverter.hh"
#include "MathTools/MathChecks.hh"
#include "Framework/VolumeCalculator.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/FilesystemException.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {

  namespace Numerics {

    namespace HessianEE {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdMetricCalc, HessEEData, HessianEEModule> stdMetricCalcProvider("StdMetricCalc");

//////////////////////////////////////////////////////////////////////////////

RealVector getTetraEdge( const CFuint& id, const vector<Node*>& tabNod)
{
  RealVector vedg( 0.0,3);
  switch ( id)
  {
  case 0:
    vedg = (*tabNod[1]->getData());
    vedg -= (*tabNod[0]->getData());
    break;

  case 1:
    vedg = (*tabNod[2]->getData());
    vedg -= (*tabNod[1]->getData());
    break;

  case 2:
    vedg = (*tabNod[0]->getData());
    vedg -= (*tabNod[2]->getData());
    break;

  case 3:
    vedg = (*tabNod[3]->getData());
    vedg -= (*tabNod[0]->getData());
    break;

  case 4:
    vedg = (*tabNod[3]->getData());
    vedg -= (*tabNod[1]->getData());
    break;

  case 5:
    vedg = (*tabNod[3]->getData());
    vedg -= (*tabNod[2]->getData());
    break;

  default:
    throw Common::BadValueException (FromHere(),"Bad id number in getTetraEdge()");
  }

  return vedg;
}

//////////////////////////////////////////////////////////////////////////////

const CFreal ZERO = 1.0e-10;

//////////////////////////////////////////////////////////////////////////////

StdMetricCalc::StdMetricCalc(const std::string& name) :
  HessEECom(name),
socket_metric("metric"),
socket_hessian("hessian"),
socket_grad("grad"),
socket_adapt_func("adapt_func"),
socket_nodes("nodes")
{
}

//////////////////////////////////////////////////////////////////////////////

StdMetricCalc::~StdMetricCalc()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdMetricCalc::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_metric);
  result.push_back(&socket_hessian);
  result.push_back(&socket_grad);
  result.push_back(&socket_adapt_func);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdMetricCalc::setup()
{
  CFAUTOTRACE;

  // first call parent method
  HessEECom::setup();
}

//////////////////////////////////////////////////////////////////////////////

void StdMetricCalc::calcMetric()
{
  CFAUTOTRACE;

  const CFuint dim = PhysicalModelStack::getActive()->getDim();


  DataHandle<RealMatrix> adapt_hess = socket_hessian.getDataHandle();
  DataHandle<RealMatrix> adapt_metric = socket_metric.getDataHandle();

  cf_assert( adapt_hess.size() == adapt_metric.size() );

  MatrixEigenSolver *pEigenSol = MatrixEigenSolver::create( dim, true);
  MatrixInverter *pInverter = MatrixInverter::create( dim, false);

  RealMatrix  mtxL( dim, dim, 0.0), mtxLI( dim, dim, 0.0), mtxTmp( dim, dim, 0.0);
  RealVector  vecD( 0.0, dim);

  CFout << "calc metric started\n";

  // calculating metric
  for (CFuint i = 0; i < adapt_hess.size(); ++i) {
    RealMatrix  mtxA( adapt_hess[i] );

    //    CFout <<"+++++++++++++++++++++++++++++++++++++++\n";
    //    CFout << "mtxA:\n" << mtxA << "\n";

    pEigenSol->eigenCalc( mtxA, mtxL, vecD);
    pInverter->invert( mtxL, mtxLI);

    for ( CFuint k=0; k<dim; ++k)
    {
      vecD[k] = ::std::abs( vecD[k] );
      if ( vecD[k] < ZERO)
        vecD[k] = ZERO;
      //else
      //  vecD[k] = 1.0 / vecD[k];
    }

//     CFout << "vecD: \n" << vecD << "\n";
//     CFout << "mtxL: \n" << mtxL << "\n";
//     CFout << "mtxLI: \n" << mtxLI << "\n";
//
//     CFout << i << "\n";
//     CFout << "metric calculated\n";

//    mtxTmp = mtxL * (vecD * mtxLI);
//    CFout << "vecTmp: \n" << mtxTmp << "\n";

//    CFout << mtxTmp.determ3() <<"\n";;

    adapt_metric[i] = mtxL * (vecD * mtxLI);
  }

  deletePtr(pInverter);
  deletePtr(pEigenSol);
}

//////////////////////////////////////////////////////////////////////////////

CFreal StdMetricCalc::calcConst()
{
  CFAUTOTRACE;

  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  cf_assert( dim == DIM_3D );

  DataHandle<RealMatrix> adapt_hess = socket_hessian.getDataHandle();
  DataHandle<RealMatrix> adapt_metric = socket_metric.getDataHandle();

  cf_assert( adapt_hess.size() == adapt_metric.size() );

  MatrixInverter *pInverter = MatrixInverter::create(dim, false);

  RealMatrix  mtxL( dim, dim, 0.0), mtxLI( dim, dim, 0.0), mtxTmp( dim, dim, 0.0);
  RealVector  vecD( 0.0, dim);
  CFreal  errmax = 0.0, err = 0.0;

  SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> > geoBuilder =
    getMethodData().getGeoWithNodesBuilder();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();
  const CFuint nbCells = trs->getLocalNbGeoEnts();

  // calculation of the constant

  for ( CFuint iCell=0; iCell < nbCells; ++iCell )
  {
    CFLogDebugMax("Cell " << iCell << "\n");

    // set the builder data and create the GeometricEntity (Cell)
    geoData.idx = iCell;
    GeometricEntity* cell = geoBuilder->buildGE();

    /// @todo only tetrahedra for now
    cf_assert( cell->getShape() == CFGeoShape::TETRA );

    vector<Node*>&  tabNod = *(cell->getNodes());

    /// @todo allocating memory all the time is ineficient
    RealMatrix mtxM(dim,dim,0.0);
    for ( CFuint i=0; i<tabNod.size(); ++i)
    {
      //if ( adapt_metric[ (*ptabNod)[i]->getLocalID() ].determ3() > ZERO )
      if ( ! MathChecks::isZero( adapt_metric[ tabNod[i]->getLocalID() ].determ3() ) )
      {
        pInverter->invert( adapt_metric[ tabNod[i]->getLocalID() ], mtxTmp);
        mtxM += mtxTmp;
      }
    }
    mtxM /= tabNod.size();

//      CFout << "iGeoEnt = " << iGeoEnt << "\n";
//      CFout << "mtxM: \n" << mtxM << "\n";
//      CFout << "det = " << mtxM.determ3() << "\n";;

    if ( ! MathChecks::isZero( mtxM.determ3() ) )
    {
      pInverter->invert( mtxM, mtxTmp);
      mtxM = mtxTmp;

      //CFout << "mtxM after invert: \n" << mtxM << "\n";

      CFreal errc = 0.0;

      for ( CFuint k=0; k<6; ++k) // for each edge int tetra
      {
        RealVector  ve = getTetraEdge( k, tabNod);
        CFreal erre = 0.0;

        for ( CFuint ki=0; ki<dim; ++ki)
        {
          CFreal sum = 0.0;

          for ( CFuint kj=0; kj<dim; ++kj)
            sum += ve[kj]*mtxM(ki,kj);

          erre += ve[ki]*sum;
        }

        if ( erre > errc)
          errc = erre;
      }

      err += errc;
      if ( errc > errmax)
        errmax = errc;
    }

    // release the GeometricEntity
    geoBuilder->releaseGE();
  }

  for ( CFuint iCell=0; iCell < nbCells; ++iCell )
  {
    CFLogDebugMax("Cell " << iCell << "\n");

    // set the builder data and create the GeometricEntity (Cell)
    geoData.idx = iCell;
    GeometricEntity* cell = geoBuilder->buildGE();

    /// @todo only tetrahedra for now
    cf_assert( cell->getShape() == CFGeoShape::TETRA );

    vector<Node*>&  tabNod = *(cell->getNodes());

    /// @todo allocating memory all the time is ineficient
    RealMatrix mtxM( dim, dim, 0.0);
    for ( CFuint i=0; i<tabNod.size(); ++i)
    {
      //if ( adapt_metric[ (*ptabNod)[i]->getLocalID() ].determ3() > ZERO )
      if ( ! MathChecks::isZero( adapt_metric[ tabNod[i]->getLocalID() ].determ3() ) )
      {
        pInverter->invert( adapt_metric[ tabNod[i]->getLocalID() ], mtxTmp);
        mtxM += mtxTmp;
      }
    }
    mtxM /= tabNod.size();

//      CFout << "iGeoEnt = " << iGeoEnt << "\n";
//      CFout << "mtxM: \n" << mtxM << "\n";
//      CFout << "det = " << mtxM.determ3() << "\n";;

    if ( ! MathChecks::isZero( mtxM.determ3() ) )
    {
      pInverter->invert( mtxM, mtxTmp);
      mtxM = mtxTmp;

      //CFout << "mtxM after invert: \n" << mtxM << "\n";

      CFreal errc = 0.0;

      for ( CFuint k=0; k<6; ++k) // for each edge int tetra
      {
        RealVector  ve = getTetraEdge( k, tabNod);
        CFreal erre = 0.0;

        for ( CFuint ki=0; ki<dim; ++ki)
        {
          CFreal sum = 0.0;

          for ( CFuint kj=0; kj<dim; ++kj)
            sum += ve[kj]*mtxM(ki,kj);

          erre += ve[ki]*sum;
        }

        if ( erre > errc)
          errc = erre;
      }

      err += errc;
      if ( errc > errmax)
        errmax = errc;
    }

    // release the GeometricEntity
    geoBuilder->releaseGE();
  }

  //  CFout << "ncell = " << ncell << "\n";
  CFout << "err = " << err << "\n";
  err /= nbCells;
  CFLog( WARN, "error av = " << err << " error max = " << errmax << " \n" );

  CFreal  c = getMethodData().getConstant() / err;
  CFLog( WARN, "Constant = " << getMethodData().getConstant() << " c = " << c << " \n" );

  deletePtr( pInverter);

  return c;
}

//////////////////////////////////////////////////////////////////////////////

void StdMetricCalc::limitMetric( const CFreal& c)
{
  CFAUTOTRACE;

  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  DataHandle<RealMatrix> adapt_hess = socket_hessian.getDataHandle();
  DataHandle<RealMatrix> adapt_metric = socket_metric.getDataHandle();

  cf_assert( adapt_hess.size() == adapt_metric.size() );

  MatrixEigenSolver *pEigenSol = MatrixEigenSolver::create( dim, true);
  MatrixInverter *pInverter = MatrixInverter::create( dim, false);

  RealMatrix  mtxL( dim, dim, 0.0), mtxLI( dim, dim, 0.0), mtxTmp( dim, dim, 0.0);
  RealVector  vecD( 0.0, dim);

  // resizing and limiting of the metric
  CFreal hmin_inv = 1.0 / (getMethodData().getMetricMinLimit() * getMethodData().getMetricMinLimit());
  CFreal hmax_inv = 1.0 / (getMethodData().getMetricMaxLimit() * getMethodData().getMetricMaxLimit());
  CFreal armax = getMethodData().getMetricMaxAR();
  CFreal lmax;

  if ( std::abs( getMethodData().getMetricMinLimit() ) < ZERO)
    hmin_inv = 1.0 / ZERO;


  CFLog( WARN, "MetricMinLimit = " << getMethodData().getMetricMinLimit() << " hmin_inv = " << hmin_inv << " \n" );
  CFLog( WARN, "MetricMaxLimit = " << getMethodData().getMetricMaxLimit() << " hmax_inv = " << hmax_inv << " \n" );
  CFLog( WARN, "MetricMaxAR = " << getMethodData().getMetricMaxAR() << " \n" );

  CFout << "MaxAR = " << armax << "\n";
  for ( CFuint i=0; i<adapt_metric.size(); ++i)
  {
    //pInverter->invert( adapt_metric[i], mtxTmp);
    //adapt_metric[i] = c * mtxTmp;

    mtxTmp = c * adapt_metric[i];
    pEigenSol->eigenCalc( mtxTmp, mtxL, vecD);
    pInverter->invert( mtxL, mtxLI);

    // limiting max and min size of eigenvalues
    for ( CFuint k=0; k<dim; ++k)
    {
      if ( vecD[k] > hmin_inv )
        vecD[k] = hmin_inv;

      if ( vecD[k] < hmax_inv )
        vecD[k] = hmax_inv;
    }

    // limiting aspect ratio of a metric
    if ( armax > 0)
    {
      lmax = vecD.max();
      //lmax = max( vecD[0], max( vecD[1], vecD[2] ) );
      lmax /= armax * armax;

      for ( CFuint k=0; k<dim; ++k)
      {
        if ( vecD[k] < lmax)
          vecD[k] = lmax;
      }
    }

    adapt_metric[i] = mtxL * (vecD * mtxLI);
  }

  deletePtr( pInverter);
  deletePtr( pEigenSol);
}

//////////////////////////////////////////////////////////////////////////////

void StdMetricCalc::execute()
{
  CFAUTOTRACE;

  calcMetric();

  CFreal c = calcConst();
  CFout << "  CONST c = " << c << "\n";

  limitMetric( c);

  writeTEC( "metric.plt" );
}

//////////////////////////////////////////////////////////////////////////////

void StdMetricCalc::writeTEC( const std::string& fname)
{
  CFAUTOTRACE;

  FILE* f = fopen( fname.c_str(), "wt");
  if(!f) throw Common::FilesystemException (FromHere(),"Cannot open file " + fname);

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<RealVector> adapt_grad = socket_grad.getDataHandle();
  DataHandle<CFreal> adapt_func = socket_adapt_func.getDataHandle();
  DataHandle<RealMatrix> adapt_hess = socket_hessian.getDataHandle();
  DataHandle<RealMatrix> adapt_metric = socket_metric.getDataHandle();

  SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> > geoBuilder =
    getMethodData().getGeoWithNodesBuilder();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();
  const CFuint nbCells = trs->getLocalNbGeoEnts();

  cf_assert( nodes.size() == adapt_func.size() );

  const CFuint dim = DIM_3D;
  MatrixEigenSolver *pEigenSol = MatrixEigenSolver::create(dim, true);
  RealMatrix  mtxL(dim, dim, 0.0), mtxTmp(dim, dim, 0.0);
  RealVector  vecD(0.0, dim);

  fprintf( f, "VARIABLES = \"X\", \"Y\", \"Z\", \"F\", \"GX\", \"GY\", \"GZ\""
              ", \"HXX\", \"HXY\", \"HXZ\", \"HYX\", \"HYY\", \"HYZ\", \"HZX\", \"HZY\", \"HZZ\""
              ", \"EX\", \"EY\", \"EZ\"\n" );
  fprintf( f, "ZONE T=\"grid\", N=%d, E=%d, F=FEPOINT, ET=BRICK C=BLACK \n", nodes.size(), nbCells );

  // nodes
  for ( CFuint ist=0; ist<nodes.size(); ++ist)
  {
    mtxTmp = adapt_metric[ist];
    pEigenSol->eigenCalc( mtxTmp, mtxL, vecD);

    fprintf( f, "%16.10e ", (*nodes[ist])[0] );
    fprintf( f, "%16.10e ", (*nodes[ist])[1] );
    fprintf( f, "%16.10e ", (*nodes[ist])[2] );
    fprintf( f, "%16.10e ", adapt_func[ist] );
    fprintf( f, "%16.10e ", (adapt_grad[ist])[0] );
    fprintf( f, "%16.10e ", (adapt_grad[ist])[1] );
    fprintf( f, "%16.10e ", (adapt_grad[ist])[2] );
    for ( CFuint ih=0; ih<3; ++ih)
      for ( CFuint jh=0; jh<3; ++jh)
        fprintf( f, "%16.10e ", (adapt_metric[ist])(ih,jh) );

    for ( CFuint ih=0; ih<dim; ++ih)
        fprintf( f, "%16.10e ", 1.0/ sqrt( vecD[ih] ) );

    fprintf( f, "\n");
  }

  // cell conectivity in BRICK form
  for ( CFuint iCell=0; iCell < nbCells; ++iCell )  {

    CFLogDebugMax("Cell " << iCell << "\n");

    // set the builder data and create the GeometricEntity (Cell)
    geoData.idx = iCell;
    GeometricEntity* cell = geoBuilder->buildGE();

    /// @todo only tetrahedra for now
    cf_assert( cell->getShape() == CFGeoShape::TETRA );

    vector<Node*>&  tabNod = *(cell->getNodes());

    fprintf( f, "%d ", tabNod[0]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[1]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[2]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[2]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[3]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[3]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[3]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[3]->getLocalID() + 1);

    fprintf( f, "\n");

    // release the GeometricEntity
    geoBuilder->releaseGE();
  }

  fclose( f);

  deletePtr( pEigenSol);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD
