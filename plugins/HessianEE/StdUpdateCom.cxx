#include "HessianEE/HessianEE.hh"
#include "HessianEE/StdUpdateCom.hh"

#include "Common/FilesystemException.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MathTools/MatrixEigenSolver.hh"
#include "MathTools/MatrixInverter.hh"
#include "MathTools/MatrixIntersect.hh"
#include "MathTools/MathChecks.hh"

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

MethodCommandProvider<StdUpdateCom, HessEEData, HessianEEModule> stdUpdateComProvider("StdUpdateCom");

//////////////////////////////////////////////////////////////////////////////

StdUpdateCom::StdUpdateCom(const std::string& name) :
HessEECom(name),
  socket_glob_metric("glob_metric"),
  socket_metric("metric"),
  socket_nodes("nodes")
{
}

//////////////////////////////////////////////////////////////////////////////

StdUpdateCom::~StdUpdateCom()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdUpdateCom::setup()
{
  CFAUTOTRACE;

  // first call parent method
  HessEECom::setup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdUpdateCom::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_glob_metric);
  result.push_back(&socket_metric);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdUpdateCom::execute()
{
  CFAUTOTRACE;

  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  DataHandle<RealMatrix> glob_metric = socket_glob_metric.getDataHandle();
  DataHandle<RealMatrix> adapt_metric = socket_metric.getDataHandle();

  cf_assert( glob_metric.size() == adapt_metric.size() );

  MatrixIntersect *pIntersect = MatrixIntersect::create( dim);


  for ( CFuint i=0; i<adapt_metric.size(); ++i)
  {
    RealMatrix  mtxA( adapt_metric[i] );
    RealMatrix  mtxB( glob_metric[i] );
                RealMatrix  mtxC( dim, dim, 0.0);

                pIntersect->intersectCalc( mtxA, mtxB, mtxC);

                glob_metric[i] = mtxC;
  }

  writeTEC( "global_metric.plt");

  deletePtr( pIntersect);
}

//////////////////////////////////////////////////////////////////////////////

void StdUpdateCom::writeTEC( const std::string& fname)
{
  FILE* f = fopen( fname.c_str(), "wt");
  if(!f) throw Common::FilesystemException (FromHere(),"Cannot open file " + fname);

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<RealMatrix> adapt_metric = socket_glob_metric.getDataHandle();

  CFuint dim = 3;
  MatrixEigenSolver *pEigenSol = MatrixEigenSolver::create( dim, true);
  RealMatrix  mtxL( dim, dim, 0.0), mtxTmp( dim, dim, 0.0);
  RealVector  vecD( 0.0, dim);

  SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> > geoBuilder =
    getMethodData().getGeoWithNodesBuilder();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  const CFuint nbCells = trs->getLocalNbGeoEnts();

  fprintf( f, "VARIABLES = \"X\", \"Y\", \"Z\""
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

    for ( CFuint ih=0; ih<3; ++ih)
      for ( CFuint jh=0; jh<3; ++jh)
        fprintf( f, "%16.10e ", (adapt_metric[ist])(ih,jh) );

    for ( CFuint ih=0; ih<dim; ++ih)
        fprintf( f, "%16.10e ", 1.0/ sqrt( vecD[ih] ) );

    fprintf( f, "\n");
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

    fprintf( f, "%d ", tabNod[0]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[1]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[2]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[2]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[3]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[3]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[3]->getLocalID() + 1);
    fprintf( f, "%d ", tabNod[3]->getLocalID() + 1);

    fprintf( f, "\n");

    geoBuilder->releaseGE();
  }

  fclose( f);

  deletePtr( pEigenSol);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD
