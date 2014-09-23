#include "RemeshMeandr/WriteMetricCom.hh"
#include "RemeshMeandr/RemeshMeandr.hh"

#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Environment/DirPaths.hh"

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

    namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WriteMetricCom, RMeshMeData, RemeshMeandrModule> WriteMetricComProvider("StdWriteControlSpc");

//////////////////////////////////////////////////////////////////////////////

WriteMetricCom::WriteMetricCom(const std::string& name) :
  RMeshMeCom(name),
  socket_nodes("nodes"),
  socket_glob_metric("glob_metric")
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

WriteMetricCom::~WriteMetricCom()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteMetricCom::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_glob_metric);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WriteMetricCom::setup()
{
  CFAUTOTRACE;

  // first call parent method
  RMeshMeCom::setup();
}

//////////////////////////////////////////////////////////////////////////////

void WriteMetricCom::execute()
{
  CFTRACEBEGIN;

  CFLog( ERROR, "WriteMetricCom::executeOnTrs() BEGIN\n" );

	WriteControlSpc( "output");

  CFLog( ERROR, "WriteMetricCom::executeOnTrs() END\n" );

  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void WriteMetricCom::WriteControlSpc( const std::string& fname)
{
  const CFuint nDim = PhysicalModelStack::getActive()->getDim();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<RealMatrix> adapt_metric = socket_glob_metric.getDataHandle();

  std::string resultsDir = Environment::DirPaths::getInstance().getResultsDir().string();
  std::string workingDir = Environment::DirPaths::getInstance().getWorkingDir().string();
  std::string meDir = Environment::DirPaths::getInstance().getWorkingDir().string() + getMethodData().getMeandrosDir();

  std::string fullname;

  const CFuint size = nDim;

  RealMatrix  mtx( size, size, 0.0);
  //char	sbuf[512];

  FILE *f;

  cout << "WriteControlSpc\n";

  fullname = meDir + fname + ".metric";

  f = fopen( fullname.c_str(), "wt");
  cf_assert( f);	// <----

  fprintf( f, "%d %d\n", nodes.size(), 3);

  // writing metric file
  for ( CFuint ist=0; ist<nodes.size(); ++ist)
  {
	mtx = adapt_metric[ist];

	fprintf( f, "%d %16.10e %16.10e %16.10e %16.10e %16.10e %16.10e\n",	ist+1, mtx(0,0), mtx(0,1), mtx(0,2), mtx(1,1), mtx(1,2), mtx(2,2) );
  }

  fclose( f);


  fullname = meDir + fname + ".node";

  f = fopen( fullname.c_str(), "wt");
  cf_assert( f);	// <----

  fprintf( f, "%d %d %d %d\n", nodes.size(), 3, 0, 1);
  for ( CFuint k=0; k<nodes.size(); ++k)
  	fprintf( f, "%d %16.10e %16.10e %16.10e %d\n", k+1, (*nodes[k])[0], (*nodes[k])[1], (*nodes[k])[2], 0 );

  fclose( f);
}


//////////////////////////////////////////////////////////////////////////////
//
void WriteMetricCom::writeTEC( const std::string& fname)
{
CF_DEBUG_STR("This function needs to be updated not to use Elements");
CF_DEBUG_ABORT;
#if 0
  FILE* f = fopen( fname.c_str(), "wt");
  if ( ! f)
  {
    CFLog( ERROR, "Could not open file \"" << fname << "\" \n" );
    cf_assert( 0);
  }

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = MeshDataStack::getActive()->getDataStorage()->getData<Node*>("Default_nodes");
  DataHandle<Element*> elements = MeshDataStack::getActive()->getDataStorage()->getData<Element*>("elements");

  DataHandle<RealMatrix> adapt_metric = MeshDataStack::getActive()->getDataStorage()->getData<RealMatrix>("Default_glob_metric");

  CFuint size = 3;
  MatrixEigenSolver *pEigenSol = MatrixEigenSolver::create( size, true);
  RealMatrix  mtxL( size, size, 0.0), mtxTmp( size, size, 0.0);
  RealVector  vecD( 0.0, size);

  CFuint csize = elements.size();

  fprintf( f, "VARIABLES = \"X\", \"Y\", \"Z\""
      ", \"HXX\", \"HXY\", \"HXZ\", \"HYX\", \"HYY\", \"HYZ\", \"HZX\", \"HZY\", \"HZZ\""
      ", \"EX\", \"EY\", \"EZ\"\n" );
  fprintf( f, "ZONE T=\"grid\", N=%d, E=%d, F=FEPOINT, ET=BRICK C=BLACK \n", nodes.size(), csize );

  // nodes
  for ( CFuint ist=0; ist<nodes.size(); ++ist)
  {
    mtxTmp = adapt_metric[ist];
    pEigenSol->eigenCalc( mtxTmp, mtxL, vecD);

    fprintf( f, "%16.10l", (*nodes[ist])[0] );
    fprintf( f, "%16.10l", (*nodes[ist])[1] );
    fprintf( f, "%16.10l", (*nodes[ist])[2] );

    for ( CFuint ih=0; ih<3; ++ih)
      for ( CFuint jh=0; jh<3; ++jh)
        fprintf( f, "%16.10l", (adapt_metric[ist])(ih,jh) );

    for ( CFuint ih=0; ih<size; ++ih)
        fprintf( f, "%16.10l", 1.0/ sqrt( vecD[ih] ) );

    fprintf( f, "\n");
  }

  // cell conectivity in BRICK form
  for ( CFuint iCell=0; iCell<elements.size(); ++iCell )
  {
    CFLogDebugMax("Cell " << iCell << "\n");

		Element* elem = elements[iCell];

		vector<Node*> tabNod( elem->getNodes()->size() );

		for ( CFuint i=0; i<elem->getNodes()->size(); ++i)
			tabNod[i] = elem->getNode( i);

		fprintf( f, "%d ", tabNod[0]->getLocalID() + 1);
		fprintf( f, "%d ", tabNod[1]->getLocalID() + 1);
		fprintf( f, "%d ", tabNod[2]->getLocalID() + 1);
		fprintf( f, "%d ", tabNod[2]->getLocalID() + 1);
		fprintf( f, "%d ", tabNod[3]->getLocalID() + 1);
		fprintf( f, "%d ", tabNod[3]->getLocalID() + 1);
		fprintf( f, "%d ", tabNod[3]->getLocalID() + 1);
		fprintf( f, "%d ", tabNod[3]->getLocalID() + 1);

		fprintf( f, "\n");

  }


  fclose( f);

  deletePtr( pEigenSol);
#endif
}

//////////////////////////////////////////////////////////////////////////////



    } // namespace RemeshMeandros

  } // namespace Numerics

} // namespace COOLFluiD
