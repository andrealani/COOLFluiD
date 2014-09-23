#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"

#include "HessianEE/HessianEE.hh"
#include "HessianEE/StdSmoothCom.hh"

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

MethodCommandProvider<StdSmoothCom, HessEEData, HessianEEModule> StdSmoothComProvider("StdSmoothCom");

//////////////////////////////////////////////////////////////////////////////

StdSmoothCom::StdSmoothCom(const std::string& name) :
HessEECom(name),
socket_adapt_matrix("adapt_matrix"),
socket_adapt_wght("adapt_wght"),
socket_hessian("hessian"),
socket_metric("metric")
{
}

//////////////////////////////////////////////////////////////////////////////

StdSmoothCom::~StdSmoothCom()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdSmoothCom::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_adapt_matrix);
  result.push_back(&socket_adapt_wght);
  result.push_back(&socket_hessian);
  result.push_back(&socket_metric);

  return result;
}


//////////////////////////////////////////////////////////////////////////////

void StdSmoothCom::setup()
{
  CFAUTOTRACE;

  // first call parent method
  HessEECom::setup();
}

//////////////////////////////////////////////////////////////////////////////

void StdSmoothCom::execute()
{
  CFAUTOTRACE;

  std::string dname = getMethodData().rSmthDataName();
  CFuint niter = getMethodData().rSmthNIter();
  CFreal wght = getMethodData().rSmthWght();

  cout << dname << " : " << niter  << " : " << wght << endl;

  const CFuint size = PhysicalModelStack::getActive()->getDim();

  DataHandle<RealMatrix> adapt_metric = socket_hessian.getDataHandle();
  if(dname == "metric"){
    adapt_metric = socket_metric.getDataHandle();
  }

  DataHandle<CFreal> adapt_wght = socket_adapt_wght.getDataHandle();
  DataHandle<RealMatrix> adapt_matrix = socket_adapt_matrix.getDataHandle();

//      CFout << adapt_matrix.size() << "amatrix\n";
//      CFout << adapt_wght.size() << "wght \n";

  cf_assert( adapt_wght.size() == adapt_matrix.size() );
  cf_assert( adapt_metric.size() == adapt_matrix.size() );

  RealMatrix initMtx( size, size, 0.0);

  SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> > geoBuilder =
    getMethodData().getGeoWithNodesBuilder();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  const CFuint nbCells = trs->getLocalNbGeoEnts();

  for ( CFuint iter=0; iter<niter; ++iter)
  {
    cout << "   smoothing iter " << iter+1 << endl;

    for ( CFuint i=0; i<adapt_metric.size(); ++i)
    {
      adapt_matrix[i] = adapt_metric[i];
      adapt_wght[i] = 1.0;
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

        for ( CFuint i=0; i<tabNod.size(); ++i)
        for ( CFuint j=0; j<tabNod.size(); ++j)
        {
          CFuint idi = tabNod[ i ]->getLocalID();
          CFuint idj = tabNod[ j ]->getLocalID();

                          if ( idi != idj)
                          {
                                  adapt_matrix[idi] = adapt_matrix[idi] + wght * adapt_metric[idj];
                                  adapt_wght[idi] += wght;
                          }
                  }
          }

          for ( CFuint i=0; i<adapt_metric.size(); ++i)
          {
            adapt_metric[i] = 1.0 / adapt_wght[i] * adapt_matrix[i];
          }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD
