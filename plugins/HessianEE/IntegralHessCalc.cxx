#include "HessianEE/HessianEE.hh"
#include "HessianEE/IntegralHessCalc.hh"

#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/VolumeCalculator.hh"
#include "Common/FilesystemException.hh"
#include "Common/BadValueException.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

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

MethodCommandProvider<IntegralHessCalc, HessEEData, HessianEEModule> integralHessCalcProvider("IntegralHessCalc");

//////////////////////////////////////////////////////////////////////////////

vector<CFuint> getTetraOppositeNodesIds( const CFuint& ibase)
{
  vector<CFuint> tab(3,0);

  switch (ibase)
  {
  case 0:
    tab[0] = 1;  tab[1] = 2;  tab[2] = 3;
    break;

  case 1:
    tab[0] = 0;  tab[1] = 3;  tab[2] = 2;
    break;

  case 2:
    tab[0] = 0;  tab[1] = 1;  tab[2] = 3;
    break;

  case 3:
    tab[0] = 0;  tab[1] = 2;  tab[2] = 1;
    break;

  default:
    throw Common::BadValueException (FromHere(),"Bad number of ibase in getTetraOppositeNodes()");
  }

  return tab;
}

//////////////////////////////////////////////////////////////////////////////

RealVector getTetraOppositeVn(const CFuint& ibase, const vector<Node*>& tabNod)
{
  vector<CFuint>  tabId = getTetraOppositeNodesIds( ibase);

  //CFout << "ibase: " << ibase << "\n";
  //for ( CFuint i=0; i<tabId.size(); ++i)
  //  CFout << tabId[i] << " ";
  //CFout << "\n";

  //CFout << "qq: " << tabNod.size() << "\n";
  //CFout << "n0: " << (*tabNod[0]->getData()) << "\n";

  RealVector v10 = (*tabNod[tabId[1]]->getData());
  RealVector v20 = (*tabNod[tabId[2]]->getData());
  RealVector vret(0,3);

  v10 -= (*tabNod[tabId[0]]->getData());
  v20 -= (*tabNod[tabId[0]]->getData());

  //CFout << "n0: " << (tabNod[tabId[0]]->getData())->size() << "\n";
  //CFout << "n1: " << (tabNod[tabId[1]]->getData())->size() << "\n";
  //CFout << "n2: " << (tabNod[tabId[2]]->getData())->size() << "\n";
  //CFout << "v10: " << v10 << "\n";
  //CFout << "v20: " << v20 << "\n";

  MathTools::MathFunctions::crossProd( v10, v20, vret);

  //CFout << "vret: " << vret << "\n";

  vret *= 0.5;
  return  vret;
}

//////////////////////////////////////////////////////////////////////////////

void IntegralHessCalc::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("RefH","Reference H.");
}

//////////////////////////////////////////////////////////////////////////////

IntegralHessCalc::IntegralHessCalc(const std::string& name) :
HessEECom(name),
socket_grad("grad"),
socket_adapt_func("adapt_func"),
socket_adapt_hess("hessian"),
socket_adapt_wght("adapt_wght"),
socket_nodes("nodes")
{
   addConfigOptionsTo(this);
  _refH = 10.0;
   setParameter("RefH",&_refH);
}

//////////////////////////////////////////////////////////////////////////////

IntegralHessCalc::~IntegralHessCalc()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
IntegralHessCalc::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_grad);
  result.push_back(&socket_adapt_func);
  result.push_back(&socket_adapt_hess);
  result.push_back(&socket_adapt_wght);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void IntegralHessCalc::setup()
{
  CFAUTOTRACE;

  // first call parent method
  HessEECom::setup();
}

//////////////////////////////////////////////////////////////////////////////

void IntegralHessCalc::calcGrad()
{
  CFAUTOTRACE;

  DataHandle<RealVector> adapt_grad = socket_grad.getDataHandle();
  DataHandle<CFreal> adapt_wght = socket_adapt_wght.getDataHandle();
  DataHandle<CFreal> adapt_func = socket_adapt_func.getDataHandle();

  cf_assert( adapt_func.size() == adapt_grad.size() && adapt_func.size() == adapt_wght.size() );

  // initializing data
  for (CFuint i=0; i<adapt_grad.size(); ++i)
  {
    adapt_grad[i] = RealVector(0.0, 3);
    adapt_wght[i] = 0.0;
  }

  // calc gradient using boundary integral
  RealVector  vtmp(0,3), vn(0,3);
  CFreal sumvol = 0.0;

  SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> > geoBuilder =
    getMethodData().getGeoWithNodesBuilder();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();
  const CFuint nbCells = trs->getLocalNbGeoEnts();
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {

    CFLogDebugMax("Cell " << iCell << "\n");

    // set the builder data and create the GeometricEntity (Cell)
    geoData.idx = iCell;
    GeometricEntity* cell = geoBuilder->buildGE();
    const CFreal volume = cell->computeVolume();
    sumvol += volume;

    vector<Node*>&  tabNod = *(cell->getNodes());

    RealVector vgrad(0,3);

    for (CFuint i=0; i<tabNod.size(); ++i)
    {
      vector<CFuint> tabId = getTetraOppositeNodesIds( i);
      vn = getTetraOppositeVn( i, tabNod );
      CFreal  avu = 0.0;

      for ( CFuint k=0; k<3; ++k)
      {
        CFuint id = tabNod[ tabId[k] ]->getLocalID();
        avu += adapt_func[id];
      }
      avu /= 3.0;

      vtmp = avu * vn;
      vgrad += vtmp;
    }

    for ( CFuint i=0; i<tabNod.size(); ++i)
    {
      adapt_grad[ tabNod[i]->getLocalID() ] += vgrad;
      adapt_wght[ tabNod[i]->getLocalID() ] += volume;
    }

    // release the GeometricEntity
    geoBuilder->releaseGE();
  }


  CFLog( WARN, "Volume = " << sumvol << "\n" );

  for ( CFuint i=0; i<adapt_grad.size(); ++i)
  {
    //cout << i << "   : " << endl << adapt_grad[i] << endl << endl;
    adapt_grad[i] /= adapt_wght[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

void IntegralHessCalc::calcHess()
{
  CFAUTOTRACE;

  cf_assert(PhysicalModelStack::getActive()->getDim() == DIM_3D);

  DataHandle<RealVector> adapt_grad = socket_grad.getDataHandle();
  DataHandle<CFreal> adapt_wght = socket_adapt_wght.getDataHandle();
  DataHandle<RealMatrix> adapt_hess = socket_adapt_hess.getDataHandle();

  SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> > geoBuilder =
    getMethodData().getGeoWithNodesBuilder();

  cf_assert( adapt_hess.size() == adapt_grad.size() && adapt_hess.size() == adapt_wght.size() );

  // initializing data
  for ( CFuint i=0; i<adapt_grad.size(); ++i)
  {
    adapt_hess[i] = RealMatrix( 3, 3, 0.0);
    adapt_wght[i] = 0.0;
  }

  // calc hessian using boundary integral
  RealVector  vtmp(0,3), vn(0,3);
  CFreal sumvol = 0.0;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();
  const CFuint nbCells = trs->getLocalNbGeoEnts();
  for ( CFuint iCell=0; iCell < nbCells; ++iCell )  {

    CFLogDebugMax("Cell " << iCell << "\n");

    // set the builder data and create the GeometricEntity (Cell)
    geoData.idx = iCell;
    GeometricEntity* cell = geoBuilder->buildGE();
    const CFreal volume = cell->computeVolume();
    sumvol += volume;

    vector<Node*>&  tabNod = *(cell->getNodes());

    /// @todo only tetrahedra for now
    cf_assert( cell->getShape() == CFGeoShape::TETRA );

    RealVector vgradx(0,3), vgrady(0,3), vgradz(0,3);

    for ( CFuint i=0; i<tabNod.size(); ++i)
    {
      vector<CFuint> tabId = getTetraOppositeNodesIds( i);
      vn = getTetraOppositeVn( i, tabNod );
      CFreal  avxu = 0.0;
      CFreal  avyu = 0.0;
      CFreal  avzu = 0.0;

      for ( CFuint k=0; k<3; ++k)
      {
        CFuint id = tabNod[ tabId[k] ]->getLocalID();
        avxu += (adapt_grad[id][0]);
        avyu += (adapt_grad[id][1]);
        avzu += (adapt_grad[id][2]);
      }
      avxu /= 3.0;
      avyu /= 3.0;
      avzu /= 3.0;

      vtmp = avxu * vn;   vgradx += vtmp;
      vtmp = avyu * vn;   vgrady += vtmp;
      vtmp = avzu * vn;   vgradz += vtmp;
    }

    for ( CFuint i=0; i<tabNod.size(); ++i)
    {
      CFuint locId = tabNod[i]->getLocalID();

      adapt_hess[ locId ](0,0) += vgradx[0];
      adapt_hess[ locId ](0,1) += 0.5* ( vgradx[1] + vgrady[0] );
      adapt_hess[ locId ](0,2) += 0.5* ( vgradx[2] + vgradz[0] );

      adapt_hess[ locId ](1,0) += 0.5* ( vgrady[0] + vgradx[1] );
      adapt_hess[ locId ](1,1) += vgrady[1];
      adapt_hess[ locId ](1,2) += 0.5* ( vgrady[2] + vgradz[1] );

      adapt_hess[ locId ](2,0) += 0.5* ( vgradz[0] + vgradx[2] );
      adapt_hess[ locId ](2,1) += 0.5* ( vgradz[1] + vgrady[2] );
      adapt_hess[ locId ](2,2) += vgradz[2];

      adapt_wght[ locId ] += volume;
    }

    // release the GeometricEntity
    geoBuilder->releaseGE();
  }

  CFLog( WARN, "Volume = " << sumvol << "\n" );

  for ( CFuint i=0; i<adapt_hess.size(); ++i) {
    adapt_hess[i] /= adapt_wght[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

void IntegralHessCalc::writeTEC( const std::string& fname)
{
  CFAUTOTRACE;

  FILE* f = fopen( fname.c_str(), "wt");
  if(!f) throw Common::FilesystemException (FromHere(),"Cannot open file " + fname);

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<RealVector> adapt_grad = socket_grad.getDataHandle();
  DataHandle<CFreal> adapt_wght = socket_adapt_wght.getDataHandle();
  DataHandle<RealMatrix> adapt_hess = socket_adapt_hess.getDataHandle();
  DataHandle<CFreal> adapt_func = socket_adapt_func.getDataHandle();

  SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> > geoBuilder =
    getMethodData().getGeoWithNodesBuilder();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();
  const CFuint nbCells = trs->getLocalNbGeoEnts();

  cf_assert( nodes.size() == adapt_func.size() );

  fprintf( f, "VARIABLES = \"X\", \"Y\", \"Z\", \"F\", "
              "\"GX\", \"GY\", \"GZ\", \"HXX\", \"HXY\", "
              "\"HXZ\", \"HYX\", \"HYY\", \"HYZ\", \"HZX\", "
              "\"HZY\", \"HZZ\"\n" );
  fprintf( f,
           "ZONE T=\"grid\", N=%d, E=%d, F=FEPOINT, ET=BRICK C=BLACK \n",
           nodes.size(),
           nbCells );

  // nodes
  for (CFuint ist=0; ist<nodes.size(); ++ist)
  {
    fprintf( f, "%16.10e ", (*nodes[ist])[0] );
    fprintf( f, "%16.10e ", (*nodes[ist])[1] );
    fprintf( f, "%16.10e ", (*nodes[ist])[2] );
    fprintf( f, "%16.10e ", adapt_func[ist] );
    fprintf( f, "%16.10e ", (adapt_grad[ist])[0] );
    fprintf( f, "%16.10e ", (adapt_grad[ist])[1] );
    fprintf( f, "%16.10e ", (adapt_grad[ist])[2] );
    for ( CFuint ih=0; ih<3; ++ih)
      for ( CFuint jh=0; jh<3; ++jh)
        fprintf( f, "%16.10e ", (adapt_hess[ist])(ih,jh) );

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

#if 0
  for ( CFuint iTR = 0; iTR < trs->getNbTRs(); ++iTR)
  {
    TopologicalRegion* tr = trs->getTopologicalRegion(iTR);
    for ( CFuint iGeoEnt = 0; iGeoEnt < tr->getNbGeomEnt(); ++iGeoEnt)
    {
      GeometricEntity* cell = tr->getGeomEnt(iGeoEnt);

      cf_assert( cell->nbNodes() == cell->nbStates() );
      cf_assert( cell->getShape() == TETRA );

      vector<Node*>* ptabNod = cell->getNodes();

      fprintf( f, "%d ", (*ptabNod)[0]->getLocalID() + 1);
      fprintf( f, "%d ", (*ptabNod)[1]->getLocalID() + 1);
      fprintf( f, "%d ", (*ptabNod)[2]->getLocalID() + 1);
      fprintf( f, "%d ", (*ptabNod)[2]->getLocalID() + 1);
      fprintf( f, "%d ", (*ptabNod)[3]->getLocalID() + 1);
      fprintf( f, "%d ", (*ptabNod)[3]->getLocalID() + 1);
      fprintf( f, "%d ", (*ptabNod)[3]->getLocalID() + 1);
      fprintf( f, "%d ", (*ptabNod)[3]->getLocalID() + 1);

      fprintf( f, "\n");
    }
  }
#endif

  fclose(f);
}


//////////////////////////////////////////////////////////////////////////////

void IntegralHessCalc::execute()
{
  CFAUTOTRACE;

  calcGrad();
  CFLog(NOTICE,"After calcGrad\n");
  calcHess();
  CFLog(NOTICE,"After calcHess\n");

  writeTEC( "grad.plt");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD
