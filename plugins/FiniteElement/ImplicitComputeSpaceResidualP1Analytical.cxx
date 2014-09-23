#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/MeshData.hh"

#include "Heat/HeatPhysicalModel.hh"  // to access _conductivity
#include "Common/NullPointerException.hh"  // guarantee the access

#include "FiniteElement/FiniteElementHeat.hh"
#include "FiniteElement/ImplicitComputeSpaceResidualP1Analytical.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ImplicitComputeSpaceResidualP1Analytical, FiniteElementMethodData, FiniteElementHeatModule> implicitComputeSpaceResidualP1AnalyticalProvider("ImplicitComputeSpaceResP1AnalyticalCom");

//////////////////////////////////////////////////////////////////////////////

ImplicitComputeSpaceResidualP1Analytical::ImplicitComputeSpaceResidualP1Analytical(const std::string& name) :
  FiniteElementMethodCom(name),
  socket_rhs("rhs")
{
}

//////////////////////////////////////////////////////////////////////////////

ImplicitComputeSpaceResidualP1Analytical::~ImplicitComputeSpaceResidualP1Analytical()
{
}

//////////////////////////////////////////////////////////////////////////////

void ImplicitComputeSpaceResidualP1Analytical::setup()
{
  CFAUTOTRACE;

  // first call parent method
  FiniteElementMethodCom::setup();

  // add here specific setup

  // link to Physics::Heat::HeatPhysicalModel to get the conductivity
  try {
    SafePtr< Physics::Heat::HeatPhysicalModel > _heatPhysicalModel =
      PhysicalModelStack::getActive()->getImplementor().
        d_castTo< Physics::Heat::HeatPhysicalModel >();
    _conductivity = _heatPhysicalModel->getConductivity();
  }
  catch (Common::FailedCastException& e) {
    CFLogError("Pointer to Physics::Heat::HeatPhysicalModel: are you actually using the Heat PhysicalModel?");
    exit(1);
  }
  CFLogInfo("Conductivity: " << _conductivity << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ImplicitComputeSpaceResidualP1Analytical::getNodes(
  GeometricEntity& cell,
  CFreal* _x, CFreal* _y, CFreal* _z)
{
  vector<Node*>& nodes = *cell.getNodes();
  for (CFuint i=0; i<4; ++i) {
    Node& node = *nodes[i];
    _x[i] = node[0];
    _y[i] = node[1];
    _z[i] = node[2];
  }
  return;
}

//////////////////////////////////////////////////////////////////////////////

void ImplicitComputeSpaceResidualP1Analytical::getStates(
  GeometricEntity& cell,
  RealVector& _states, CFuint* _statesIDs, bool* _statesIPU)
{
  const CFuint nbEqs = 1;
  vector< State* >& vStates = *cell.getStates();
  for (CFuint iState=0; iState<4; ++iState) {
    State& state = *vStates[iState];
    _statesIDs[iState] = state.getLocalID();
    _statesIPU[iState] = state.isParUpdatable();
    for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
      _states[iState*nbEqs + iEq] = state[iEq];
    }
  }
  return;
}

//////////////////////////////////////////////////////////////////////////////

void ImplicitComputeSpaceResidualP1Analytical::calcVolume(
  const CFreal* x, const CFreal* y, const CFreal* z, CFreal& _V )
{
  _V =
    ((x[1]-x[0])*((y[2]-y[0])*(z[3]-z[0])-(y[3]-y[0])*(z[2]-z[0]))-
     (x[2]-x[0])*((y[1]-y[0])*(z[3]-z[0])-(y[3]-y[0])*(z[1]-z[0]))+
     (x[3]-x[0])*((y[1]-y[0])*(z[2]-z[0])-(y[2]-y[0])*(z[1]-z[0])))/6.;
  return;
}

//////////////////////////////////////////////////////////////////////////////

void ImplicitComputeSpaceResidualP1Analytical::calcTetrahedraKeFe(
  const CFreal V, const CFreal* b, const CFreal* c, const CFreal* d,
  RealMatrix& Ke, RealVector& Fe)
{
  // contribution matrix calculation
  const CFreal t = V * _conductivity;
  CFreal r;
  for (CFuint i=0; i<4; ++i) {
    for (CFuint j=i; j<4; ++j) { // matrix is symetric, start from i
      r = ( b[i]*b[j]+c[i]*c[j]+d[i]*d[j] ) * t;
      Ke(i,j) = r;
      Ke(j,i) = r;
    }
  }
//  Fe source terms?
  return;
}

//////////////////////////////////////////////////////////////////////////////

void ImplicitComputeSpaceResidualP1Analytical::calcTetrahedraShapeFunctions(
  const CFreal& V, const CFreal* x, const CFreal* y, const CFreal* z,
  CFreal* _b, CFreal* _c, CFreal* _d)
{
  // shape functions coefficients
  // (works with absolute and centroid-relative nodal coordinates)
  CFuint j, k, l;
  CFreal xj, yj, zj,  xk, yk, zk,  xl, yl, zl;

  const CFuint jj[4]={1,0,0,0};
  const CFuint kk[4]={3,2,3,1};
  const CFuint ll[4]={2,3,1,2};
  CFreal inv6V = 1./(6.*V);

  for (CFuint i=0;i<4;++i) // loop over nodes of tetrahedra
  {
    j = jj[i];
    k = kk[i];
    l = ll[i];

    xj = x[j];
    xk = x[k];
    xl = x[l];

    yj = y[j];
    yk = y[k];
    yl = y[l];

    zj = z[j];
    zk = z[k];
    zl = z[l];

    _b[i] = ((yk-yj)*(zl-zj) - (yl-yj)*(zk-zj)) * inv6V;
    _c[i] = ((xl-xj)*(zk-zj) - (xk-xj)*(zl-zj)) * inv6V;
    _d[i] = ((xk-xj)*(yl-yj) - (xl-xj)*(yk-yj)) * inv6V;
  }
  return;
}

//////////////////////////////////////////////////////////////////////////////

void ImplicitComputeSpaceResidualP1Analytical::executeOnTrs()
{
  CFAUTOTRACE;

  // rhs storage
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  /// @todo computation of source terms
  /// @todo make it work for multi-scalar problems (with coupling)

  CFreal x[4], y[4], z[4];    // coordinates of nodes
  CFreal b[4], c[4], d[4];    // shape function coefficients
  CFreal V;                   // element volume

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  const CFuint nbGeos = trs->getLocalNbGeoEnts();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = 4;

  // values of the cell states and their IDs
  RealVector states(nbStates*nbEqs);
  CFuint statesID[nbStates];
  bool statesIsParUpdatable[nbStates];

  // block accumulator to jacobian matrix
  BlockAccumulator *acc = getMethodData().getLinearSystemSolver()[0]->
    createBlockAccumulator(nbStates,nbStates,nbEqs);

  // element contributions
  RealMatrix Ke(nbStates*nbEqs,nbStates*nbEqs);
  RealVector Fe(nbStates*nbEqs);
  RealVector Re(nbStates*nbEqs);

  // get system matrix
  SafePtr<LSSMatrix> jacobMatrix =
    getMethodData().getLinearSystemSolver()[0]->getMatrix();

  // provide access to the GeometricEntities
  SafePtr< GeometricEntityPool< StdTrsGeoBuilder > >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;


  // for all the Tetrahedra
  for(CFuint iGeoEnt=0; iGeoEnt<nbGeos; ++iGeoEnt) {

    CFLogDebugMax("Assembling Tetrahedra " << iGeoEnt << "\n");

    // compute the jacobian matrix if it is not frozen
    if(!getMethodData().isSysMatrixFrozen()) {

      // build the GeometricEntity
      geoData.idx = iGeoEnt;
      GeometricEntity& cell = *geoBuilder->buildGE();

      // get nodal coordinates and cell states
      getNodes(cell, x,y,z);
      getStates(cell, states,statesID,statesIsParUpdatable);

      // calculate volume and shape functions
      calcVolume(x,y,z, V);
      calcTetrahedraShapeFunctions(V,x,y,z, b,c,d);

      // calculate tetrahedra contributions to jacobian and RHS vector
      calcTetrahedraKeFe(V,b,c,d, Ke,Fe);
      Re = Ke*states-Fe;

      // for all the states in this cell
      for (CFuint iState=0; iState<nbStates; ++iState) {

        // add the contribution to the RHS vector
        for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
          rhs(statesID[iState], iEq, nbEqs) -=
            Re[iState*nbEqs + iEq];
        }

        // set the index of the block corresponding to the current
        // state in the jacobian matrix
        acc->setRowColIndex(iState, statesID[iState]);
      }

      // add the contribution to the jacobian matrix
      acc->setValuesM(Ke);
      jacobMatrix->addValues( *acc );

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }

  } // for each tetrahedra
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
ImplicitComputeSpaceResidualP1Analytical::needsSockets()
{
  std::vector< SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
