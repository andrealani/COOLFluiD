#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"      // create a block accumulator
#include "Framework/LSSMatrix.hh"             // access the system matrix
#include "Framework/SubSystemStatus.hh"       // access the simulation time

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/NeumannBCImplicitP1Analytical.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NeumannBCImplicitP1Analytical, FiniteElementMethodData, FiniteElementModule> NeumannBCImplicitP1AnalyticalProvider("NeumannBCImplicitP1Analytical");

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicitP1Analytical::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< bool >("Average","Use element-averaged boundary condition contribution as opposed to nodal contribution (default).");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

NeumannBCImplicitP1Analytical::NeumannBCImplicitP1Analytical(const std::string& name) :
  FiniteElementMethodCom(name),
  socket_rhs("rhs")
{
   addConfigOptionsTo(this);
  _functions = std::vector<std::string>();
   setParameter("Def",&_functions);



  _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);



  _useAverage = false;
   setParameter("Average",&_useAverage);



}

//////////////////////////////////////////////////////////////////////////////

NeumannBCImplicitP1Analytical::~NeumannBCImplicitP1Analytical()
{
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicitP1Analytical::setup()
{
  // first call parent method
  FiniteElementMethodCom::setup();
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicitP1Analytical::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicitP1Analytical::getNodes(
  GeometricEntity& cell, CFuint nbNodes,
  CFreal* _x, CFreal* _y, CFreal* _z)
{
  vector<Node*>& nodes = *cell.getNodes();
  for (CFuint i=0;i<nbNodes;i++) {
    Node& node = *nodes[i];
    _x[i] = node[0];
    _y[i] = node[1];
    _z[i] = node[2];
  }
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicitP1Analytical::getStates(
  GeometricEntity& cell, CFuint nbStates,
  RealVector& _states, CFuint* _statesIDs, bool* _statesIPU)
{
  const CFuint nbEqs = 1;
  vector<State*>& vStates = *cell.getStates();
  for (CFuint iState=0; iState<nbStates; iState++) {
    State& state = *vStates[iState];
    _statesIDs[iState] = state.getLocalID();
    _statesIPU[iState] = state.isParUpdatable();
    for (CFuint iEq=0; iEq<nbEqs; iEq++) {
      _states[iState*nbEqs + iEq] = state[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicitP1Analytical::calcArea(
    const CFreal* x, const CFreal* y, const CFreal* z,
    CFreal& A)
{
  CFreal a = ( (y[1]-y[0]) * (z[2]-z[0]) - (z[1]-z[0]) * (y[2]-y[0]) );
  CFreal b = ( (z[1]-z[0]) * (x[2]-x[0]) - (x[1]-x[0]) * (z[2]-z[0]) );
  CFreal c = ( (x[1]-x[0]) * (y[2]-y[0]) - (y[1]-y[0]) * (x[2]-x[0]) );
  A = sqrt(a*a + b*b + c*c)/2.;
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicitP1Analytical::executeOnTrs()
{
  CFAUTOTRACE;

  if(getMethodData().isDirichletBCApplied())
  {
    CFLog(WARN," #!!!!!!!!!!!!!!!\n You are applying a NeumannBC after a DirichletBC...\n expect strange results!!!\n #!!!!!!!!!!!!!!!\n");
  }

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  /// @todo work for multi-scalar problems (with coupling)

  // PhysicalModel properties
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim   = PhysicalModelStack::getActive()->getDim();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  const CFuint nbGeos = trs->getLocalNbGeoEnts();
  const CFuint nbStates = 3;  // triangles hardcoding

  // nodes coordinates, area, state values and IDs
  CFreal x[nbStates], y[nbStates], z[nbStates], A;
  RealVector states(nbStates*nbEqs);
  CFuint statesID[nbStates];
  bool statesIsParUpdatable[nbStates];

  // block accumulator to jacobian matrix
  BlockAccumulator *acc = getMethodData().getLinearSystemSolver()[0]->
    createBlockAccumulator(nbStates,nbStates,nbEqs);

  // element contributions averaged, unperturbed and perturbed and its derivative
  RealVector JnAvg(nbEqs);
  RealVector JpAvg(nbEqs);
  RealVector dJdU(nbEqs);

  // variables to evaluate the boundary contributions
  // number of variables: dimensions, time and PhysicalModel dimensions
  RealVector variables(dim+1+nbEqs);
  variables[dim] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

  // integral result, unperturbed and perturbed and its derivative
  RealVector _cNormal(nbEqs);
  RealVector _cPertbd(nbEqs);
  RealVector _dInt(nbEqs);

  // get system matrix and strategy
  SafePtr<LSSMatrix> jacobMatrix = getMethodData().getLinearSystemSolver()[0]->
    getMatrix();
  NumericalJacobian& numericalJacob = getMethodData().
    getNumericalJacobian();

  // provide access to the GeometricEntities
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;


  // for all the Triangles
  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) {

    CFLogDebugMax("Assembling BC Triangle " << iGeoEnt << "\n");

    // build the GeometricEntity
    geoData.idx = iGeoEnt;
    GeometricEntity& cell = *geoBuilder->buildGE();

    // get nodal coordinates and cell states (isoparametric hardcoding)
    getNodes(cell,nbStates, x,y,z);
    getStates(cell,nbStates, states,statesID,statesIsParUpdatable);

    // calculate area
    calcArea(x,y,z, A);

    // reset block accumulator
    acc->reset();


    if (!_useAverage) {


      // use nodal boundary condition contribution

      // for all the states in this cell
      for (CFuint iState=0; iState<nbStates; ++iState) {

        // only calculate if this state is updatable in current computing node
        if ( statesIsParUpdatable[iState] ) {

          // get necessary variables for evaluation of VectorialFunction,
          // variables space, time (omitted) and states
          variables[0] = x[iState];
          variables[1] = y[iState];
          variables[2] = z[iState];
          for(CFuint iEq=0; iEq<nbEqs; iEq++) {
            variables[dim+1+iEq] = states[iState*nbEqs+iEq];
          }

          // integrate non-perturbed bc contribution
          // (hardcodes the weight of the shape function for triangles, 1./3.)
          _vFunction.evaluate(variables,_cNormal);
          _cNormal *= A/3.;

          // set the index of the block corresponding to the current state in the jacobian matrix
          acc->setRowColIndex(iState, statesID[iState]);

          // perturb one component at a time
          for (CFuint iEq=0; iEq<nbEqs; ++iEq) {

            // add rhs contribution, distributed per node
            rhs(statesID[iState], iEq, nbEqs) += _cNormal[iEq];

            // perturb state and integrate bc contribution
            // (hardcodes the weight of the shape function for triangles, 1./3.)
            numericalJacob.perturb(iEq,variables[dim+1+iEq]);
            _vFunction.evaluate(variables,_cPertbd);
            numericalJacob.restore(variables[dim+1+iEq]);
            _cPertbd *= A/3.;

            // jacobian contribution (accumulator)
            // add vector to rows [iState-1 to iState-nbEqs], column [iState-iEq]
            numericalJacob.computeDerivative(_cPertbd,_cNormal,_dInt);
            acc->addValues(iState, iState, iEq, &_dInt[0] );

          } // iEq

        } // statesIsParUpdatable[iState]
      }  // for iState


    } else {


      // use element-averaged boundary condition contribution

      // loop to calculate the average contributions (normal and perturbed)
      JnAvg=0.;
      JpAvg=0.;
      for (CFuint iState=0; iState<nbStates; ++iState) {

        // get necessary variables for evaluation of VectorialFunction
        variables[0] = x[iState];
        variables[1] = y[iState];
        variables[2] = z[iState];
        for(CFuint iEq=0; iEq<nbEqs; iEq++) {
          variables[dim+1+iEq] = states[iState*nbEqs+iEq];
        }

        // integrate nodal contributions (normal and perturbed) and update averages
        // (hardcodes the weight of the shape function for triangles, 1./3.)
        _vFunction.evaluate(variables,_cNormal);
        _cNormal *= A/3.;
        for (CFuint iEq=0; iEq<nbEqs; iEq++) {
          numericalJacob.perturb(iEq,variables[dim+1+iEq]);
          _vFunction.evaluate(variables,_cPertbd);
          numericalJacob.restore(variables[dim+1+iEq]);
          _cPertbd *= A/3.;

          JnAvg[iEq] = ( JnAvg[iEq]*iState+_cNormal[iEq] )/( iState+1. );
          JpAvg[iEq] = ( JpAvg[iEq]*iState+_cPertbd[iEq] )/( iState+1. );
        };
      }

      // calculate averaged jacobian contribution, distributed per node
      numericalJacob.computeDerivative(JpAvg,JnAvg,dJdU);

      // for all the states in this cell
      for (CFuint iState=0; iState<nbStates; ++iState) {

        // only calculate if this state is updatable in current computing node
        if ( statesIsParUpdatable[iState] ) {

          // add to rhs: averaged contribution, distributed per node
          // add to accumulator: vector to rows [iState-1 to iState-nbEqs], column [iState-iEq]
          for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
            rhs(statesID[iState], iEq, nbEqs) += JnAvg[iEq];
            acc->addValues(iState, iState, iEq, &dJdU[0] );
          }

          // set the index of the block corresponding to the current state in the jacobian matrix
          acc->setRowColIndex(iState, statesID[iState]);

        } // statesIsParUpdatable[iState]
      }  // for iState


    } // _useAverage?


    // add the values in the jacoban matrix
    jacobMatrix->addValues( (* acc) );

    //release the GeometricEntity
    geoBuilder->releaseGE();

  } // for iGeoEnt

/*
MathTools::ConstantFunctor<COOLFluiD::DIM_1D> _functor;
_functor.setValue(1000./3.);
getMethodData().getVolumeIntegrator()->
integrateConstantFunctorOnGeoEnt<ConstantFunctor<DIM_1D> >( &cell ,_functor,_dInt);
CF_DEBUG_OBJ(state);
*/

}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicitP1Analytical::configure ( Config::ConfigArgs& args )
{
  FiniteElementMethodCom::configure(args);

  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);
  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NeumannBCImplicitP1Analytical::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

