#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/VolumeIntegrator.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/NeumannEntity.hh"
#include "FiniteElement/NeumannBCImplicit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NeumannBCImplicit, FiniteElementMethodData, FiniteElementModule> NeumannBCImplicitProvider("NeumannBCImplicit");

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicit::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

NeumannBCImplicit::NeumannBCImplicit(const std::string& name) :
  FiniteElementMethodCom(name),
  _localElemData(CFNULL),
  socket_rhs("rhs")
{
   addConfigOptionsTo(this);
  _functions = std::vector<std::string>();
   setParameter("Def",&_functions);



  _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);


}

//////////////////////////////////////////////////////////////////////////////

NeumannBCImplicit::~NeumannBCImplicit()
{
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicit::setup()
{
  // first call parent method
  FiniteElementMethodCom::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // set correct size for integration result
  _intNormal.resize(nbEqs);
  _intPertbd.resize(nbEqs);
  _dInt.resize(nbEqs);


  _localElemData = &(getMethodData().getLocalElementData());

  std::string neumannEntityStr = "Galerkin" + getMethodData().getDiffusiveVar()->getName();
  Common::SafePtr<
  BaseMethodStrategyProvider<FiniteElementMethodData,NeumannEntity> > provNeumannEnt =
    Environment::Factory<NeumannEntity>::getInstance().getProvider(neumannEntityStr);
  cf_assert(provNeumannEnt.isNotNull());
  _neumannEntity = provNeumannEnt->create(neumannEntityStr,SharedPtr<FiniteElementMethodData>(&getMethodData()));
  configureNested ( _neumannEntity.getPtr(), m_stored_args );

  _neumannEntity->setup();

}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicit::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicit::executeOnTrs()
{
  CFAUTOTRACE;

  if(getMethodData().isDirichletBCApplied())
  {
    CFLog(WARN," #!!!!!!!!!!!!!!!\n You are applying a NeumannBC after a DirichletBC...\n expect strange results!!!\n #!!!!!!!!!!!!!!!\n");
  }

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "NeumannBCImplicit::execute() called for TRS: " <<
                 trs->getName() << "\n");

  // accumulators and pointer
  std::vector<bool> vHaveThisSize;
  std::vector<BlockAccumulator*> vFaceAcc;
  BlockAccumulator* faceAcc;

  // jacobian properties
  SafePtr<LSSMatrix> jacobMatrix =
    getMethodData().getLinearSystemSolver()[0]->getMatrix();
  NumericalJacobian& numericalJacob = getMethodData().getNumericalJacobian();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbGeos = trs->getLocalNbGeoEnts();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;
  _localElemData->trs = trs;
  _neumannEntity->setVectorialFunction(&_vFunction);

  // for all the faces
  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) {

    // build the GeometricEntity
    geoData.idx = iGeoEnt;

    GeometricEntity *const face = geoBuilder->buildGE();
    std::vector<State*>* faceStates = face->getStates();
    CFuint nbStatesInFace = face->getNbNodesSolutionShapeFunction();
    _localElemData->cell  = face;
    // faceStates.size() //?

    CFLogDebugMin("Face " << iGeoEnt << "\n");

    // get appropriate accumulator, creating new of required dimension if needed
    if (vHaveThisSize.size()+1 < nbStatesInFace) { // vectors are zero-base indexed
      vFaceAcc.resize(nbStatesInFace);
      vHaveThisSize.resize(nbStatesInFace,false);
    }
    if (!vHaveThisSize[nbStatesInFace-1]) {
      vFaceAcc[nbStatesInFace-1] = getMethodData().getLinearSystemSolver()[0]->createBlockAccumulator(nbStatesInFace,nbStatesInFace,nbEqs);
      vHaveThisSize[nbStatesInFace-1] = true;
    }
    faceAcc = vFaceAcc[nbStatesInFace-1];
    faceAcc->reset();


    // for all the states in this face
    for (CFuint iState = 0; iState < nbStatesInFace; ++iState) {

      // get state and its ID
      State * currState = (*faceStates)[iState];
      CFuint currStateID = currState->getLocalID();
      _localElemData->iState = iState;

      // only calculate if this state is updatable in current computing node
      if ( currState->isParUpdatable() ) {

        // Integrate non perturbed BC contribution (f_iState)
        getMethodData().getFEMVolumeIntegrator()->
          integrateFaceFEMEntityOnGeoEnt<NeumannEntity>(*_neumannEntity, _intNormal);

        // rhs is set to
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(currStateID, iEq, nbEqs) += _intNormal[iEq];
        }

        // set the index of the block corresponding to the current state in the jacobian matrix
        faceAcc->setRowColIndex(iState,currStateID);

        // loop over the variables in the state vector to perturb one component at a time
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {

          CFLogDebugMax( "NeumannBCImplicit: Perturbing iState = " << iState << " at iEq = " << iEq << "\n");

                // calculate perturbed BC contribution, f_(iState,iEq)*
          // CFreal AntesDaPert = (*currState)[0];
          numericalJacob.perturb(iEq, (*currState)[iEq]);
          getMethodData().getFEMVolumeIntegrator()->
            integrateFaceFEMEntityOnGeoEnt<NeumannEntity>(*_neumannEntity, _intPertbd);
          numericalJacob.restore((*currState)[iEq]);


          // jacobian contribution -df_iState/dU_iState (vector 1 to nbEqs)
          numericalJacob.computeDerivative(_intPertbd,_intNormal,_dInt);

          // add vector to row [iState-1 to iState-nbEqs], column [iState-iEq]
          faceAcc->addValues(iState, iState, iEq, &_dInt[0] );

        } // iEq

      } // currState->isParUpdatable()
    } // for iState

    // add the values in the jacoban matrix
    jacobMatrix->addValues( (* faceAcc) );

    //release the GeometricEntity
    geoBuilder->releaseGE();
  } // iGeoEnt

}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCImplicit::configure ( Config::ConfigArgs& args )
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
  m_stored_args = args;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NeumannBCImplicit::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
