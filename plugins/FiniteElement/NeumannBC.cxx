#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/Node.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/NormalsCalculator.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/NeumannBC.hh"
#include "FiniteElement/NeumannEntity.hh"

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

MethodCommandProvider<NeumannBC, FiniteElementMethodData, FiniteElementModule> NeumannBCProvider("NeumannBC");

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

NeumannBC::NeumannBC(const std::string& name) :
  FiniteElementMethodCom(name),
  _localElemData(CFNULL),
  socket_rhs("rhs"),
  socket_isUpdated("isUpdated"),
  socket_updateCoeff("updateCoeff"),
  socket_states("states"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  _integResult(0)
{
   addConfigOptionsTo(this);
  _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);
}

//////////////////////////////////////////////////////////////////////////////

NeumannBC::~NeumannBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::setup()
{
  // first call parent method
  FiniteElementMethodCom::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // settin correct size for integration result
  _integResult.resize(nbEqs);

  // add specific configuration here
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

void NeumannBC::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::executeOnTrs()
{
  CFAUTOTRACE;

  if(getMethodData().isDirichletBCApplied())
  {
     CFLog(WARN," #!!!!!!!!!!!!!!!\n You are applying a NeumannBC after a DirichletBC...\n expect strange results!!!\n #!!!!!!!!!!!!!!!\n");
  }

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "NeumannBC::execute() called for TRS: " << trs->getName() << "\n");

  // Get the NeumannTerm
  // and setup the VectorialFunction in the NeumannTerm

  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbGeos = trs->getLocalNbGeoEnts();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;
  _localElemData->trs = trs;
  _neumannEntity->setVectorialFunction(&_vFunction);

  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) {
    // build the GeometricEntity
    geoData.idx = iGeoEnt;
    GeometricEntity *const face = geoBuilder->buildGE();
    _localElemData->cell  = face;

    vector<State*>* faceStates = face->getStates();

    CFLogDebugMax("Face " << iGeoEnt << "\n");

    // Integration of Boundary condition functor
    const CFuint nbStatesInFace = face->getNbNodesSolutionShapeFunction();
    for (CFuint iState = 0; iState < nbStatesInFace; ++iState) {

      State *const currState = (*faceStates)[iState];
      _localElemData->iState = iState;

      // Integrate the Equation Term for the stiffness matrix
      getMethodData().getFEMVolumeIntegrator()->
          integrateFaceFEMEntityOnGeoEnt<NeumannEntity>(*_neumannEntity, _integResult);

        const CFuint stateID = currState->getLocalID();
        // rhs is set to
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(stateID, iEq, nbEqs) += _integResult[iEq];
        }

        isUpdated[stateID] = true; // flagging is important!!!!!
    }

    //release the GeometricEntity
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

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
NeumannBC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_states);
  result.push_back(&socket_bStatesNeighbors);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
