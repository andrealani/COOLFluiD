#include "Common/PE.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/Node.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/NormalsCalculator.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/CoupledNeumannBC.hh"
#include "FiniteElement/CoupledNeumannEntity.hh"

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

MethodCommandProvider<CoupledNeumannBC, FiniteElementMethodData, FiniteElementModule> CoupledNeumannBCProvider("CoupledNeumannBC");

//////////////////////////////////////////////////////////////////////////////

void CoupledNeumannBC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Interface","Name of the Interface.");
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");

   options.addConfigOption< bool >("AlternateBC","Alternate between different BCs");
/*   options.addConfigOption< CFuint >("AlternateRate","Alternate Rate between different BCs");*/
   options.addConfigOption< bool >("AlternateStart","Start with this BC when alternating between different BCs");
}

//////////////////////////////////////////////////////////////////////////////

CoupledNeumannBC::CoupledNeumannBC(const std::string& name) :
  FiniteElementMethodCom(name),
  _sockets(),
  socket_rhs("rhs"),
  socket_isUpdated("isUpdated"),
  socket_bStatesNeighbors("bStatesNeighbors"),
  _integResult(0)
{
   addConfigOptionsTo(this);
  _interfaceName = "";
   setParameter("Interface",&_interfaceName);

  _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);

  _alternateBC = false;
  setParameter("AlternateBC",&_alternateBC);

//   _alternateRate = 2;
//   setParameter("AlternateRate",&_alternateRate);

  _alternateStart = true;
  setParameter("AlternateStart",&_alternateStart);

  _isRobinBC = false;
}

//////////////////////////////////////////////////////////////////////////////

CoupledNeumannBC::~CoupledNeumannBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void CoupledNeumannBC::setup()
{
  CFAUTOTRACE;

  // first call parent method
  FiniteElementMethodCom::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // settin correct size for integration result
  _integResult.resize(nbEqs);

  _currentAlternateRun = _alternateStart;
  if(!_alternateBC) _currentAlternateRun = true;


  // add specific configuration here
  _localElemData = &(getMethodData().getLocalElementData());

  std::string neumannEntityStr = "Galerkin" + getMethodData().getDiffusiveVar()->getName();

  Common::SafePtr<
  BaseMethodStrategyProvider<FiniteElementMethodData,CoupledNeumannEntity> > provNeumannEnt =
    Environment::Factory<CoupledNeumannEntity>::getInstance().getProvider(neumannEntityStr);
  cf_assert(provNeumannEnt.isNotNull());
  _coupledNeumannEntity = provNeumannEnt->create(neumannEntityStr,SharedPtr<FiniteElementMethodData>(&getMethodData()));

  configureNested ( _coupledNeumannEntity.getPtr(), m_stored_args );

  _coupledNeumannEntity->setup();

}

//////////////////////////////////////////////////////////////////////////////

void CoupledNeumannBC::unsetup()
{
  CFAUTOTRACE;

  // then call parent method
  FiniteElementMethodCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void CoupledNeumannBC::executeOnTrs()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "CoupledNeumannBC::execute() called for TRS: " << trs->getName() << "\n");

  if(_currentAlternateRun)
  {
    std::cout << "Running NeumannBC for TRS: " << trs->getName() << "(Interface: " << _interfaceName <<")" << std::endl;

/*    if(getMethodData().isDirichletBCApplied())
    {
      CFLog(WARN," #!!!!!!!!!!!!!!!\n You are applying a NeumannBC after a DirichletBC...\n expect strange results!!!\n #!!!!!!!!!!!!!!!\n");
    }*/

    // Get all the datahandles
    DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
    DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();

    // Get the datahandles of accepted data
    const std::string trsName = trs->getName();
    SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
    std::vector<CFuint>::iterator itd;

    const std::string currentSubSystem = SubSystemStatusStack::getActive()->getSubSystemName();
    const std::string nameSpace = getMethodData().getNamespace();

    const std::string baseName =
        "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Gauss_";
    const std::string socketDataName = baseName + "DATA";
    const std::string socketAcceptName = baseName + "ISACCEPTED";

    DataHandle<RealVector>
      interfaceData = _sockets.getSocketSink<RealVector>(socketDataName)->getDataHandle();
    DataHandle<CFreal>
      isAccepted = _sockets.getSocketSink<CFreal>(socketAcceptName)->getDataHandle();

    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    const CFuint nbGeos = trs->getLocalNbGeoEnts();

    Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
      geoBuilder = getMethodData().getStdTrsGeoBuilder();

    StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
    geoData.trs = trs;
    _localElemData->trs = trs;

    _coupledNeumannEntity->setDataHandle(interfaceData);
    _coupledNeumannEntity->setIsAcceptedDataHandle(isAccepted);
    _coupledNeumannEntity->setVectorialFunction(&_vFunction);

    CFuint idxData(0);
    CFuint idxAccepted(0);
    for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) {

      // build the GeometricEntity
      geoData.idx = iGeoEnt;

      GeometricEntity *const face = geoBuilder->buildGE();
      vector<State*> *const faceStates = face->getStates();
      _localElemData->cell  = face;

      CFLogDebugMax("Face " << iGeoEnt << "\n");

      // Integration of Boundary condition functor
      const CFuint nbStatesInFace = face->getNbNodesSolutionShapeFunction();
      for (CFuint iState = 0; iState < nbStatesInFace; ++iState) {
        State *const currState = (*faceStates)[iState];
        _localElemData->iState = iState;

        // Integrate the Equation Term for the stiffness matrix
        _coupledNeumannEntity->setDataHandleIndex(idxData);
        _coupledNeumannEntity->setIsAcceptedDataHandleIndex(idxAccepted);

        getMethodData().getFEMVolumeIntegrator()->
          integrateFaceFEMEntityOnGeoEnt<CoupledNeumannEntity>(*_coupledNeumannEntity, _integResult);

        const CFuint stateID = currState->getLocalID();
        // rhs is set to
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(stateID, iEq, nbEqs) += _integResult[iEq];
        }

        isUpdated[stateID] = true; // flagging is important!!!!!
      }

      //Compute the Nb of Accepted Points
      const CFuint nbQuadPoints =
        getMethodData().getVolumeIntegrator()->
        getSolutionIntegrator(face)->getIntegratorPattern().totalNbPts();

      //increase the counters
      DataHandle<CFreal>
        isAccepted = _sockets.getSocketSink<CFreal>(socketAcceptName)->getDataHandle();

      // Update the counters for the nbQuadPoints
      CFuint nbAcceptedPointsInFace = 0;
      for(CFuint iQuad = 0; iQuad < nbQuadPoints ;iQuad++) {
        if(isAccepted[idxAccepted + iQuad] >= 0.) nbAcceptedPointsInFace++;
      }

      idxData += nbAcceptedPointsInFace;
      idxAccepted += nbQuadPoints;

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }

  //update currentAlternateRun
  const CFuint currentIteration = SubSystemStatusStack::getActive()->getNbIter();
  if((currentIteration > 0) &&(_alternateBC)) _currentAlternateRun = !_currentAlternateRun;

}

//////////////////////////////////////////////////////////////////////////////

void CoupledNeumannBC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  FiniteElementMethodCom::configure(args);

  // Function used for the NeumannBC
  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);
  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

  // Sockets to be used
  const std::string nameSpace = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(nameSpace);
  Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

  const std::string currentSubSystem = subsystemStatus->getSubSystemName();
  const std::vector<std::string>& trsNames = getTrsNames();

  for(CFuint iTRS = 0; iTRS < trsNames.size(); iTRS++)
  {
    const std::string trsName = trsNames[iTRS];

    const std::string baseName = "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Gauss_";
    const std::string socketDataName = baseName + "DATA";
    const std::string socketAcceptName = baseName + "ISACCEPTED";

    _sockets.createSocketSink<CFreal>(socketAcceptName);
    _sockets.createSocketSink<RealVector>(socketDataName);

  }

  m_stored_args = args;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CoupledNeumannBC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();

  result.push_back(&socket_rhs);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_bStatesNeighbors);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
