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
#include "FiniteElement/CoupledNeumannImplicitBC.hh"

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

MethodCommandProvider<CoupledNeumannImplicitBC, FiniteElementMethodData, FiniteElementModule>
CoupledNeumannImplicitBCProvider("CoupledNeumannImplicitBC");

//////////////////////////////////////////////////////////////////////////////

void CoupledNeumannImplicitBC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Interface","Name of the Interface.");
   options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");

   options.addConfigOption< bool >("AlternateBC","Alternate between different BCs");
/*   options.addConfigOption< CFuint >("AlternateRate","Alternate Rate between different BCs");*/
   options.addConfigOption< bool >("AlternateStart","Start with this BC when alternating between different BCs");
}

//////////////////////////////////////////////////////////////////////////////

CoupledNeumannImplicitBC::CoupledNeumannImplicitBC(const std::string& name) :
  FiniteElementMethodCom(name),
  _sockets(),
  socket_rhs("rhs")
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

CoupledNeumannImplicitBC::~CoupledNeumannImplicitBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void CoupledNeumannImplicitBC::setup()
{
  CFAUTOTRACE;

  // first call parent method
  FiniteElementMethodCom::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // set correct size for integration result
  _intNormal.resize(nbEqs);
  _intPertbd.resize(nbEqs);
  _dInt.resize(nbEqs);

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

void CoupledNeumannImplicitBC::unsetup()
{
  CFAUTOTRACE;

  // then call parent method
  FiniteElementMethodCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void CoupledNeumannImplicitBC::executeOnTrs()
{
  CFAUTOTRACE;

  Common::SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "CoupledNeumannImplicitBC::execute() called for TRS: " << trs->getName() << "\n");

  if(_currentAlternateRun)
  {
    // accumulators and pointer
    std::vector<bool> vHaveThisSize;
    std::vector<BlockAccumulator*> vFaceAcc;
    BlockAccumulator* faceAcc;

    // jacobian properties
    SafePtr<LSSMatrix> jacobMatrix =
      getMethodData().getLinearSystemSolver()[0]->getMatrix();
    NumericalJacobian& numericalJacob = getMethodData().getNumericalJacobian();


//    std::cout << "Running NeumannBC" << std::endl;

/*    if(getMethodData().isDirichletBCApplied())
    {
      CFLog(WARN," #!!!!!!!!!!!!!!!\n You are applying a NeumannBC after a DirichletBC...\n expect strange results!!!\n #!!!!!!!!!!!!!!!\n");
    }*/

    // Get the datahandles
    DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

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
    for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt)
    {

      // build the GeometricEntity
      geoData.idx = iGeoEnt;

      GeometricEntity *const face = geoBuilder->buildGE();
      vector<State*> *const faceStates = face->getStates();
      _localElemData->cell  = face;

      CFLogDebugMax("Face " << iGeoEnt << "\n");

      // Integration of Boundary condition functor
      const CFuint nbStatesInFace = face->getNbNodesSolutionShapeFunction();

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

      for (CFuint iState = 0; iState < nbStatesInFace; ++iState)
      {
        State *const currState = (*faceStates)[iState];
        CFuint currStateID = currState->getLocalID();

        // only calculate if this state is updatable in current computing node
        if ( currState->isParUpdatable() ) {

          _localElemData->iState = iState;

          // Integrate the Equation Term for the stiffness matrix
          _coupledNeumannEntity->setDataHandleIndex(idxData);
          _coupledNeumannEntity->setIsAcceptedDataHandleIndex(idxAccepted);
          getMethodData().getFEMVolumeIntegrator()->
            integrateFaceFEMEntityOnGeoEnt<CoupledNeumannEntity>(*_coupledNeumannEntity, _intNormal);

          // rhs is set to
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) { rhs(currStateID, iEq, nbEqs) += _intNormal[iEq];  }

          // set the index of the block corresponding to the current state in the jacobian matrix
          faceAcc->setRowColIndex(iState,currStateID);

          // loop over the variables in the state vector to perturb one component at a time
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {

            CFLogDebugMax( "NeumannBCImplicit: Perturbing iState = " << iState << " at iEq = " << iEq << "\n");

                  // calculate perturbed BC contribution, f_(iState,iEq)*
            // CFreal AntesDaPert = (*currState)[0];
            numericalJacob.perturb(iEq, (*currState)[iEq]);

            _localElemData->iState = iState;

            // Integrate the Equation Term for the stiffness matrix
            _coupledNeumannEntity->setDataHandleIndex(idxData);
            _coupledNeumannEntity->setIsAcceptedDataHandleIndex(idxAccepted);
            getMethodData().getFEMVolumeIntegrator()->
              integrateFaceFEMEntityOnGeoEnt<CoupledNeumannEntity>(*_coupledNeumannEntity, _intPertbd);

            numericalJacob.restore((*currState)[iEq]);

            // jacobian contribution -df_iState/dU_iState (vector 1 to nbEqs)
            numericalJacob.computeDerivative(_intPertbd,_intNormal,_dInt);

            // add vector to row [iState-1 to iState-nbEqs], column [iState-iEq]
            faceAcc->addValues(iState, iState, iEq, &_dInt[0] );

          } // iEq
        } // currState->isParUpdatable()
      } //iState

      // add the values in the jacoban matrix
      jacobMatrix->addValues( (* faceAcc) );

      //release the GeometricEntity
      geoBuilder->releaseGE();

      //Compute the Nb of Accepted Points
      const CFuint nbQuadPoints =
        getMethodData().getVolumeIntegrator()->
        getSolutionIntegrator(face)->getIntegratorPattern().totalNbPts();

      //increase the counters
      DataHandle<CFreal>
        isAccepted = _sockets.getSocketSink<CFreal>(socketAcceptName)->getDataHandle();

      CFuint nbAcceptedPointsInFace = 0;
      for(CFuint iQuad = 0; iQuad < nbQuadPoints ;iQuad++) {
        if(isAccepted[idxAccepted + iQuad] >= 0.) nbAcceptedPointsInFace++;
      }

      idxData += nbAcceptedPointsInFace;
      idxAccepted += nbQuadPoints;
    }
  }

  //update currentAlternateRun
  const CFuint currentIteration = SubSystemStatusStack::getActive()->getNbIter();
  if((currentIteration > 0) &&(_alternateBC)) _currentAlternateRun = !_currentAlternateRun;

}

//////////////////////////////////////////////////////////////////////////////

void CoupledNeumannImplicitBC::configure ( Config::ConfigArgs& args )
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
  const std::vector<std::string> trsNames = getTrsNames();

  for(CFuint iTRS = 0; iTRS < trsNames.size(); iTRS++)
  {
    const std::string trsName = trsNames[iTRS];

    const std::string baseName =
        "COUPLING_" + _interfaceName + "_" + trsName + "_" + nameSpace + "_" + currentSubSystem + "_Gauss_";
    const std::string socketDataName = baseName + "DATA";
    const std::string socketAcceptName = baseName + "ISACCEPTED";

    _sockets.createSocketSink<CFreal>(socketAcceptName);
    _sockets.createSocketSink<RealVector>(socketDataName);
  }

  m_stored_args = args;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CoupledNeumannImplicitBC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();

  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
