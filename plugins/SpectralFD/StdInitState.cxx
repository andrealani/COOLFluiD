#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "SpectralFD/StdInitState.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdInitState, SpectralFDMethodData, SpectralFDModule>
StdInitStateProvider("StdInitState");

//////////////////////////////////////////////////////////////////////////////

void StdInitState::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::string >("InputVar","Input variables.");
}

//////////////////////////////////////////////////////////////////////////////

StdInitState::StdInitState(const std::string& name) :
  SpectralFDMethodCom(name),
  m_varSet(CFNULL),
  m_inputToUpdateVar(),
  m_initPntsStates(),
  m_inputState(),
  m_nbrEqs(),
  m_dim(),
  m_initPntCoords()
{
  addConfigOptionsTo(this);
  m_functions = vector<std::string>();
  setParameter("Def",&m_functions);

  m_vars = vector<std::string>();
  setParameter("Vars",&m_vars);

  m_inputVarStr = "";
  setParameter("InputVar",&m_inputVarStr);
}

//////////////////////////////////////////////////////////////////////////////

StdInitState::~StdInitState()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    StdInitState::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdInitState::setup()
{

  CFAUTOTRACE;
  SpectralFDMethodCom::setup();

  // get number of physical variables
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim ();

  // resize variables
  m_initPntCoords.resize(m_dim);

  // set maxNbStatesInCell for variable transformer
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_inputToUpdateVar->setup(maxNbStatesInCell);

  m_initPntsStates.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState)
  {
    m_initPntsStates[iState] = new State();
  }

  // create state for input variables
  m_inputState = new State();

  // get updateVarSet
  m_varSet = getMethodData().getUpdateVar();
}

//////////////////////////////////////////////////////////////////////////////

void StdInitState::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iState = 0; iState < m_initPntsStates.size(); ++iState)
  {
    deletePtr(m_initPntsStates[iState]);
  }
  m_initPntsStates.resize(0);

  deletePtr(m_inputState);

  SpectralFDMethodCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void StdInitState::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  SpectralFDMethodCom::configure(args);

  // get the physical model that we are dealing with
  // to pass it to the variable transformer
  std::string namespc = getMethodData().getNamespace();
  SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(namespc);
  SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  // get the name of the update variable set
  std::string _updateVarStr = getMethodData().getUpdateVarStr();

  // create the transformer from input to update variables
  if (m_inputVarStr.empty()) {
    m_inputVarStr = _updateVarStr;
  }

  std::string provider =
    VarSetTransformer::getProviderName(physModel->getNameImplementor(), m_inputVarStr, _updateVarStr);

  m_inputToUpdateVar =
    Environment::Factory<VarSetTransformer>::getInstance()
    .getProvider(provider)->create(physModel->getImplementor());
  cf_assert(m_inputToUpdateVar.isNotNull());

  // parsing the functions that the user inputed
  m_vFunction.setFunctions(m_functions);
  m_vFunction.setVariables(m_vars);
  try {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdInitState::executeOnTrs()
{
  CFAUTOTRACE;

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");
  CFLogDebugMin("StdInitState::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  // prepares to loop over cells by getting the GeometricEntityPool
  SafePtr< GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

  // loop over elements
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get the number of elements
    const CFuint nbrElems = (*elemType)[iElemType].getNbElems();

    // get start index of this element type in global element list
    CFuint cellIdx = (*elemType)[iElemType].getStartIdx();

    // local coordinates of initialization points and matrix used for initialization
    SafePtr< vector< RealVector > > locInitPntCoords = sdLocalData[iElemType]->getInitPntsCoords();
    SafePtr< RealMatrix >           initTransfMatrix = sdLocalData[iElemType]->getInitTransfMatrix();

    // get number of states
    const CFuint nbrStates = sdLocalData[iElemType]->getNbrOfSolPnts();

    // loop over elements
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
    {
      // build the GeometricEntity
      geoData.idx = cellIdx;
      GeometricEntity *const cell = geoBuilder->buildGE();

      // get the states
      vector<State*>* solPntStates = cell->getStates();
      cf_assert(solPntStates->size() == nbrStates);

      // loop over initialization points
      for (CFuint iPnt = 0; iPnt < nbrStates; ++iPnt)
      {
        // compute initialization node global coordinate
        m_initPntCoords = cell->computeCoordFromMappedCoord((*locInitPntCoords)[iPnt]);

        // evaluate the function at the state coordinate
        m_vFunction.evaluate(m_initPntCoords,*m_inputState);
        *m_initPntsStates[iPnt] = *m_inputToUpdateVar->transform(m_inputState);
      }

      // transform initialization points solutions to solution point solutions
      for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
      {
        // variable for a dimensional state
        State dimState(RealVector(0.0,m_nbrEqs));

        for (CFuint iPnt = 0; iPnt < nbrStates; ++iPnt)
        {
          dimState += (*initTransfMatrix)(iSol,iPnt)*(*m_initPntsStates[iPnt]);
        }

        // adimensionalize the value if needed and store
        m_varSet->setAdimensionalValues(dimState, *(*solPntStates)[iSol]);
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
