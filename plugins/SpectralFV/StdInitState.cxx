#include "SpectralFV/SpectralFV.hh"


#include "SpectralFV/StdInitState.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdInitState, SpectralFVMethodData, SpectralFVModule>
StdIQunitStateProvider("StdInitState");

//////////////////////////////////////////////////////////////////////////////

void StdInitState::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::string >("InputVar","Input variables.");
}

//////////////////////////////////////////////////////////////////////////////

StdInitState::StdInitState(const std::string& name) :
  SpectralFVMethodCom(name),
  m_varSet(CFNULL),
  m_inputToUpdateVar(),
  m_initPntsStates(),
  m_inputState(),
  m_nbrEqs()
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

void StdInitState::setup()
{

  CFAUTOTRACE;
  SpectralFVMethodCom::setup();

  // get number of physical variables
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

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

  SpectralFVMethodCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void StdInitState::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  SpectralFVMethodCom::configure(args);

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

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

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

  // only element types with a P1 geometry at most can be treated in this loop
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get the number of elements
    const CFuint nbrElems = (*elemType)[iElemType].getNbElems();

    // get start index of this element type in global element list
    CFuint cellIdx = (*elemType)[iElemType].getStartIdx();

    // check if geometric shape function is not higher order than 1
    if ((*elemType)[iElemType].getGeoOrder() > 1) throw Common::ShouldNotBeHereException (FromHere(),"Only elements with a P0 or a P1 mapping to a reference element can have a linear transformation to this reference element");

    // variables for cv volume fractions, local cv - local node connectivity, local node coordinates
    SafePtr< vector< RealVector > > locInitPntCoords = svLocalData[iElemType]->getInitPntsCoords();
    SafePtr< RealMatrix >           initTransfMatrix = svLocalData[iElemType]->getInitTransfMatrix();

    // get number of control volumes
    const CFuint nbrCVs = svLocalData[iElemType]->getNbrOfCVs();

    // get dimensionality
    const CFDim dimensionality = svLocalData[iElemType]->getDimensionality();

    // vector for initialization node coordinates
    RealVector initPntCoords(dimensionality);

    // loop over elements
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
    {
      // build the GeometricEntity
      geoData.idx = cellIdx;
      GeometricEntity *const cell = geoBuilder->buildGE();

      // get the SV nodes
      vector<Node*>* svNodes = cell->getNodes();

      // get the states
      vector<State*>* svStates = cell->getStates();
      cf_assert(svStates->size() == nbrCVs);

      // loop over initialization points
      for (CFuint iNode = 0; iNode < nbrCVs; ++iNode)
      {
        // compute CV node coordinates
        initPntCoords = (1.0 - (*locInitPntCoords)[iNode].sum())*(*(*svNodes)[0]);
        for (CFuint iCoor = 0; iCoor < (CFuint) dimensionality; ++iCoor)
        {
          initPntCoords += (*locInitPntCoords)[iNode][iCoor]*(*(*svNodes)[iCoor+1]);
        }

        // evaluate the function at the state coordinate
        m_vFunction.evaluate(initPntCoords,*m_inputState);
        *m_initPntsStates[iNode] = *m_inputToUpdateVar->transform(m_inputState);
      }

      // transform points solutions to CV-averaged solutions
      for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
      {
        // variable for a dimensional state
        State dimState(RealVector(0.0,m_nbrEqs));

        for (CFuint iNode = 0; iNode < nbrCVs; ++iNode)
        {
          dimState += (*initTransfMatrix)(iCV,iNode)*(*m_initPntsStates[iNode]);
        }

        // adimensionalize the value if needed and store
        m_varSet->setAdimensionalValues(dimState, *(*svStates)[iCV]);
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD
