#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

#include "LinEuler/LinearizedEuler.hh"
#include "LinEuler/LinEuler3DCons.hh"
#include "LinEuler/LinEuler3DVarSet.hh"
#include "LinEuler/LinEulerTerm.hh"
#include "LinEuler/LinEuler3DLinearCons.hh"

#include "SpectralFDLinEuler/LEEInitState3D.hh"

#include "Framework/MeshData.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SubSystemStatus.hh"

////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearizedEuler;

////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

////////////////////////////////////////////////////////////////////////////
//

MethodCommandProvider<LEEInitState3D, SpectralFDMethodData, SpectralFDModule>
LEEInitState3DProvider("LEEInitState3D");

////////////////////////////////////////////////////////////////////////////
//

std::vector<Common::SafePtr<BaseDataSocketSink> >
LEEInitState3D::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_meanflow);

  return result;
}

////////////////////////////////////////////////////////////////////////////
//

void LEEInitState3D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::string >("InputVar","Input variables.");
}

////////////////////////////////////////////////////////////////////////////

LEEInitState3D::LEEInitState3D(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_meanflow("meanflow"),
  m_varSet(CFNULL),
  m_inputToUpdateVar(),
  m_initPntsStates(),
  m_inputState(),
  m_nbrEqs(),
  m_dim(),
  m_initPntCoords(),
  _varSet()
{
  addConfigOptionsTo(this);
  m_functions = vector<std::string>();
  setParameter("Def",&m_functions);

  m_vars = vector<std::string>();
  setParameter("Vars",&m_vars);

  m_inputVarStr = "";
  setParameter("InputVar",&m_inputVarStr);
}

////////////////////////////////////////////////////////////////////////////

LEEInitState3D::~LEEInitState3D()
{
}

////////////////////////////////////////////////////////////////////////////

void LEEInitState3D::setup()
{

  CFAUTOTRACE;
  SpectralFDMethodCom::setup();

  // get number of physical variables
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim ();

  _varSet->setup();

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

  getMethodData().getUpdateVar().d_castTo<LinEuler3DVarSet>();

  //Get mean flow socket
  socket_meanflow.setParentNamespace( getMethodData().getNamespace() );

}

////////////////////////////////////////////////////////////////////////////
//

void LEEInitState3D::unsetup()
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

////////////////////////////////////////////////////////////////////////////
//

void LEEInitState3D::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  SpectralFDMethodCom::configure(args);

  // get the physical model that we are dealing with
  // to pass it to the variable transformer
  std::string namespc = getMethodData().getNamespace();
  SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(namespc);
  SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  // get the name of the update variable set
  std::string _updateVarStr = getMethodData().getUpdateVarStr();

  // create the transformer from input to update variables
  if (m_inputVarStr.empty()) {
    m_inputVarStr = _updateVarStr;
  }

  std::string provider =
    VarSetTransformer::getProviderName(physModel->getNameImplementor(), m_inputVarStr, _updateVarStr);

  m_inputToUpdateVar = Environment::Factory<VarSetTransformer>::getInstance().getProvider(provider)->create(physModel->getImplementor());
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

/***************************************************************************
******************************/

  std::string varSetName = "LinEuler3DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::LinearizedEuler::LinEuler3DCons>());

/***************************************************************************
******************************/

}

////////////////////////////////////////////////////////////////////////////
//

void LEEInitState3D::executeOnTrs()
{
  CFAUTOTRACE;

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");
  CFLogDebugMin("LEEInitState3D::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  // prepares to loop over cells by getting the GeometricEntityPool
  SafePtr< GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

  // get the physical model that we are dealing with
  // to pass it to the variable transformer
  std::string namespc = getMethodData().getNamespace();
  SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(namespc);
  SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();


/*
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
xxxxxxxxxxxx */
cout << "********** number of states:     " << trs->getNbStatesInTrs() << "\n" << flush;
cout << "********** number of nodes:      " << trs->getNbNodesInTrs() << "\n" << flush;
cout << "********** number of elements:   " << trs->getLocalNbGeoEnts() << "\n" << flush;
cout << "********** size of meanflow:     " << meanflow.size() << "\n" << flush;
cout << "********** sub size of meanflow: " << meanflow[0].size() << "\n" << flush;
// 0-2 xyz
// 3   localid
// 4   isnodal, 1 if states is node, -1 if only quadrature point
// 5-9 meanflow
std::vector< std::vector<CFreal> > dumpvec(trs->getNbStatesInTrs());
for(CFuint i=0; i<trs->getNbStatesInTrs(); i++) {
  dumpvec[i].resize(10,-1.);
  dumpvec[i][5]=meanflow[i][0];
  dumpvec[i][6]=meanflow[i][1];
  dumpvec[i][7]=meanflow[i][2];
  dumpvec[i][8]=meanflow[i][3];
  dumpvec[i][9]=meanflow[i][4];
}

for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
{
  const CFuint nbrElems = (*elemType)[iElemType].getNbElems();
  CFuint cellIdx = (*elemType)[iElemType].getStartIdx();
  SafePtr< vector< RealVector > > locInitPntCoords = sdLocalData[iElemType]->getInitPntsCoords();
  const CFuint nbrStates = sdLocalData[iElemType]->getNbrOfSolPnts();
  for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
  {
    geoData.idx = cellIdx;
    GeometricEntity *const cell = geoBuilder->buildGE();
    vector<State*>* solPntStates = cell->getStates();
    vector<Node*>* elmnodes=cell->getNodes();

    for (CFuint iPnt = 0; iPnt < nbrStates; ++iPnt)
    {
      CFuint lid=(*solPntStates)[iPnt]->getLocalID();
      m_initPntCoords = cell->computeCoordFromMappedCoord((*locInitPntCoords)[iPnt]);
      dumpvec[lid][0]=m_initPntCoords[0];
      dumpvec[lid][1]=m_initPntCoords[1];
      dumpvec[lid][2]=m_initPntCoords[2];
      dumpvec[lid][3]=(CFreal)lid;
      lid=(*elmnodes)[iPnt]->getLocalID();
      dumpvec[lid][4]=1.;
    }
    geoBuilder->releaseGE();
  }
}

ofstream fos("data.plt");
fos << "VARIABLES= x y z lid isnodal mA mB mC Md\n";

for(CFuint i=0; i<trs->getNbStatesInTrs(); i++) {
  fos << dumpvec[i][0] << " " << dumpvec[i][1] << " " << dumpvec[i][2] << " " << dumpvec[i][3] << " ";
  fos << dumpvec[i][4] << " " << dumpvec[i][5] << " " << dumpvec[i][6] << " " << dumpvec[i][7] << " ";
  fos << dumpvec[i][8] << " " << dumpvec[i][9] << "\n";
}

/*
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
xxxxxxxxxxxx */

  // loop over elements
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get the number of elements
    const CFuint nbrElems = (*elemType)[iElemType].getNbElems();

    // get start index of this element type in global element list
    CFuint cellIdx = (*elemType)[iElemType].getStartIdx();

    // local coordinates of initialization points and matrix used for initialization
    SafePtr< vector< RealVector > > locInitPntCoords = sdLocalData[iElemType]->getInitPntsCoords();
    SafePtr< RealMatrix > initTransfMatrix = sdLocalData[iElemType]->getInitTransfMatrix();

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


/***************************************************************************
**************************************/

  RealVector& linearData = _varSet->getModel()->getPhysicalData();

  RealVector _sumZ(m_nbrEqs);
  RealVector _avZ(m_nbrEqs);

  _sumZ = 0.0;

  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) {
    for (CFuint iState = 0; iState < nbrStates; ++iState) {
       CFuint IDstate = (*solPntStates)[iState]->getLocalID();
       RealVector meanflow_state = meanflow[IDstate];
       //CF_DEBUG_OBJ(meanflow_state);
      _sumZ[iEq] += meanflow_state[iEq];
    }
  }

  _avZ = _sumZ/static_cast<CFreal>(nbrStates);

  linearData[LinEulerTerm::GAMMA]  = _varSet->getModel()->getgamma();
  linearData[LinEulerTerm::rho0]   = _avZ[0];
  linearData[LinEulerTerm::U0]     = _avZ[1];
  linearData[LinEulerTerm::V0]     = _avZ[2];
  linearData[LinEulerTerm::W0]     = _avZ[3];
  linearData[LinEulerTerm::P0]     = _avZ[4];
  linearData[LinEulerTerm::c]      = sqrt(linearData[LinEulerTerm::GAMMA]*linearData[LinEulerTerm::P0]/linearData[LinEulerTerm::rho0]);

//CF_DEBUG_OBJ(linearData[LinEulerTerm::c]);

/***************************************************************************
**************************************/



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

////////////////////////////////////////////////////////////////////////////
//

  } // namespace SpectralFD

} // namespace COOLFluiD
