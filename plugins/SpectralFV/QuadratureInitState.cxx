#include "SpectralFV/SpectralFV.hh"


#include "SpectralFV/QuadratureInitState.hh"
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

MethodCommandProvider<QuadratureInitState, SpectralFVMethodData, SpectralFVModule>
QuadratureInitStateProvider("QuadratureInitState");

//////////////////////////////////////////////////////////////////////////////

void QuadratureInitState::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::string >("InputVar","Input variables.");
}

//////////////////////////////////////////////////////////////////////////////

QuadratureInitState::QuadratureInitState(const std::string& name) :
  SpectralFVMethodCom(name),
  socket_states("states"),
  m_varSet(CFNULL),
  m_inputToUpdateVar(),
  m_inputVars(CFNULL),
  m_nbrEqs()
{
  addConfigOptionsTo(this);
  m_functions = std::vector<std::string>();
  setParameter("Def",&m_functions);

  m_vars = std::vector<std::string>();
  setParameter("Vars",&m_vars);

  m_inputVarStr = "";
  setParameter("InputVar",&m_inputVarStr);
}

//////////////////////////////////////////////////////////////////////////////

QuadratureInitState::~QuadratureInitState()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
QuadratureInitState::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void QuadratureInitState::setup()
{

  CFAUTOTRACE;
  SpectralFVMethodCom::setup();

  // get number of physical variables
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  m_inputVars = new State();

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  m_inputToUpdateVar->setup(maxNbStatesInCell);

  m_varSet = getMethodData().getUpdateVar();
}

//////////////////////////////////////////////////////////////////////////////

void QuadratureInitState::unsetup()
{
  CFAUTOTRACE;

  deletePtr(m_inputVars);

  SpectralFVMethodCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void QuadratureInitState::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  SpectralFVMethodCom::configure(args);

  // get the physical model that we are dealing with
  // to pass it to the variable transformer
  std::string namespc = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(namespc);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

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

void QuadratureInitState::executeOnTrs()
{
  CFAUTOTRACE;

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs = MeshDataStack::getActive()->getTrs("InnerCells");
  CFLogDebugMin("QuadratureInitState::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  // prepares to loop over cells by getting the GeometricEntityPool
  SafePtr< GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

  // get the cell - state connectivity
  SafePtr<MeshData::ConnTable> cellStates = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

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
    SafePtr< vector< CFreal > >           volFracCV = svLocalData[iElemType]->getVolFracCV();
    SafePtr< vector< RealVector > >       localNodeCoord = svLocalData[iElemType]->getLocalNodeCoord();
    SafePtr< vector< vector< CFuint > > > localCVNodeConn = svLocalData[iElemType]->getLocalCVNodeConn();

    // get element shape
    const CFGeoShape::Type shape = svLocalData[iElemType]->getShape();

    // get dimensionality
    const CFDim dimensionality = svLocalData[iElemType]->getDimensionality();

    // get polynomial order
    const CFPolyOrder::Type polyOrder = svLocalData[iElemType]->getPolyOrder();

    // get number of control volumes
    const CFuint nbrCVs = svLocalData[iElemType]->getNbrOfCVs();

    // set the simplex integrator
    m_sIntegrator.setDimensionality(dimensionality);
    m_sIntegrator.setIntegratorOrder(polyOrder);
    const CFuint nbrSimplexQuadPnts = m_sIntegrator.getNbrQuadPnts();

    switch (shape)
    {
      case CFGeoShape::LINE:
      {
        for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
        {
          // build the GeometricEntity
          geoData.idx = cellIdx;
          GeometricEntity *const cell = geoBuilder->buildGE();

          // get the SV nodes
          vector<Node*>* svNodes = cell->getNodes();

          // check if element is a line
          cf_assert(cell->getShape() == shape);
          cf_assert(svNodes->size() == 2);

          // loop over control volumes in this spectral volume
          for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
          {
            // number of nodes in this cv
            const CFuint nbrCVNodes = (*localCVNodeConn)[iCV].size();

            // compute coordinates of the nodes in this CV
            vector< RealVector > cvNodeCoords(nbrCVNodes);
            for (CFuint iCVNode = 0; nbrCVNodes; ++iCVNode)
            {
              // resize
              cvNodeCoords[iCV].resize(dimensionality);

              // get local node ID
              const CFuint localNodeID = (*localCVNodeConn)[iCV][iCVNode];

              // compute CV node coordinates
              cvNodeCoords[iCV] = (1.0 - (*localNodeCoord)[localNodeID].sum())*(*(*svNodes)[0]);
              for (CFuint iCoor = 0; iCoor < (CFuint) dimensionality; ++iCoor)
              {
                cvNodeCoords[iCV] += (*localNodeCoord)[localNodeID][iCoor]*(*(*svNodes)[iCoor+1]);
              }
            }

            // get quadrature data
            vector< RealVector >  qNodeCoord  = m_sIntegrator.getQuadPntsCoords(cvNodeCoords);
            vector< CFreal >      qWheights   = m_sIntegrator.getQuadPntsWheights(cvNodeCoords);

            // variable for a dimensional state
            State dimState(RealVector(0.0,m_nbrEqs));

            // compute integral over control volumes
            for (CFuint iQPnt = 0; iQPnt < nbrSimplexQuadPnts; ++iQPnt)
            {
              // evaluate the function at the state coordinate
              m_vFunction.evaluate(qNodeCoord[iQPnt],*m_inputVars);

              // transform from input variables to update variables
              dimState += qWheights[iQPnt]*(*m_inputToUpdateVar->transform(m_inputVars));
            }

            // divide by control volume volume
            dimState /= (*volFracCV)[iCV]*cell->computeVolume();

            // get the state pointer
            State *const state = states[(*cellStates)(cellIdx,iCV)];

            // adimensionalize the value if needed and store
            m_varSet->setAdimensionalValues(dimState, *state);

          }

          //release the GeometricEntity
          geoBuilder->releaseGE();
        }
      } break;
      case CFGeoShape::TRIAG:
      {
        for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
        {
          // build the GeometricEntity
          geoData.idx = cellIdx;
          GeometricEntity *const cell = geoBuilder->buildGE();

          // get the SV nodes
          vector<Node*>* svNodes = cell->getNodes();

          // check if element is a line
          cf_assert(cell->getShape() == shape);
          cf_assert(svNodes->size() == 3);

          // loop over control volumes in this spectral volume
          for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
          {
            // number of nodes in this cv
            const CFuint nbrCVNodes = (*localCVNodeConn)[iCV].size();
            vector< RealVector > cvNodeCoords(nbrCVNodes);
            for (CFuint iCVNode = 0; iCVNode < nbrCVNodes; ++iCVNode)
            {
              // resize
              cvNodeCoords[iCVNode].resize(dimensionality);

              // get local node ID
              const CFuint localNodeID = (*localCVNodeConn)[iCV][iCVNode];

              // compute CV node coordinates
              cvNodeCoords[iCVNode] = (1.0 - (*localNodeCoord)[localNodeID].sum())*(*(*svNodes)[0]);
              for (CFuint iCoor = 0; iCoor < (CFuint) dimensionality; ++iCoor)
              {
                cvNodeCoords[iCVNode] += (*localNodeCoord)[localNodeID][iCoor]*(*(*svNodes)[iCoor+1]);
              }
            }

            // variables for quadrature node coordinates and wheights
            vector< RealVector >  qNodeCoord;
            vector< CFreal >      qWheights;

            // compute quadrature data for this this CV
            const CFuint nbrTriangles = nbrCVNodes-2;
            for (CFuint iTriangle = 0; iTriangle < nbrTriangles; ++iTriangle)
            {
              // triangle node coordinates variable
              vector< RealVector > triagNodeCoords(3);
              for (CFuint iTriagNode = 0; iTriagNode < 3; ++iTriagNode)
              {
                triagNodeCoords[iTriagNode].resize(dimensionality);
              }

              // assign triangle nodes
              triagNodeCoords[0] = cvNodeCoords[0];
              triagNodeCoords[1] = cvNodeCoords[iTriangle+1];
              triagNodeCoords[2] = cvNodeCoords[iTriangle+2];

              // get triangle quadrature nodes and wheights
              vector< RealVector > qNodeCoordTriangle = m_sIntegrator.getQuadPntsCoords(triagNodeCoords);
              vector< CFreal > qWheightsTriangle = m_sIntegrator.getQuadPntsWheights(triagNodeCoords);

              // add triangle quadrature nodes and wheights to global list
              qNodeCoord.insert(qNodeCoord.end(),qNodeCoordTriangle.begin(),qNodeCoordTriangle.end());
              qWheights.insert(qWheights.end(),qWheightsTriangle.begin(),qWheightsTriangle.end());
            }

            // number of quadrature points
            const CFuint nQuadPnts = qWheights.size();

            // variable for a dimensional state
            State dimState(RealVector(0.0,m_nbrEqs));

            // compute integral over control volumes
            for (CFuint iQPnt = 0; iQPnt < nQuadPnts; ++iQPnt)
            {
              // evaluate the function at the state coordinate
              m_vFunction.evaluate(qNodeCoord[iQPnt],*m_inputVars);

              // transform from input variables to update variables
              dimState += qWheights[iQPnt]*(*m_inputToUpdateVar->transform(m_inputVars));
            }

            // divide by control volume volume
            dimState /= (*volFracCV)[iCV]*cell->computeVolume();

            // get the state pointer
            State *const state = states[(*cellStates)(cellIdx,iCV)];

            // adimensionalize the value if needed and store
            m_varSet->setAdimensionalValues(dimState, *state);
          }

          //release the GeometricEntity
          geoBuilder->releaseGE();
        }
      } break;
      default:
        throw Common::NotImplementedException (FromHere(),"Only linear and triangular cells have been implemented for spectral FV!");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD
