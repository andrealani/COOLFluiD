// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include<numeric>

#include "Common/PE.hh"
#include "Common/CFLog.hh"
#include "Framework/MeshDataBuilder.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/LocalConnectionData.hh"
#include "Common/BadValueException.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/SetElementStateCoord.hh"
#include "Framework/GeometricEntityRegister.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void MeshDataBuilder::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("PolynomialType","Polynomial type of the Geometric Entities to be built.");
}

//////////////////////////////////////////////////////////////////////////////

MeshDataBuilder::MeshDataBuilder(const std::string& name) :
  Common::OwnedObject(),
  Config::ConfigObject(name),
  m_cfmeshData(CFNULL),
  m_currBFaceID(0),
  m_mapGeoProviderNameToType()
{
  addConfigOptionsTo(this);

  m_polynomialTypeName = "Lagrange";
  setParameter( "PolynomialType", &m_polynomialTypeName );
}

//////////////////////////////////////////////////////////////////////////////

MeshDataBuilder::~MeshDataBuilder()
{
  releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void MeshDataBuilder::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);

  if (m_polynomialTypeName.empty())
    throw BadValueException (FromHere(),"No polynomial type selected in configuration of mesh data builder.\nSpecify 'PolynomialType' variable.");
}

//////////////////////////////////////////////////////////////////////////////

void MeshDataBuilder::releaseMemory()
{
  m_mapGeoProviderNameToType.clear();
}

//////////////////////////////////////////////////////////////////////////////

void MeshDataBuilder::setCFmeshData(SafePtr<CFmeshReaderSource> data)
{
  cf_assert(data.isNotNull());
  m_cfmeshData = data;
  m_cfmeshData->consistencyCheck();
}

//////////////////////////////////////////////////////////////////////////////

CFuint MeshDataBuilder::getNbElements() const
{
  return getCFmeshData().getNbElements();
}

//////////////////////////////////////////////////////////////////////////////

CFuint MeshDataBuilder::getNbElementTypes() const
{
  return getCFmeshData().getNbElementTypes();
}

//////////////////////////////////////////////////////////////////////////////

CFuint MeshDataBuilder::getNbNodes() const
{
  return getCFmeshData().getTotalNbNodes();
}

//////////////////////////////////////////////////////////////////////////////

CFuint MeshDataBuilder::getNbStates() const
{
  return getCFmeshData().getTotalNbStates();
}

//////////////////////////////////////////////////////////////////////////////

CFPolyOrder::Type MeshDataBuilder::getGeometricPolyOrder() const
{
  return getCFmeshData().getGeometricPolyOrder();
}

//////////////////////////////////////////////////////////////////////////////

CFPolyOrder::Type MeshDataBuilder::getSolutionPolyOrder() const
{
  return getCFmeshData().getSolutionPolyOrder();
}

//////////////////////////////////////////////////////////////////////////////

CFPolyForm::Type MeshDataBuilder::getGeometricPolyType() const
{
  return getCFmeshData().getGeometricPolyType();
}

//////////////////////////////////////////////////////////////////////////////

CFPolyForm::Type MeshDataBuilder::getSolutionPolyType() const
{
  return getCFmeshData().getSolutionPolyType();
}

//////////////////////////////////////////////////////////////////////////////

std::string MeshDataBuilder::makeGeomEntName(const CFGeoShape::Type& shape,
                            const CFPolyForm::Type&  geomPolyType,
                            const CFPolyOrder::Type& geomPolyOrder,
                            const CFPolyForm::Type&  solPolyType,
                            const CFPolyOrder::Type& solPolyOrder) const
{
  return
    CFGeoShape::Convert::to_str(shape) +
    CFPolyForm::Convert::to_str(geomPolyType) +
    CFPolyOrder::Convert::to_str(geomPolyOrder) +
    CFPolyForm::Convert::to_str(solPolyType) +
    CFPolyOrder::Convert::to_str(solPolyOrder);
}

//////////////////////////////////////////////////////////////////////////////

void MeshDataBuilder::createTopologicalRegionSets()
{
  CFAUTOTRACE;

  createInnerCells();

  setCoordInCellStates();

  createBoundaryTRS();

  setMapGeoToTrs();
}

//////////////////////////////////////////////////////////////////////////////

void MeshDataBuilder::createBoundaryTRS()
{
  const CFuint nbTRSs = getCFmeshData().getNbTRSs();

  SafePtr< vector<CFuint> > nbTRs =
    getCFmeshData().getNbTRs();

  SafePtr< vector< vector<CFuint> > > nbGeomEntsPerTR =
    getCFmeshData().getNbGeomEntsPerTR();

  SafePtr< vector<std::string> > nameTRS =
    getCFmeshData().getNameTRS();

  for(CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {

    cf_assert((*nameTRS)[iTRS] != "InnerCells");
    cf_assert((*nameTRS)[iTRS] != "InnerFaces");

    CFLogDebugMed("Nb TRs in TRS: " << (*nbGeomEntsPerTR)[iTRS].size() << "\n");

    Common::SafePtr<Framework::TopologicalRegionSet> ptrs =
    createTopologicalRegionSet((*nbGeomEntsPerTR)[iTRS],
                               (*nameTRS)[iTRS],
                               getCFmeshData().getTRGeoConn(iTRS),
                               iTRS);
    ptrs->attachTag("writable");
    ptrs->attachTag("boundary");
  }
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<TopologicalRegionSet> MeshDataBuilder::createTopologicalRegionSet
(const vector<CFuint>& nbGeomEntsPerTR,
 const std::string& name,
 const TRGeoConn& trGeoConn,
 const CFuint iTRS)
{
  CFAUTOTRACE;

  /// @todo TRGeoConn must be replaced by TopologicalRegion everywhere
  CFLog(INFO,"TRS Name: " << name << "\n");
  DataHandle < Framework::State*, Framework::GLOBAL > states =
    MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes =
    MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  // Create the TopologicalRegion storage
  const CFuint nbTRs = nbGeomEntsPerTR.size();
  vector<TopologicalRegion*>* storeTR = new vector<TopologicalRegion*>(nbTRs);

  // compute the total number of geometric entities
  const CFuint totalNbGeos = std::accumulate(nbGeomEntsPerTR.begin(),
                                        nbGeomEntsPerTR.end(),
                                        static_cast<CFuint>(0));
  // fill in two arrays specifying the number of states and nodes
  // per geometric entity in this TR
  std::valarray<CFuint> nbNodesInGeo(totalNbGeos);
  std::valarray<CFuint> nbStatesInGeo(totalNbGeos);
  setNbNodesNbStatesInGeo(nbGeomEntsPerTR, trGeoConn,
                          nbNodesInGeo, nbStatesInGeo);

  // allocate connectivity tables (you could allocate just 1 if you were in a case
  // when the node-connectivity is the same as the state connectivity)
  ConnectivityTable<CFuint>* geoNodes = new ConnectivityTable<CFuint>(nbNodesInGeo);
  ConnectivityTable<CFuint>* geoStates = new ConnectivityTable<CFuint>(nbStatesInGeo);

  // store the connectivities in the MeshData
  /// @todo think about having SharedPtr<ConnectivityTable<CFuint> > and not put
  /// the TRS connectivities in MeshData

  cf_always_assert (!name.empty());

  MeshDataStack::getActive()->storeConnectivity(name + "Nodes", geoNodes);
  MeshDataStack::getActive()->storeConnectivity(name + "States", geoStates);

  // array with all the IDs of all the geometric entities in this TRS
  // ownership of this array will be given to the TR
  vector<CFuint>* localGeoIDs = new vector<CFuint>(totalNbGeos);
  vector<CFuint>* globalGeoIDs = new vector<CFuint>(totalNbGeos);
  vector<CFuint>* geoTypes = new vector<CFuint>(totalNbGeos);
  CFuint nbProcessedGeoEnts = 0;

  const std::string faceProviderName = "Face";

  SafePtr<vector<vector<vector<CFuint> > > > trsGlobalIDs =
    MeshDataStack::getActive()->getGlobalTRSGeoIDs();

  const bool hasGlobalIDs = (trsGlobalIDs->size() > 0);

  // let's create the required number of TR's
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
    const GeoConn& geoConn = trGeoConn[iTR];
    const CFuint nbGeomEnts = nbGeomEntsPerTR[iTR];

    CFLog(VERBOSE, "iTR= " << iTR << ", nbProcessedGeoEnts= "
    << nbProcessedGeoEnts << "\n");

    // allocate the TopologicalRegion
    (*storeTR)[iTR] = new TopologicalRegion();
    (*storeTR)[iTR]->setLocalNbGeoEnts(nbGeomEnts);

    cf_assert(nbGeomEnts <= totalNbGeos);

    // check that the number of GeometricEntities in this TR is > 0
    // if not just skip all this
    if (nbGeomEnts > 0) {
      cf_assert(   (nbProcessedGeoEnts < localGeoIDs->size()) ||
                 ( (nbProcessedGeoEnts == localGeoIDs->size()) && ( localGeoIDs->size() == 0) ) );

      (*storeTR)[iTR]->setGeoEntsLocalIdx(&(*localGeoIDs)[nbProcessedGeoEnts],
            nbProcessedGeoEnts);

      // set the connectivity tables of the TRS
      for (CFuint iGeo = 0; iGeo < nbGeomEnts; ++iGeo, ++nbProcessedGeoEnts) {
  const CFuint nbGeoNodes = geoNodes->nbCols(nbProcessedGeoEnts);
  for (CFuint iNode = 0; iNode < nbGeoNodes; ++iNode) {
    const CFuint nodeID =  geoConn[iGeo].first[iNode];
    (*geoNodes)(nbProcessedGeoEnts, iNode) = nodeID;

    if (nodeID >= nodes.size() ) {
          CFLogDebugMax( "NodeID: " << nodeID << "\n");
          throw BadFormatException (FromHere(),"CFmesh had bad node index in GeometricEntity");
    }
  }

      const CFuint nbGeoStates = geoStates->nbCols(nbProcessedGeoEnts);

      for (CFuint iState = 0; iState < nbGeoStates; ++iState) {
        const CFuint stateID =  geoConn[iGeo].second[iState];
        (*geoStates)(nbProcessedGeoEnts, iState) = stateID;
        if (stateID >= states.size() ) {
          throw BadFormatException (FromHere(),"CFmesh had bad state index in GeometricEntity");
        }

  CFLogDebugMed(name << ", iTR= " << iTR << ", iGeo= " << iGeo
                << ", iState= " << iState << ", globalID= " <<
          states[stateID]->getGlobalID() << "\n");
        //CFLogDebugMed("state coord= " << states[stateID]->getCoordinates() << "\n");
      }

      // assign the local ID of the current geometric entity
      (*localGeoIDs)[nbProcessedGeoEnts] = m_currBFaceID++;

  // assign the global ID of the current geometric entity
  (*globalGeoIDs)[nbProcessedGeoEnts] = (hasGlobalIDs) ?
    (*trsGlobalIDs)[iTRS][iTR][iGeo] : iGeo; // watch out here !!!!

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
 const std::string faceGeoTypeName = makeGeomEntName (getGeoShape(CFGeoEnt::FACE, dim, nbGeoNodes),
                                                    getGeometricPolyType(),
                                                    getGeometricPolyOrder(),
                                                    getSolutionPolyType(),
                                                    getSolutionPolyOrder());

      const std::string geoProviderName = faceProviderName + faceGeoTypeName;
      (*geoTypes)[nbProcessedGeoEnts] = m_mapGeoProviderNameToType.find(geoProviderName);
      }
    }
  }

  // create TopologicalRegionSet
  TopologicalRegionSet* ptrs = new TopologicalRegionSet(name, storeTR);
  // we know it is a TRS of faces
  ptrs->attachTag("face");

  // set local GeometricEntity IDs
  ptrs->setGeoEntsLocalIdx(localGeoIDs);
  // set global GeometricEntity IDs
  ptrs->setGeoEntsGlobalIdx(globalGeoIDs);
  // set the connectivity GeometricEntity to nodeIDs
  ptrs->setGeo2NodesConn(geoNodes);
  // set the connectivity GeometricEntity to stateIDs
  ptrs->setGeo2StatesConn(geoStates);
  // set the GeometricEntity types
  ptrs->setGeoTypes(geoTypes);
    // create the cached list of state indexes in the all TRS
  ptrs->createStatesList();
  // put it into MeshData TRS list
  MeshDataStack::getActive()->addTrs(ptrs);

  cf_assert(nbProcessedGeoEnts == totalNbGeos);

  // increment nb faces in statistics
  const CFuint new_nb_faces = MeshDataStack::getActive()->Statistics().getNbFaces() + nbProcessedGeoEnts;
  MeshDataStack::getActive()->Statistics().setNbFaces(new_nb_faces);

  CFLog(INFO, "Created Boundary TRS \'" << name << "\' with " << nbProcessedGeoEnts << " faces\n");

  return ptrs;
}

//////////////////////////////////////////////////////////////////////////////

void MeshDataBuilder::setElementTypeData()
{
  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  // deallocate any used memory for the element type data
  vector<ElementTypeData>().swap(*elementType);

  SafePtr<vector<ElementTypeData> > elementTypeRead =
    getCFmeshData().getElementTypeData();

  // copy to MeshData the element type data read from file (stored in CFmeshData)
  *elementType = *elementTypeRead;

  CFuint nbElemTypes = getNbElementTypes();
  cf_assert(nbElemTypes == elementType->size());

  CFuint startIdx = 0;
  for (CFuint iType = 0; iType < nbElemTypes; ++iType)
  {
    (*elementType)[iType].setStartIdx( startIdx );
    (*elementType)[iType].setGeoOrder( getCFmeshData().getGeometricPolyOrder() );
    (*elementType)[iType].setSolOrder( getCFmeshData().getSolutionPolyOrder()  );
    startIdx += (*elementType)[iType].getNbElems();
  }
}

//////////////////////////////////////////////////////////////////////////////

void MeshDataBuilder::setNbNodesNbStatesInGeo
(const vector<CFuint>& nbGeomEntsPerTR,
 const TRGeoConn& trGeoConn,
 std::valarray<CFuint>& nbNodesInGeo,
 std::valarray<CFuint>& nbStatesInGeo)
{
  cf_assert(nbStatesInGeo.size() == nbNodesInGeo.size());

  CFuint geoID = 0;
  const CFuint nbTRs = nbGeomEntsPerTR.size();
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
    const CFuint nbGeomEnts = nbGeomEntsPerTR[iTR];
    const GeoConn& geoConn = trGeoConn[iTR];
    for (CFuint iGeo = 0; iGeo < nbGeomEnts; ++iGeo, ++geoID) {
      nbNodesInGeo[geoID] = geoConn[iGeo].first.size();
      nbStatesInGeo[geoID] = geoConn[iGeo].second.size();

      if(nbNodesInGeo[geoID] < 1) {
        throw BadFormatException
          (FromHere(),"MeshDataBuilder::setNbNodesNbStatesInGeo()=> wrong nbNodesInGeo");
      }
    }
  }

  cf_assert(geoID == nbStatesInGeo.size());
}

//////////////////////////////////////////////////////////////////////////////

void MeshDataBuilder::computeGeoTypeInfo()
{
  CFAUTOTRACE;

  // list of the geometric polynomial order types in the elements of the mesh
  vector<CFPolyOrder::Type> geomOrders ( getNbElementTypes(), getGeometricPolyOrder());

  // list of the solution polynomial order types in the elements of the mesh
  vector<CFPolyOrder::Type> solOrders ( getNbElementTypes(), getSolutionPolyOrder());

  // copy the element type data in the MeshData
  setElementTypeData();

  SafePtr<vector<ElementTypeData> > elementType = getCFmeshData().getElementTypeData();

  const CFuint nbElemTypes = elementType->size();
  vector<CFGeoShape::Type> elementGeoShapes(nbElemTypes);
  for (CFuint i = 0; i < nbElemTypes; ++i) {
    elementGeoShapes[i] = (*elementType)[i].getGeoShape();
  }

  // here add face types and "register" them

  vector< vector<CFGeoShape::Type> > faceShapesPerElemType(nbElemTypes);
  LocalConnectionData::getInstance().setFaceShapesPerElemType(faceShapesPerElemType);

  set<std::string> geoProviderNames;

  for (CFuint iType = 0; iType < nbElemTypes; ++iType)
  {
    const std::string cellGeoTypeName = makeGeomEntName
      ((*elementType)[iType].getGeoShape(),
       getGeometricPolyType(),
       getGeometricPolyOrder(),
       getSolutionPolyType(),
       getSolutionPolyOrder());

    const std::string cellProviderName = "Cell" + cellGeoTypeName;

    geoProviderNames.insert(cellProviderName);

    const CFuint nbElemFaces = faceShapesPerElemType[iType].size();
    for (CFuint iFace = 0; iFace < nbElemFaces; ++iFace)
    {
      const std::string faceGeoTypeName =
        makeGeomEntName(faceShapesPerElemType[iType][iFace],
                        getGeometricPolyType(),
                        getGeometricPolyOrder(),
                        getSolutionPolyType(),
                        getSolutionPolyOrder());

      const std::string faceProviderName = "Face" + faceGeoTypeName;

      geoProviderNames.insert(faceProviderName);
    }
  }

  // regist all the found GeometricEntity's provider names that will be used by the simulation
  // store locally (in this class) the mapping GeometricEntity's provider
  // name to corresponding geometric type (which is an integer) to reuse
  // it later while constructing the TRSs
  set<std::string>::const_iterator it;
  for (it = geoProviderNames.begin(); it != geoProviderNames.end(); ++it)
  {
    m_mapGeoProviderNameToType.insert(*it, GeometricEntityRegister::getInstance().regist(*it));
  }
  m_mapGeoProviderNameToType.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void MeshDataBuilder::setCoordInCellStates()
{
  CFAUTOTRACE;

  // gets the nodes and states
  DataHandle<State*,GLOBAL> states = getCFmeshData().getStatesHandle();
  DataHandle<Node*,GLOBAL> nodes = getCFmeshData().getNodesHandle();

  // get the TRS of Cells
  std::vector< Common::SafePtr<TopologicalRegionSet> > cells = MeshDataStack::getActive()->getTrsList();
  for(CFuint iTrs=0; iTrs < cells.size(); ++iTrs)
  {
    if(cells[iTrs]->hasTag("inner") && cells[iTrs]->getName() != "InnerFaces"){
      //       cf_assert(getNbElements() == cells[iTrs]->getLocalNbGeoEnts());
      
      // get element type data
      SafePtr<vector<ElementTypeData> > elementType =
        MeshDataStack::getActive()->getElementTypeData(cells[iTrs]->getName());

      const CFuint nbElemTypes = getNbElementTypes();
      cf_assert(nbElemTypes == elementType->size());
      
      // loop over the element types
      CFuint elemID = 0;
      for (CFuint iType = 0; iType < nbElemTypes; ++iType) {

        std::string elemName = makeGeomEntName
          ((*elementType)[iType].getGeoShape(),
          getGeometricPolyType(),
          getGeometricPolyOrder(),
          getSolutionPolyType(),
          getSolutionPolyOrder());
	
        SelfRegistPtr<SetElementStateCoord>* setStateCoord =
          FACTORY_GET_PROVIDER(getFactoryRegistry(), SetElementStateCoord, elemName)->
	  createPtr();
	
	(*setStateCoord)->setFactoryRegistry(getFactoryRegistry());
        const CFuint nbNodesPerElem = (*elementType)[iType].getNbNodes();
        const CFuint nbStatesPerElem = (*elementType)[iType].getNbStates();
	
        vector<Node*> eNodes(nbNodesPerElem);
        vector<State*> eStates(nbStatesPerElem);

        // loop on all the elements in this type
        const CFuint nbElemPerType = (*elementType)[iType].getNbElems();

        for (CFuint iElem = 0; iElem < nbElemPerType; ++iElem, ++elemID)
        {
          cf_assert(nbNodesPerElem == cells[iTrs]->getNbNodesInGeo(elemID));
          cf_assert(nbStatesPerElem == cells[iTrs]->getNbStatesInGeo(elemID));

          // copy the Node pointers into the vector
          for (CFuint iNode = 0; iNode < nbNodesPerElem; ++iNode)
          {
            const CFuint nodeID = cells[iTrs]->getNodeID(elemID,iNode);
//     CFout << "NodeID: " << nodeID << "\n";
            Node* nodePtr = nodes[nodeID];
            eNodes[iNode] = nodePtr;
          }

          // copy the State pointers into the vector
          for (CFuint iState = 0; iState < nbStatesPerElem; ++iState) {
            const CFuint stateID = cells[iTrs]->getStateID(elemID,iState);
//     CFout << "StateID: " << stateID << "\n";
            State* statePtr = states[stateID];
            eStates[iState] = statePtr;
          }

          // finally set the coordinates in the State accordingly
          // to the algorithm of the element type
          (*(*setStateCoord))(eNodes, eStates);
        }
        delete setStateCoord;
      }

//       cf_assert(elemID == getNbElements());
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFGeoShape::Type MeshDataBuilder::getGeoShape(CFGeoEnt::Type  geomType,
            const CFuint& dim,
                              const CFuint& nbGeoNodes) const
{
  switch(dim)
  {
  case(DIM_2D):
  switch(nbGeoNodes)
  {
  case(2): return CFGeoShape::LINE;
  case(3):
    if (geomType == CFGeoEnt::FACE) return CFGeoShape::LINE;
    if (geomType == CFGeoEnt::CELL) return CFGeoShape::TRIAG;
    throw BadValueException (FromHere(),"MeshDataBuilder::getGeoShape() bad geomType [" + Common::StringOps::to_str(geomType) + "]" );
  case(4):
    if (geomType == CFGeoEnt::FACE) return CFGeoShape::LINE;
    if (geomType == CFGeoEnt::CELL) return CFGeoShape::QUAD;
    throw BadValueException (FromHere(),"MeshDataBuilder::getGeoShape() bad geomType [" + Common::StringOps::to_str(geomType) + "]" );
  default:
    throw BadValueException (FromHere(),"MeshDataBuilder::getGeoShape() nb nodes [" + Common::StringOps::to_str(nbGeoNodes) + "]" );
  }
    break;
  case(DIM_3D):
    switch(nbGeoNodes) {
    case(3):
      if (geomType == CFGeoEnt::FACE) return CFGeoShape::TRIAG;
      throw BadValueException (FromHere(),"MeshDataBuilder::getGeoShape() bad geomType [" + Common::StringOps::to_str(geomType) + "]" );
    case(4):
      if (geomType == CFGeoEnt::FACE) return CFGeoShape::QUAD;
      if (geomType == CFGeoEnt::CELL) return CFGeoShape::TETRA;
    case(5): return CFGeoShape::PYRAM;
    case(6):
	if (geomType == CFGeoEnt::FACE) return CFGeoShape::TRIAG;
	if (geomType == CFGeoEnt::CELL)  return CFGeoShape::PRISM;
    case(8): return CFGeoShape::HEXA;
    default:
      throw BadValueException (FromHere(),"MeshDataBuilder::getGeoShape() unknown geometric entity type");
    }
    break;
  case(DIM_1D):
    switch(nbGeoNodes)
    {
    case(1):
      if (geomType == CFGeoEnt::FACE) return CFGeoShape::POINT;
      throw BadValueException (FromHere(),"MeshDataBuilder::getGeoShape() bad geomType [" + Common::StringOps::to_str(geomType) + "]" );
    case(2):
      if (geomType == CFGeoEnt::CELL) return CFGeoShape::LINE;
      throw BadValueException (FromHere(),"MeshDataBuilder::getGeoShape() bad geomType [" + Common::StringOps::to_str(geomType) + "]" );
    default:
      throw BadValueException (FromHere(),"MeshDataBuilder::getGeoShape() nb nodes [" + Common::StringOps::to_str(nbGeoNodes) + "]" );
    }
    break;
  default:
    throw BadValueException (FromHere(),"MeshDataBuilder::getGeoShape() bad dimension [" + Common::StringOps::to_str(dim) + "]" );
  }
}

//////////////////////////////////////////////////////////////////////////////

void MeshDataBuilder::createInnerCells()
{
  CFAUTOTRACE;

  const CFuint nbGroups = getCFmeshData().getNbGroups();

  Common::SafePtr< vector<std::string> > groupNames = getCFmeshData().getGroupNames();

  Common::SafePtr<MeshData::ConnTable> cellNodesConn = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
  Common::SafePtr<MeshData::ConnTable> cellStatesConn = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");
//
  // create the list of TopologicalRegions
  vector<TopologicalRegion*>* cellTrList(CFNULL);

  // create the TopologicalRegionSet "InnerCells"
  TopologicalRegionSet* ptrTRS(CFNULL);
  for(CFuint iGroup = 0; iGroup < nbGroups; ++iGroup)
  {
    const CFuint nbElem = (*getCFmeshData().getGroupSizes())[iGroup];

    vector<CFuint>* cellLocalIDs  = new vector<CFuint>(nbElem);
    vector<CFuint>* cellGlobalIDs = new vector<CFuint>(nbElem);
    vector<CFuint>* cellGeoTypes  = new vector<CFuint>(nbElem);

    cellTrList = new vector<TopologicalRegion*>(1);
    (*cellTrList)[0] = new TopologicalRegion();

    const std::string TRSname = (*groupNames)[iGroup];
    ptrTRS = new TopologicalRegionSet (TRSname, cellTrList);

    // this is a inner TRS
    ptrTRS->attachTag("inner");
    // we know it is a TRS of cells
    ptrTRS->attachTag("cell");

    // allocate element type data
    MeshDataStack::getActive()->createElementTypeData(TRSname);
    // get element type data
    SafePtr<vector<ElementTypeData> > trsElementType =
      MeshDataStack::getActive()->getElementTypeData(TRSname);

    SafePtr< vector<ElementTypeData> > elementType = MeshDataStack::getActive()->getElementTypeData();

    // copy the global element type data to the one local to this TRS
    // to be modified later in this function
    *trsElementType = *elementType;

    SafePtr< vector<CFuint> > globalElementIDs = MeshDataStack::getActive()->getGlobalElementIDs();
    const bool hasGlobalID = (globalElementIDs->size() > 0);

    const CFuint nbElemTypes = trsElementType->size();
    for (CFuint iType = 0; iType < nbElemTypes; ++iType)
    {
      (*trsElementType)[iType].setNbElems(0);
    }


    std::valarray<CFuint> elemNodePattern(nbElem);
    std::valarray<CFuint> elemStatePattern(nbElem);

    CFuint elemID = 0;
    for(CFuint iElem=0; iElem < nbElem; ++iElem, ++elemID)
    {

      const CFuint localCellID = (*getCFmeshData().getGroupElementLists())[iGroup][iElem];

      //find what is the type of the element
      CFuint rangeMin = 0;
      CFuint rangeMax = 0;
      CFuint cellType = 0;
      bool found(false);
      for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
        const CFuint nbElemPerType = (*elementType)[iType].getNbElems();
        rangeMax += nbElemPerType;
// CFout << "Checking cell iCell: " << iCell << " - LocalID: " << localCellID <<  " <<"\n"
        if((localCellID < rangeMax) && (localCellID >= rangeMin)) {
          cellType = iType;
          found = true;
        }
        rangeMin += nbElemPerType;
      }

      cf_assert(found);
      (*trsElementType)[cellType].setNbElems((*trsElementType)[cellType].getNbElems() + 1);

      const std::string cellGeoTypeName = makeGeomEntName
         ((*elementType)[cellType].getGeoShape(),
         getGeometricPolyType(),
         getGeometricPolyOrder(),
         getSolutionPolyType(),
         getSolutionPolyOrder());
      const std::string cellProviderName = "Cell" + cellGeoTypeName;

      // place the provider ID in a map that the GeometricEntityBuilder will latter use
      (*cellLocalIDs) [iElem]  = localCellID;
      (*cellGlobalIDs)[iElem]  = (hasGlobalID) ? (*globalElementIDs)[localCellID] : localCellID;
      (*cellGeoTypes) [iElem]  = m_mapGeoProviderNameToType.find(cellProviderName);

      elemNodePattern[iElem] = cellNodesConn->nbCols(localCellID);
      elemStatePattern[iElem] = cellStatesConn->nbCols(localCellID);

    }

    CFuint startidx = 0;
    for(CFuint iType=0; iType < trsElementType->size(); ++iType){
      (*trsElementType)[iType].setStartIdx(startidx);
      startidx += (*trsElementType)[iType].getNbElems();

/*  CFout << "nbElements[" <<iType <<": " << (*trsElementType)[iType].getNbElems() <<"\n";
  CFout << "getNbNodes[" <<iType <<": " << (*trsElementType)[iType].getNbNodes() <<"\n";
  CFout << "getNbStates[" <<iType <<": " << (*trsElementType)[iType].getNbStates() <<"\n";*/
    }

    // place the data in the TRS "InnerCells"
    ptrTRS->setGeoTypes(cellGeoTypes);
    ptrTRS->setGeoEntsLocalIdx(cellLocalIDs);

    // place the data for the only exixting TR
    (*cellTrList)[0]->setLocalNbGeoEnts(nbElem);
    (*cellTrList)[0]->setGeoEntsLocalIdx(&(*cellLocalIDs)[0], 0);


    if(nbGroups > 1)
    {
      //take care of the new connectivity
      ConnectivityTable<CFuint>* newCellNodes = new
        ConnectivityTable<CFuint>(elemNodePattern);

      ConnectivityTable<CFuint>* newCellStates = new
        ConnectivityTable<CFuint>(elemStatePattern);

      MeshDataStack::getActive()->storeConnectivity("cellNodes_" + ptrTRS->getName(), newCellNodes);
      MeshDataStack::getActive()->storeConnectivity("cellStates_" + ptrTRS->getName(), newCellStates);

      ptrTRS->setGeo2NodesConn(MeshDataStack::getActive()->getConnectivity("cellNodes_" + ptrTRS->getName()));
      ptrTRS->setGeo2StatesConn(MeshDataStack::getActive()->getConnectivity("cellStates_" + ptrTRS->getName()));

      for(CFuint iElem=0; iElem < nbElem; ++iElem){
        const CFuint localCellID = ptrTRS->getLocalGeoID(iElem);
        const CFuint nbNodesInCell = elemNodePattern[iElem];
        const CFuint nbStatesInCell = elemStatePattern[iElem];
        cf_assert(nbNodesInCell == cellNodesConn->nbCols(localCellID));
        cf_assert(nbStatesInCell == cellStatesConn->nbCols(localCellID));

        // renumber the node connectivity of the current cell
        for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
          const CFuint nodeID = (*cellNodesConn)(localCellID, iNode);
          ptrTRS->setNodeID(iElem, iNode, nodeID);
        }
        for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
          const CFuint stateID = (*cellStatesConn)(localCellID, iState);
          ptrTRS->setStateID(iElem, iState, stateID);
        }
      }
    }
    else
    {
      ptrTRS->setGeo2NodesConn(cellNodesConn);
      ptrTRS->setGeo2StatesConn(cellStatesConn);
    }

    /// @todo this should be checked if it is working in paralell
    ptrTRS->setGlobalNbGeoEnts(cellGlobalIDs->size());
    ptrTRS->setGeoEntsGlobalIdx(cellGlobalIDs);
    ptrTRS->createStatesList();

    // add the TRS to the MeshData
    MeshDataStack::getActive()->addTrs(ptrTRS);

  }

  getCFmeshData().setNbGroups(0);

  if(nbGroups > 1)
  {
    MeshDataStack::getActive()->removeConnectivity("cellNodes_InnerCells");
    MeshDataStack::getActive()->removeConnectivity("cellStates_InnerCells");
  }
//   // check that the number of global element IDs match the number of elements
//   if (hasGlobalID)
//   {
//     if (globalElementIDs->size() != nbElem)
//     {
//       throw BadFormatException
//          ("MeshDataBuilder::createInnerCells() : globalElementIDs->size() != nbElements");
//     }
//   }
//     else
//   {
//     CFLogDebugMed("No global element IDs!\n");
//
//   }

}

//////////////////////////////////////////////////////////////////////////////

CFuint MeshDataBuilder::getNbBoundaryFaces() const
{
  const CFuint nbTRSs = getCFmeshData().getNbTRSs();
  SafePtr< const vector<CFuint> > nbTRs = getCFmeshData().getNbTRs();

  SafePtr< const vector< vector<CFuint> > > nbGeomEntsPerTR = getCFmeshData().getNbGeomEntsPerTR();

  CFuint nbBFaces = 0;
  for(CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    nbBFaces += std::accumulate((*nbGeomEntsPerTR)[iTRS].begin(),
                           (*nbGeomEntsPerTR)[iTRS].end(),
                           static_cast<CFuint>(0));
  }
  return nbBFaces;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

