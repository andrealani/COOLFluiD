// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <numeric>

#include "Framework/CFmeshReaderSource.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

CFmeshReaderSource::CFmeshReaderSource() :
  socket_nodes("Null"),
  socket_states("Null"),
  _elementNode(CFNULL),
  _elementState(CFNULL),
  _nodalExtraDataSourceSockets(),
  _stateExtraDataSourceSockets(),
  _extraDataSourceSockets(),
  _totalNbStateExtraVars(0),
  _totalNbNodalExtraVars(0),
  _totalNbExtraVars(0)
{
}

//////////////////////////////////////////////////////////////////////////////

CFmeshReaderSource::~CFmeshReaderSource()
{
  releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::releaseMemory()
{
  BaseCFMeshFileSource::releaseMemory();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::setDataSockets(DataSocketSink<State*,GLOBAL> statesSocket,
					DataSocketSink<Node*,GLOBAL> nodesSocket,
                                        Common::SafePtr<DynamicDataSocketSet<> > dynamicSocket)
{
  // set the sockets for the needed storages
  socket_nodes = nodesSocket;
  socket_states = statesSocket;
  dynamicSockets = dynamicSocket;

  // gets the connectivity of nodes and states
  _elementNode = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
  _elementState = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");
}

//////////////////////////////////////////////////////////////////////////////

///@todo this is to be removed later but for compatibility with new mesh reader (hdf5)
void CFmeshReaderSource::setDataSockets(DataSocketSink<State*,GLOBAL> statesSocket,
					DataSocketSink<Node*,GLOBAL> nodesSocket)
{
  // set the sockets for the needed storages
  socket_nodes = nodesSocket;
  socket_states = statesSocket;

  // gets the connectivity of nodes and states
  _elementNode = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
  _elementState = MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::resizeNodes(const CFuint nbNodes)
{
  if (socket_nodes.getDataHandle().size() > 0) {
    const CFuint nbNodes = socket_nodes.getDataHandle().size();
    for (CFuint i = 0; i < nbNodes; ++i) {
      deletePtr(socket_nodes.getDataHandle()[i]);
    }
    // here resizing to 0 can lead to seg fault with some compilers
    socket_nodes.getDataHandle().resize(1);
  }
  
  cf_assert(socket_nodes.getDataHandle().size() == 0 || 
	    socket_nodes.getDataHandle().size() == 1 );
  socket_nodes.getDataHandle().resize(nbNodes);
  IndexList<Node>::getList().reset();
  cf_assert(socket_nodes.getDataHandle().size() == nbNodes);
  _NodeLocalToGlobal.resize (nbNodes);
  _NodeOwnership.resize(nbNodes);
  
  //Resize the pastNodes Datahandle
  if(_storePastNodes){
    std::string socketName = "pastNodes";
    DataHandle<Node*> pastNodes = dynamicSockets->getSocketSink<Node*>(socketName)->getDataHandle();
    pastNodes.resize(nbNodes);
  }

  //Setup for the Extra Vars associated with the nodes
  const CFuint nbExtraVarsInMesh = getNbExtraNodalVars();
  if((nbExtraVarsInMesh > 0) && (_extraNVarMap.isNotNull()))
  {
    const CFuint nbExtraVarsInCFcase = _extraNVarMap->size();
    cf_assert(nbExtraVarsInCFcase <= nbExtraVarsInMesh);

    ///resizing the nodal extra vars datahandles
    CFuint nbVarFound = 0;
    _extraNVarExists.resize(getNbExtraNodalVars());
    for(CFuint iVar = 0; iVar < getNbExtraNodalVars();iVar++)
    {
      _extraNVarExists[iVar] = _extraNVarMap->exists((*(getExtraNodalVarNames()))[iVar]);
      if(_extraNVarExists[iVar])
      {
        std::string socketName = _extraNVarMap->find((*(getExtraNodalVarNames()))[iVar]).first;
	DataHandle<CFreal> extraVars = dynamicSockets->getSocketSource<CFreal>(socketName)->getDataHandle();
        extraVars.resize(nbNodes*(*(getExtraNodalVarStrides()))[iVar]);
        nbVarFound++;
      }
      else{
        CFLog(WARN, "Extra Nodal Variable : " << (*(getExtraNodalVarNames()))[iVar] << " present in CFmesh but not read...\n");
      }
    }

    if(nbVarFound != nbExtraVarsInCFcase)
    {
      CFLog(WARN, "Some extra nodal variables requested in the CFcase could be found...\n");
      cf_assert(nbVarFound == nbExtraVarsInCFcase);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::resizeStates(const CFuint nbStates)
{
  if (socket_states.getDataHandle().size() > 0) {
    const CFuint nbStates = socket_states.getDataHandle().size();
    for (CFuint i = 0; i < nbStates; ++i) {
      deletePtr(socket_states.getDataHandle()[i]);
    }
    // here resizing to 0 can lead to seg fault with some compilers
    socket_states.getDataHandle().resize(1);
  }
  
  cf_assert(socket_states.getDataHandle().size() == 0 ||
	    socket_states.getDataHandle().size() == 1);
  socket_states.getDataHandle().resize(nbStates);
  IndexList<State>::getList().reset();
  cf_assert(socket_states.getDataHandle().size() == nbStates);
  _StateLocalToGlobal.resize (nbStates);
  _StateOwnership.resize(nbStates);

  if(_storePastStates){
    //Resize the pastStates Datahandle
    std::string socketName = "pastStates";
    DataHandle<State*> pastStates = dynamicSockets->getSocketSink<State*>(socketName)->getDataHandle();
    pastStates.resize(nbStates);
  }

  if(_storeInterStates){
    //Resize the pastStates Datahandle
    std::string socketName = "interStates";
    DataHandle<State*> interStates = dynamicSockets->getSocketSink<State*>(socketName)->getDataHandle();
    interStates.resize(nbStates);
  }

  //Setup for the Extra Vars associated to the states
  const CFuint nbExtraVarsInMesh = getNbExtraStateVars();
  if((nbExtraVarsInMesh > 0) && (_extraSVarMap.isNotNull()))
  {
    const CFuint nbExtraVarsInCFcase = _extraSVarMap->size();
    cf_assert(nbExtraVarsInCFcase <= nbExtraVarsInMesh);

    ///resizing the state extra vars datahandles
    CFuint nbVarFound = 0;
    _extraSVarExists.resize(getNbExtraStateVars());
    for(CFuint iVar = 0; iVar < getNbExtraStateVars();iVar++)
    {
      _extraSVarExists[iVar] = _extraSVarMap->exists((*(getExtraStateVarNames()))[iVar]);
      if(_extraSVarExists[iVar])
      {
	std::string socketName = _extraSVarMap->find((*(getExtraStateVarNames()))[iVar]).first;
	DataHandle<CFreal> extraVars = dynamicSockets->getSocketSource<CFreal>(socketName)->getDataHandle();
	extraVars.resize(nbStates*(*(getExtraStateVarStrides()))[iVar]);
        nbVarFound++;
      }
      else {
        CFLog(WARN, "Extra State Variable : " << (*(getExtraStateVarNames()))[iVar] << " present in CFmesh but not read...\n");
      }
    }

    if(nbVarFound != nbExtraVarsInCFcase)
    {
      CFLog(WARN, "Some extra state variables requested in the CFcase could be found...\n");
      cf_assert(nbVarFound == nbExtraVarsInCFcase);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::resizeExtraVars()
{
  //Setup for the Extra Vars associated with the nodes
  const CFuint nbExtraVarsInMesh = getNbExtraVars();
  if((nbExtraVarsInMesh > 0) && (_extraVarMap.isNotNull()))
  {
    const CFuint nbExtraVarsInCFcase = _extraVarMap->size();
    cf_assert(nbExtraVarsInCFcase <= nbExtraVarsInMesh);
    ///resizing the extra vars datahandles
    CFuint nbVarFound = 0;
    _extraVarExists.resize(getNbExtraVars());
    for(CFuint iVar = 0; iVar < getNbExtraVars();iVar++)
    {
      _extraVarExists[iVar] = _extraVarMap->exists((*(getExtraVarNames()))[iVar]);
      if(_extraVarExists[iVar])
      {
        std::string socketName = _extraVarMap->find((*(getExtraVarNames()))[iVar]).first;

        DataHandle<CFreal> extraVars = dynamicSockets->getSocketSource<CFreal>(socketName)->getDataHandle();
        extraVars.resize((*(getExtraVarStrides()))[iVar]);
        nbVarFound++;
      }
      else{
        CFLog(WARN, "Extra Variable : " << (*(getExtraVarNames()))[iVar] << " present in CFmesh but not read...\n");
      }
    }
    
    if(nbVarFound != nbExtraVarsInCFcase)
    {
      CFLog(WARN, "Some extra variables requested in the CFcase could be found...\n");
      cf_assert(nbVarFound == nbExtraVarsInCFcase);
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::setPastNode(const CFuint nodeID, const RealVector& value)
{
  if(_storePastNodes){
    std::string socketName = "pastNodes";
    DataHandle<Node*> pastNodes = dynamicSockets->getSocketSink<Node*>(socketName)->getDataHandle();

    pastNodes[nodeID] = new Node();
    pastNodes[nodeID]->setIsOnMesh(false);
    *(pastNodes[nodeID]) = value;
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::setPastState(const CFuint stateID, const RealVector& value)
{
  if(_storePastStates){
    std::string socketName = "pastStates";
    DataHandle<State*> pastStates = dynamicSockets->getSocketSink<State*>(socketName)->getDataHandle();

    pastStates[stateID] = new State();
    *(pastStates[stateID]) = value;
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::setInterNode(const CFuint nodeID, const RealVector& value)
{
  if(_storeInterNodes){
    std::string socketName = "interNodes";
    DataHandle<Node*> interNodes = dynamicSockets->getSocketSink<Node*>(socketName)->getDataHandle();

    interNodes[nodeID] = new Node();
    interNodes[nodeID]->setIsOnMesh(false);
    *(interNodes[nodeID]) = value;
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::setInterState(const CFuint stateID, const RealVector& value)
{
  if(_storeInterStates){
    std::string socketName = "interStates";
    DataHandle<State*> interStates = dynamicSockets->getSocketSink<State*>(socketName)->getDataHandle();

    interStates[stateID] = new State();
    *(interStates[stateID]) = value;
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::setNodalExtraVar(const CFuint nodeID, const RealVector& value)
{
  cf_assert(value.size() == _totalNbNodalExtraVars);
  SafePtr<vector<CFuint> > strides = getExtraNodalVarStrides();
  CFuint count = 0;

  for(CFuint iVar = 0; iVar < _nbExtraNodalVars; ++iVar) {
    if(_extraNVarExists[iVar]) {
      DataHandle<CFreal> extraVars = _nodalExtraDataSourceSockets[iVar]->getDataHandle();
      const CFuint currStride = (*strides)[iVar];
      const CFuint start = nodeID*currStride;
      for (CFuint i = 0; i < currStride; ++i, ++count) {
	extraVars[start+i] = value[count];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::setStateExtraVar(const CFuint stateID, const RealVector& value)
{
  cf_assert(value.size() == _totalNbStateExtraVars);
  SafePtr<vector<CFuint> > strides = getExtraStateVarStrides();
  CFuint count = 0;

  for(CFuint iVar = 0; iVar < _nbExtraStateVars; ++iVar) {
    if(_extraSVarExists[iVar]) {
      DataHandle<CFreal> extraVars = _stateExtraDataSourceSockets[iVar]->getDataHandle();
      const CFuint currStride = (*strides)[iVar];
      const CFuint start = stateID*currStride;
      for (CFuint i = 0; i < currStride; ++i, ++count) {
        extraVars[start+i] = value[count];
      }
    }
  }
}
  
//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::setExtraVar(const RealVector& value)
{
  cf_assert(value.size() == _totalNbExtraVars);
  SafePtr<vector<CFuint> > strides = getExtraVarStrides();
  CFuint count = 0;
  
  for(CFuint iVar = 0; iVar < _nbExtraVars; ++iVar) {
    if(_extraVarExists[iVar]) {
      DataHandle<CFreal> extraVars = _extraDataSourceSockets[iVar]->getDataHandle();
      const CFuint currStride = (*strides)[iVar];
      const CFuint start = 0;
      for (CFuint i = 0; i < currStride; ++i, ++count) {
        extraVars[start+i] = value[count];
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::setNode(const CFuint nodeID, const RealVector& value)
{
  createNode(nodeID,value);
}

//////////////////////////////////////////////////////////////////////////////

Node* CFmeshReaderSource::createNode(const CFuint nodeID, const RealVector& value)
{
  const bool ownedByMesh = true;
  cf_assert(nodeID < socket_nodes.getDataHandle().size());

  Node* nodePtr;
  if (value.size() == PhysicalModelStack::getActive()->getDim()) {
    nodePtr = new Node(value,ownedByMesh);
  }
  else {
    nodePtr = new Node(ownedByMesh);
  }

  socket_nodes.getDataHandle()[nodeID] = nodePtr;
  IndexList<Node>::getList().createID(nodePtr);
  return nodePtr;
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::setNode(const CFuint nodeID, CFreal * Mem, const RealVector& D, bool IsUpdatable)
{
  createNode(nodeID,Mem,D,IsUpdatable);
}

//////////////////////////////////////////////////////////////////////////////

Node* CFmeshReaderSource::createNode(const CFuint nodeID, CFreal * Mem, const RealVector& D, bool IsUpdatable)
{
  const bool ownedByMesh = true;
  cf_assert(nodeID < socket_nodes.getDataHandle().size());

/// TODO: change this
  Node* nodePtr = new Node (Mem, ownedByMesh);

  nodePtr->setParUpdatable(IsUpdatable);

  socket_nodes.getDataHandle()[nodeID] = nodePtr;

  // This assumed the mesh would have the same state size ha
  //  cf_assert (statePtr->size()==D.size());

  if(nodePtr->size()==D.size()) {
    for (CFuint i=0; i<D.size(); ++i) {
      (*nodePtr)[i] = D[i];
    }
  }

  IndexList<Node>::getList().createID(nodePtr);
  return nodePtr;
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::setState(const CFuint stateID, const RealVector& value)
{
  createState(stateID,value);
}

//////////////////////////////////////////////////////////////////////////////

State* CFmeshReaderSource::createState(const CFuint stateID, const RealVector& value)
{
  cf_assert(stateID < socket_states.getDataHandle().size());
  
  const bool Ghost = false;

  State* statePtr;
  if (value.size() == PhysicalModelStack::getActive()->getNbEq()) {
    statePtr = new State(value);
  }
  else {
    statePtr = new State();
  }
  
  statePtr->setGhost(Ghost);
  
  socket_states.getDataHandle()[stateID] = statePtr;
  IndexList<State>::getList().createID(statePtr);
  return statePtr;
}
    
//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::setState(const CFuint stateID, CFreal * Mem, const RealVector& D, bool IsUpdatable)
{
  createState(stateID,Mem,D,IsUpdatable);
}

//////////////////////////////////////////////////////////////////////////////

State* CFmeshReaderSource::createState(const CFuint stateID, CFreal * Mem, const RealVector& D, bool IsUpdatable)
{
  cf_assert(stateID < socket_states.getDataHandle().size());

  const bool Ghost = false;

  /// TODO: change this
  State* statePtr = new State (Mem);

  statePtr->setGhost(Ghost);

  statePtr->setParUpdatable(IsUpdatable);
  
  socket_states.getDataHandle()[stateID] = statePtr;
  
  CFLog(DEBUG_MAX, stateID << "CFmeshReaderSource::createState() => socket_states.size() " 
	<< socket_states.getDataHandle().size() << " \n"); 
  
  // This assumed the mesh would have the same state size ha
  //  cf_assert (statePtr->size()==D.size());

  if(statePtr->size()==D.size())
  {
    for (CFuint i=0; i<D.size(); ++i)
    {
      (*statePtr)[i] = D[i];
    }
  }

  IndexList<State>::getList().createID(statePtr);
  CFLog(DEBUG_MAX, "CFmeshReaderSource::createState() => IndexList<State>::getList().size() "
	<< IndexList<State>::getList().size() << " \n");
  
  return statePtr;
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource:: getElementData (const CFuint EleNum, PairNodeState & Out) const
{
    cf_assert (_elementNode->nbRows() > 0);
    cf_assert (_elementState->nbRows() > 0);

    const unsigned int NC = _elementNode->nbCols(EleNum);
    const unsigned int SC = _elementState->nbCols(EleNum);

    cf_assert (NC > 0);
    cf_assert (SC > 0);

    PairNodeState PNS;

    Out.first.resize (NC);
    Out.second.resize (SC);

    for (unsigned int i=0; i<Out.first.size(); ++i) {
      Out.first[i] = (*_elementNode)(EleNum,i);
    }

    for (unsigned int i=0; i<Out.second.size(); ++i) {
      Out.second[i] = (*_elementState)(EleNum,i);
    }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::prepareNodalExtraVars()
{
  const vector<CFuint>& strides = *getExtraNodalVarStrides();
  _totalNbNodalExtraVars = std::accumulate(strides.begin(), strides.end(), 0);

  // here the order of extra vars is the order imposed by the input mesh file
  _nodalExtraDataSourceSockets.resize(_nbExtraNodalVars);
  for(CFuint iVar = 0; iVar < _nbExtraNodalVars;iVar++) {
    if(_extraNVarExists[iVar]) {
      std::string socketName = _extraNVarMap->find((*(getExtraNodalVarNames()))[iVar]).first;
      _nodalExtraDataSourceSockets[iVar] = dynamicSockets->getSocketSource<CFreal>(socketName);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::prepareStateExtraVars()
{
  const vector<CFuint>& strides = *getExtraStateVarStrides();
  _totalNbStateExtraVars = std::accumulate(strides.begin(), strides.end(), 0);

  // here the order of extra vars is the order imposed by the input mesh file
  _stateExtraDataSourceSockets.resize(_nbExtraStateVars);
  for(CFuint iVar = 0; iVar < _nbExtraStateVars;iVar++) {
    if(_extraSVarExists[iVar]) {
      std::string socketName = _extraSVarMap->find((*(getExtraStateVarNames()))[iVar]).first;
      _stateExtraDataSourceSockets[iVar] = dynamicSockets->getSocketSource<CFreal>(socketName);
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderSource::prepareExtraVars()
{
  const vector<CFuint>& strides = *getExtraVarStrides();
  _totalNbExtraVars = std::accumulate(strides.begin(), strides.end(), 0);
  
  // here the order of extra vars is the order imposed by the input mesh file
  _extraDataSourceSockets.resize(_nbExtraVars);
  for(CFuint iVar = 0; iVar < _nbExtraVars;iVar++) {
    if(_extraVarExists[iVar]) {
      std::string socketName = _extraVarMap->find((*(getExtraVarNames()))[iVar]).first;
      _extraDataSourceSockets[iVar] = dynamicSockets->getSocketSource<CFreal>(socketName);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
