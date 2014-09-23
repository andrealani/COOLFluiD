#include "FluctSplit/SpaceTimeRD_2LayerSplitStrategy.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/SpaceTime_Splitter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<SpaceTimeRD_2LayerSplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitSpaceTimeModule>
spaceTime2LayerStrategyProvider("ST2_RDS");

//////////////////////////////////////////////////////////////////////////////

SpaceTimeRD_2LayerSplitStrategy::SpaceTimeRD_2LayerSplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  socket_pastNormals("pastNormals",false),
  socket_interNormals("interNormals",false),
  socket_pastStates("pastStates"),
  socket_interStates("interStates"),
  socket_pastCellVolume("pastCellVolume",false),
  socket_cellSpeed("cellSpeed",false),
  socket_interUpdateCoeff("interUpdateCoeff"),
  m_splitter(CFNULL),
  _pastStates(0),
  _interStates(0),
  _null(),
  _pastResiduals(0),
  _interResiduals(0),
  _interVolume()
{
}

//////////////////////////////////////////////////////////////////////////////

SpaceTimeRD_2LayerSplitStrategy::~SpaceTimeRD_2LayerSplitStrategy()
{
  // deallocating data for the temporary local residual
  for (CFuint i = 0; i < _pastStates.size(); ++i) {
    deletePtr(_pastStates[i]);
  }

  // allocating data for the temporary local residual
  for (CFuint i = 0; i < _pastStates.size(); ++i) {
    deletePtr(_interStates[i]);
  }
 }

//////////////////////////////////////////////////////////////////////////////

void SpaceTimeRD_2LayerSplitStrategy::computeFluctuation(vector<RealVector>& residual)
{
  (getMethodData().getDistributionData().isPerturb) ?
    doComputeFluct(residual) : doComputeFluctAndUpdateCoeff(residual);  
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTimeRD_2LayerSplitStrategy::doComputeFluct(vector<RealVector>& residual)
{

  ///@todo NV:I think that this should still be modified to do blending in Space-Time 2 Layer
  cf_assert(getMethodData().getDistributionData().isPerturb);
  
  DistributionData& ddata = getMethodData().getDistributionData();
  const CFuint cellID = ddata.cellID;
  const CFuint nbStatesInCell = ddata.states->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle< CFreal> volumes = socket_volumes.getDataHandle();
  
  if (SubSystemStatusStack::getActive()->isMovingMesh() != false){
    DataHandle< CFreal> pastCellVolume = socket_pastCellVolume.getDataHandle();
    DataHandle< RealVector> cellSpeed = socket_cellSpeed.getDataHandle();
    
    m_splitter->setPastCellVolume(pastCellVolume[cellID]);
    
    // build the GeometricEntity linked to cellID
    _stdTrsGeoBuilder.getDataGE().trs = MeshDataStack::getActive()->getTrs("InnerCells");
    _stdTrsGeoBuilder.getDataGE().idx = cellID;

    GeometricEntity& cell = *_stdTrsGeoBuilder.buildGE();
    
    volumes[cellID] = cell.computeVolume();
    m_splitter->setCellVolume(volumes[cellID]);
    m_splitter->setCellSpeed(cellSpeed[cellID]);
    _stdTrsGeoBuilder.releaseGE();
  }
  else{
    m_splitter->setPastCellVolume(volumes[cellID]);
    m_splitter->setCellVolume(volumes[cellID]);
    m_splitter->setCellSpeed(_null);
  }
  
  m_splitter->setDT(SubSystemStatusStack::getActive()->getInnerDT(0));
  
  /// Compute residual from K1
  // If it is the first step, then first
  // compute the residuals due to the past
  if (SubSystemStatusStack::getActive()->isFirstStep()){
    DataHandle<State*> pastStatesStorage = socket_pastStates.getDataHandle();
    
    // get the paststates in this cell
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      const CFuint stateID = (*ddata.states)[i]->getLocalID();
      *_pastStates[i] = *pastStatesStorage[stateID];
    }
    
    // linearize the states in the cell
    // includes the transformation from update to linearization
    // variables to evaluate the jacobians in the average state
    // and then do the transformation from update to consistent variables
    vector<State*> *const tPastStates = computeConsistentStates(&_pastStates);
    
    // Compute the upwind parameters k in this cell at past states
    if (SubSystemStatusStack::getActive()->isMovingMesh()){
      DataHandle< InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();
      m_splitter->computeK(_pastStates, pastNormals[cellID]);
    }
    else{
      m_splitter->computeK(_pastStates, normals[cellID]);
    }
    
    // set the conservative past states
    m_splitter->setConsStates(_pastStates);
 
   //We point the past_residual from the DistributeData to the past residual of cellID
    getMethodData().getDistributionData().past_residuals = &_pastResiduals[ddata.cellID];
    // Compute the part of the residual from past states
    m_splitter->distributePast(*tPastStates);
  }

  /// Add the past contributions to the residual of intermediate layer
  for (CFuint iState=0; iState < nbStatesInCell; ++iState){
    for (CFuint j = 0; j < nbEqs; ++j) {
      residual[iState][j] = _pastResiduals[cellID][iState*nbEqs + j];
    }
  }
  
  DataHandle<State*> interStatesStorage = socket_interStates.getDataHandle();
  // K1: contribution from intermediate states
  // get the intermediate States in this cell
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    const CFuint stateID = (*ddata.states)[i]->getLocalID();
    *_interStates[i] = *interStatesStorage[stateID];
  }
  
  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  vector<State*> *const tInterStates = computeConsistentStates(&_interStates);
  
  // set the conservative states
  m_splitter->setConsStates(_interStates);
  
  // Compute the upwind parameters k in this cell at intermediate states
  if (SubSystemStatusStack::getActive()->isMovingMesh()){
    DataHandle< InwardNormalsData*> interNormals = socket_interNormals.getDataHandle();
    m_splitter->computeK(_interStates, interNormals[cellID]);
  }
  else{
    m_splitter->computeK(_interStates, normals[cellID]);
  }
  // compute the residual from cell K1
  m_splitter->distributeInterK1(*tInterStates, residual);
  
  // compute the residual from cell K2
  m_splitter->setDT(SubSystemStatusStack::getActive()->getInnerDT(1));
  // K2: contribution from intermediate states
  m_splitter->distributeInterK2(*tInterStates, residual);
  
  // Compute residual phi_n+1 (in newresidual) from K2
  ddata.tStates = computeConsistentStates(ddata.states);
  
  // set the conservative states
  m_splitter->setConsStates(*ddata.states);
  m_splitter->setInterStates(*tInterStates);
  
  // compute the residual and the upwind parameters k in this cell at current states
  m_splitter->computeK(*ddata.states, normals[cellID]);
  
  m_splitter->distribute(residual);
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTimeRD_2LayerSplitStrategy::doComputeFluctAndUpdateCoeff(vector<RealVector>& residual)
{
   ///@todo NV:I think that this should still be modified to do blending in Space-Time 2 Layer

  cf_assert(!getMethodData().getDistributionData().isPerturb);
  
  DistributionData& ddata = getMethodData().getDistributionData();
  const CFuint cellID = ddata.cellID;
  const CFuint nbStatesInCell = ddata.states->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle< CFreal> volumes = socket_volumes.getDataHandle();
    
  if (SubSystemStatusStack::getActive()->isMovingMesh() != false){
    DataHandle< CFreal> pastCellVolume = socket_pastCellVolume.getDataHandle();
    DataHandle< RealVector> cellSpeed = socket_cellSpeed.getDataHandle();
    
    m_splitter->setPastCellVolume(pastCellVolume[cellID]);
    
    // build the GeometricEntity linked to cellID
    _stdTrsGeoBuilder.getDataGE().trs = MeshDataStack::getActive()->getTrs("InnerCells");
    _stdTrsGeoBuilder.getDataGE().idx = cellID;

    GeometricEntity& cell = *_stdTrsGeoBuilder.buildGE();
    
    volumes[cellID] = cell.computeVolume();
    _interVolume = (pastCellVolume[cellID]*SubSystemStatusStack::getActive()->getInnerDTRatio(1)) + 
      (volumes[cellID]*SubSystemStatusStack::getActive()->getInnerDTRatio(0));
    m_splitter->setCellVolume(_interVolume);
    m_splitter->setCellSpeed(cellSpeed[cellID]);
    _stdTrsGeoBuilder.releaseGE();
  }
  else{
    m_splitter->setPastCellVolume(volumes[cellID]);
    m_splitter->setCellVolume(volumes[cellID]);
    m_splitter->setCellSpeed(_null);
  }
  
  m_splitter->setDT(SubSystemStatusStack::getActive()->getInnerDT(0));
  // Compute residual from K1
  // If it is the first step, then first
  // compute the residuals due to the past
  if (SubSystemStatusStack::getActive()->isFirstStep()){
    DataHandle<State*> pastStatesStorage = socket_pastStates.getDataHandle();
    
    // get the paststates in this cell
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      const CFuint stateID = (*ddata.states)[i]->getLocalID();
      *_pastStates[i] = *pastStatesStorage[stateID];
    }
    
    // linearize the states in the cell
    // includes the transformation from update to linearization
    // variables to evaluate the jacobians in the average state
    // and then do the transformation from update to consistent variables
    vector<State*> *const tPastStates = computeConsistentStates(&_pastStates);
    
    // back up the flag telling if you are perturbing and set it to true
    // so that the update coefficient is not computed
    const bool isPerturbBkp = getMethodData().getDistributionData().isPerturb;
    getMethodData().getDistributionData().isPerturb = true;
    
    // Compute the upwind parameters k in this cell at past states
    if (SubSystemStatusStack::getActive()->isMovingMesh()){
      DataHandle< InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();
      m_splitter->computeK(_pastStates, pastNormals[cellID]);
    }
    else{
      m_splitter->computeK(_pastStates, normals[cellID]);
    }
    
    // restore the backed up flag
    getMethodData().getDistributionData().isPerturb = isPerturbBkp;
    
    // set the conservative past states
    m_splitter->setConsStates(_pastStates);

    //We point the past_residual from the DistributeData to the past residual of cellID
      getMethodData().getDistributionData().past_residuals = &_pastResiduals[ddata.cellID];

    // Compute the part of the residual from past states
    m_splitter->distributePast(*tPastStates);
  }

  /// Add the past contributions to the residual of intermediate layer
  for (CFuint iState=0; iState < nbStatesInCell; ++iState){
    //  residual[iState].slice(0) = _pastResiduals[cellID].slice(iState*nbEqs);
    for (CFuint j = 0; j < nbEqs; ++j) {
      residual[iState][j] = _pastResiduals[cellID][iState*nbEqs + j];
    }
  }
  
  DataHandle<State*> interStatesStorage = socket_interStates.getDataHandle();
  // K1: contribution from intermediate states
  // get the intermediate States in this cell
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    const CFuint stateID = (*ddata.states)[i]->getLocalID();
    *_interStates[i] = *interStatesStorage[stateID];
  }
  
  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  vector<State*> *const tInterStates = computeConsistentStates(&_interStates);

  // set the conservative states
  m_splitter->setConsStates(_interStates);
   
  // back up the flag telling if you are perturbing and set it to true
  // so that the update coefficient is not computed
  const bool isPerturbBkp = getMethodData().getDistributionData().isPerturb;
  getMethodData().getDistributionData().isPerturb = true;
  
  // Compute the upwind parameters k in this cell at intermediate states
  if (SubSystemStatusStack::getActive()->isMovingMesh()){
    DataHandle< InwardNormalsData*> interNormals = socket_interNormals.getDataHandle();
    m_splitter->computeK(_interStates, interNormals[cellID]);
  }
  else{
    m_splitter->computeK(_interStates, normals[cellID]);
  }
  
  // restore the backed up flag
  getMethodData().getDistributionData().isPerturb = isPerturbBkp;
  
  // compute the residual from cell K1
  m_splitter->distributeInterK1(*tInterStates, residual);
     
  // Now we work on cell K2
  // Set the cell Volumes and Speed
  if (SubSystemStatusStack::getActive()->isMovingMesh() != false){
    DataHandle< CFreal> pastCellVolume = socket_pastCellVolume.getDataHandle();
    DataHandle< RealVector> cellSpeed = socket_cellSpeed.getDataHandle();
    
    // pastVolume is an average between past and current.
    _interVolume = (pastCellVolume[cellID]*SubSystemStatusStack::getActive()->getInnerDTRatio(1)) 
      + (volumes[cellID]*SubSystemStatusStack::getActive()->getInnerDTRatio(0));
    m_splitter->setPastCellVolume(_interVolume);
    m_splitter->setCellVolume(volumes[cellID]);
    m_splitter->setCellSpeed(cellSpeed[cellID]);
  }
  
  // compute the residual from cell K2
  // Set the timeStep to dt2
  m_splitter->setDT(SubSystemStatusStack::getActive()->getInnerDT(1));
  
  // K2: contribution from intermediate states
  m_splitter->distributeInterK2(*tInterStates, residual);
  
  // Compute the Update coefficient for intermediate states
  DataHandle< CFreal> interUpdateCoeff = socket_interUpdateCoeff.getDataHandle();
  m_splitter->setUpdateCoeff(interUpdateCoeff);
  
  /// @todo use _interStates and not states!! (for consistency, no difference)
  // here we use states and not interStates because ??????
  if (SubSystemStatusStack::getActive()->isMovingMesh()){
    DataHandle< InwardNormalsData*> interNormals = socket_interNormals.getDataHandle();
    m_splitter->computeK(*ddata.states, interNormals[cellID]);
  }
  else{
    m_splitter->computeK(*ddata.states, normals[cellID]);
  }
  
  // K2: contribution from current states
  ddata.tStates = computeConsistentStates(ddata.states);
  
  // set the conservative states
  m_splitter->setConsStates(*ddata.states);
  // set the intermediate states
  m_splitter->setInterConsStates(_interStates);
  m_splitter->setInterStates(*tInterStates);
  
  // compute the residual and the upwind parameters k in this cell at current states
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  m_splitter->setUpdateCoeff(updateCoeff);
  m_splitter->computeK(*ddata.states, normals[cellID]);
  m_splitter->distribute(residual);
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTimeRD_2LayerSplitStrategy::unsetup()
{

  for (CFuint i = 0; i < _pastStates.size(); ++i) {
    deletePtr(_pastStates[i]);
  }

  for (CFuint i = 0; i < _interStates.size(); ++i) {
    deletePtr(_interStates[i]);
  }

  FluctuationSplitStrategy::unsetup();

}

//////////////////////////////////////////////////////////////////////////////

void SpaceTimeRD_2LayerSplitStrategy::setup()
{
  FluctuationSplitStrategy::setup();

  _null.resize(PhysicalModelStack::getActive()->getDim());
  _null = 0.;

  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
  // get the splitter
  m_splitter = getMethodData().getSplitter().d_castTo<SpaceTime_Splitter>();
  m_splitter->setUpdateCoeff(updateCoeff);

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  _pastStates.resize(maxNbStatesInCell);
  // Resizing pastStates
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _pastStates[i] = new State();
  }

  _interStates.resize(maxNbStatesInCell);
  // Resizing pastStates
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _interStates[i] = new State();
  }

  // if the mesh is moving, then get the necessary information
  if (SubSystemStatusStack::getActive()->isMovingMesh()) {
    cf_assert(socket_cellSpeed.isConnected());
    cf_assert(socket_pastCellVolume.isConnected());
    cf_assert(socket_pastNormals.isConnected());
    cf_assert(socket_interNormals.isConnected());
  }

  Common::SafePtr<TopologicalRegionSet> innerCells = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = innerCells->getLocalNbGeoEnts();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  CFuint nbStatesInCell;

  /// Resizing the storage for the full past Residuals
  _pastResiduals.resize(nbCells);
  // Resize the RealVector
  _stdTrsGeoBuilder.getDataGE().trs = innerCells;

  for(CFuint iGeoEnt = 0; iGeoEnt < nbCells; ++iGeoEnt) {

    CFLogDebugMax("Cell " << iGeoEnt << "\n");

    // build the GeometricEntity
    _stdTrsGeoBuilder.getDataGE().idx = iGeoEnt;
    GeometricEntity& cell = *_stdTrsGeoBuilder.buildGE();

    nbStatesInCell = cell.nbStates();
    _pastResiduals[iGeoEnt].resize(nbStatesInCell*nbEqs);

    //release the GeometricEntity
    _stdTrsGeoBuilder.releaseGE();
  }

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
SpaceTimeRD_2LayerSplitStrategy::needsSockets()
{
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = FluctuationSplitStrategy::needsSockets();
   
   result.push_back(&socket_updateCoeff);
   result.push_back(&socket_pastNormals);
   result.push_back(&socket_interNormals);
   result.push_back(&socket_pastStates);
   result.push_back(&socket_interStates);
   result.push_back(&socket_pastCellVolume);
   result.push_back(&socket_cellSpeed);
   result.push_back(&socket_interUpdateCoeff);

   return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
