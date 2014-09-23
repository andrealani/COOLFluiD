#include "FluctSplit/STM_RDS_SplitStrategy.hh"

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

namespace COOLFluiD
{

		namespace FluctSplit
		{

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STM_RDS_SplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitSpaceTimeModule>
			  aSTM_RDS_SplitStrategyProvider("STM_RDS");

//////////////////////////////////////////////////////////////////////////////

STM_RDS_SplitStrategy::STM_RDS_SplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  socket_pastNormals("pastNormals",false),
  socket_pastStates("pastStates"),
  socket_pastCellVolume("pastCellVolume",false),
  socket_cellSpeed("cellSpeed",false),
  m_splitter(CFNULL),
  _pastStates(0),
  m_tmpvec(),
  _pastResiduals(0),
temp_residual(0)
{
}

//////////////////////////////////////////////////////////////////////////////

STM_RDS_SplitStrategy::~STM_RDS_SplitStrategy()
{
  // deallocating data for the temporary local residual
  for (CFuint i = 0; i < _pastStates.size(); ++i) {
    deletePtr(_pastStates[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void STM_RDS_SplitStrategy::computeFluctuation(vector<RealVector>& residual)
{

  DistributionData& ddata = getMethodData().getDistributionData();
  const CFuint nbStatesInCell = ddata.states->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_splitter->setDT(SubSystemStatusStack::getActive()->getDT());
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle< CFreal> volumes = socket_volumes.getDataHandle();
  setCurrentCell();
  if (SubSystemStatusStack::getActive()->isMovingMesh() != false)
  {
    DataHandle< RealVector> cellSpeed = socket_cellSpeed.getDataHandle();
    DataHandle< CFreal> pastCellVolume = socket_pastCellVolume.getDataHandle();

    m_splitter->setPastCellVolume(pastCellVolume[ddata.cellID]);
    // build the GeometricEntity linked to ddata.cellID
    _stdTrsGeoBuilder.getDataGE().trs = MeshDataStack::getActive()->getTrs("InnerCells");
    _stdTrsGeoBuilder.getDataGE().idx = ddata.cellID;

    GeometricEntity& cell = *_stdTrsGeoBuilder.buildGE();

    volumes[ddata.cellID] = cell.computeVolume();
    m_splitter->setCellVolume(volumes[ddata.cellID]);
    m_splitter->setCellSpeed(cellSpeed[ddata.cellID]);
    _stdTrsGeoBuilder.releaseGE();
  }
  else
  {
    m_splitter->setPastCellVolume(volumes[ddata.cellID]);
    m_splitter->setCellVolume(volumes[ddata.cellID]);
    m_splitter->setCellSpeed(m_tmpvec);
  }
  // If it is the first step, then first
  // compute the residuals due to the past
  if (SubSystemStatusStack::getActive()->isFirstStep()){
    DataHandle<State*> pastStatesStorage = socket_pastStates.getDataHandle();
    // get the paststates in this cell
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      const CFuint stateID = (*ddata.states)[i]->getLocalID();
       State & temp =  *_pastStates[i];
    //We need the local ID when we do Linearized Euler
    //to access the meanflow
    temp.clone(*(*ddata.states)[i]);
    *_pastStates[i] = *pastStatesStorage[stateID];

    }
    // linearize the states in the cell
    // includes the transformation from update to linearization
    // variables to evaluate the jacobians in the average state
    // and then do the transformation from update to consistent variables
    vector<State*> *const tPastStates = computeConsistentStates(&_pastStates);
    // back up the flag telling if you are perturbing and set it to true
    // so that the update coefficient is not computed
    bool backUpPerturb = ddata.isPerturb;
    ddata.isPerturb = true;
    // Compute the upwind parameters k in this cell at past states
    if (SubSystemStatusStack::getActive()->isMovingMesh())
    {
      DataHandle< InwardNormalsData*> pastNormals = socket_pastNormals.getDataHandle();
      m_splitter->computeK(_pastStates, pastNormals[ddata.cellID]);
    }
    else
    {
      m_splitter->computeK(_pastStates, normals[ddata.cellID]);
    }
    // restore the backed up flag
    ddata.isPerturb = backUpPerturb;
    //We point the past_residual from the DistributeData to the past residual of cellID
    ddata.past_residuals = &_pastResiduals[ddata.cellID];
//CF_DEBUG_OBJ((ddata.past_residuals_order1).size());
    //We point the past_residual of order1 from the DistributeData to the past residual of order 1 of cellID
    // If we are not unsing blending scheme these are vector of dimension 0
    ddata.past_residuals_order1 = &_pastResiduals_order1[ddata.cellID];
    // set the conservative past states
    m_splitter->setConsStates(_pastStates);
    ddata.time = SubSystemStatusStack::getActive()->getCurrentTime();
    const CFuint nbST = getMethodData().getSourceTermSplitter()->size();


    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
       for (CFuint iEq = 0; iEq <  nbEqs; ++iEq){
        temp_residual[iState][iEq]=0.0;
      }
    }

   if (nbST == 1) {
      ddata.sourceTermID = 0;
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        ddata.cell->setState(iState, (*tPastStates)[iState]);
      }
      getMethodData().getSourceTermSplitter(0)->computeSourceTerm(*normals[ddata.cellID]);
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        ddata.cell->setState(iState, m_statesBkp[iState]);
      }
      m_splitter->distributePast(*tPastStates);
      getMethodData().getSourceTermSplitter(0)->distribute(temp_residual);
   }

   else {
      m_splitter->distributePast(*tPastStates);
      for (CFuint i = 0; i < nbST; ++i) {
           ddata.sourceTermID = i;
            for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
               ddata.cell->setState(iState, _pastStates[iState]);
            }
            getMethodData().getSourceTermSplitter(i)->computeSourceTerm(*normals[ddata.cellID]);
            for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
              ddata.cell->setState(iState, m_statesBkp[iState]);
            }
            getMethodData().getSourceTermSplitter(i)->distribute(temp_residual);
        }
   }
  RealVector& past_residuals = *getMethodData().getDistributionData().past_residuals;
   for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
     for (CFuint j=0; j< nbEqs; ++j){
        past_residuals[(iState*nbEqs)+j] += temp_residual[iState][j];

     }
   }
}
  // Add the past contributions to the residual
  //for (CFuint iState=0; iState < nbStatesInCell; ++iState){
    //  residual[iState].slice(0) = _pastResiduals[ddata.cellID].slice(iState*nbEqs);
   // for (CFuint j = 0; j < nbEqs; ++j) {
     // (residual[iState])[j] = _pastResiduals[ddata.cellID][iState*nbEqs + j];
   // }
 // }
  if (!SubSystemStatusStack::getActive()->isFirstStep())
  {
    ddata.past_residuals = &_pastResiduals[ddata.cellID];
    ddata.past_residuals_order1 = &_pastResiduals_order1[ddata.cellID];
  }
  ddata.tStates = computeConsistentStates(ddata.states);
  // set the conservative states
  m_splitter->setConsStates(*ddata.states);
  // compute the residual and the upwind parameters k in this cell at current states
  m_splitter->computeK(*ddata.states,normals[ddata.cellID]);
  // Distribute the residual to the states of the present layer
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  ddata.time = SubSystemStatusStack::getActive()->getCurrentTime()+dt;
        const CFuint nbST = getMethodData().getSourceTermSplitter()->size();
        if (nbST == 1) {
          ddata.sourceTermID = 0;
          getMethodData().getSourceTermSplitter(0)->computeSourceTerm(*normals[ddata.cellID]);
          m_splitter->distribute(residual);
          getMethodData().getSourceTermSplitter(0)->distribute(residual);
        }
        else {
          m_splitter->distribute(residual);

          for (CFuint i = 0; i < nbST; ++i) {
            ddata.sourceTermID = i;
            getMethodData().getSourceTermSplitter(i)->computeSourceTerm(*normals[ddata.cellID]);
            getMethodData().getSourceTermSplitter(i)->distribute(residual);
          }
        }


			}


//////////////////////////////////////////////////////////////////////////////

void STM_RDS_SplitStrategy::unsetup()
{
  CFAUTOTRACE;

  for (CFuint i = 0; i < _pastStates.size(); ++i) {
    deletePtr(_pastStates[i]);
  }

  FluctuationSplitStrategy::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void STM_RDS_SplitStrategy::setup()
{
  CFAUTOTRACE;

  FluctuationSplitStrategy::setup();

  m_tmpvec.resize(PhysicalModelStack::getActive()->getDim());
  m_tmpvec = 0.;

  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // get the splitter
  m_splitter = getMethodData().getSplitter().d_castTo<SpaceTime_Splitter>();
  m_splitter->setUpdateCoeff(updateCoeff);

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

// back up cell state pointers
  m_statesBkp.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());
  _pastStates.resize(maxNbStatesInCell);
  // Resizing pastStates
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _pastStates[i] = new State();
  }

  // if the mesh is moving, then get the necessary information
  if (SubSystemStatusStack::getActive()->isMovingMesh() != false){
    cf_assert(socket_cellSpeed.isConnected());
    cf_assert(socket_pastCellVolume.isConnected());
    cf_assert(socket_pastNormals.isConnected());
  }

  Common::SafePtr<TopologicalRegionSet> innerCells = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = innerCells->getLocalNbGeoEnts();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  CFuint nbStatesInCell;

 // Resizing the storage for the full past Residuals of scheme of order 1 (if we use ST blending)
  std::string name = m_splitter->getName();
  // Resizing the storage for the full past Residuals
  _pastResiduals_order1.resize(nbCells);
  // Resize the RealVector
  _stdTrsGeoBuilder.getDataGE().trs = innerCells;

  for(CFuint iGeoEnt = 0; iGeoEnt < nbCells; ++iGeoEnt) {

    CFLogDebugMax("Cell " << iGeoEnt << "\n");

    // build the GeometricEntity
    _stdTrsGeoBuilder.getDataGE().idx = iGeoEnt;
    GeometricEntity& cell = *_stdTrsGeoBuilder.buildGE();

    nbStatesInCell = cell.nbStates();
    _pastResiduals_order1[iGeoEnt].resize(nbStatesInCell*nbEqs);

    //release the GeometricEntity
    _stdTrsGeoBuilder.releaseGE();
   }
 
  // Resizing the storage for the full past Residuals
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
temp_residual.resize(maxNbStatesInCell);

  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState)
    temp_residual[iState].resize(nbEqs);

getMethodData().getDistributionData().isHO = false;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
STM_RDS_SplitStrategy::needsSockets()
{
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = FluctuationSplitStrategy::needsSockets();

   result.push_back(&socket_updateCoeff);
   result.push_back(&socket_pastNormals);
   result.push_back(&socket_pastStates);
   result.push_back(&socket_pastCellVolume);
   result.push_back(&socket_cellSpeed);

   return result;
}
//////////////////////////////////////////////////////////////////////////////

void STM_RDS_SplitStrategy::setCurrentCell()
{
  DistributionData& ddata = getMethodData().getDistributionData();

  // back up the update states and sets them in the linearizer
  // AL: don't touch this! update states have to be backed up
  // because linearizing states are temporarily stored
  // in the cell and used to contour integrate
  const CFuint nbCellStates = ddata.states->size();

  for (CFuint iState = 0; iState < nbCellStates; ++iState) {

    m_statesBkp[iState] = (*ddata.states)[iState];
  }
  getMethodData().getLinearizer()->setUpdateStates(&m_statesBkp);

  ddata.tStates = computeConsistentStates(ddata.states);
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD
