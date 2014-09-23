#include "FluctSplit/STU_CRD_SplitStrategy.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/SpaceTime_Splitter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluctSplit {

/////////////////////////////////////////////////////////////////////////////
MethodStrategyProvider<STU_CRD_SplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitSpaceTimeModule>
spaceTimeucrdFluctSplitStrategyProvider("STU_CRD");

//////////////////////////////////////////////////////////////////////////////

STU_CRD_SplitStrategy::STU_CRD_SplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  socket_pastNormals("pastNormals",false),
  socket_pastStates("pastStates"),
  socket_pastCellVolume("pastCellVolume",false),
  socket_cellSpeed("cellSpeed",false),
  m_splitter(CFNULL),
  _pastStates(0),
  tPastStates(0),
  _null(),
  m_flux(),
  m_contourIntegrator(CFNULL),
  m_unitFaceNormals(),
  m_statesBkp(0),
  m_phi()
{
}

//////////////////////////////////////////////////////////////////////////////

STU_CRD_SplitStrategy::~STU_CRD_SplitStrategy()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void STU_CRD_SplitStrategy::computeFluctuation(vector<RealVector>& residual)
{
  ///@todo Nv: There is no moving mesh implemented
  
  DistributionData& ddata = getMethodData().getDistributionData();
  m_splitter->setDT(SubSystemStatusStack::getActive()->getDT());
  DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  // Since there is no moving mesh here, we set the past cell volume to the same as
  // the non moving one
  CellVolume = volumes[ddata.cellID];
  m_splitter->setCellVolume(volumes[ddata.cellID]);
  m_splitter->setPastCellVolume(volumes[ddata.cellID]);
  m_splitter->setCellSpeed(_null);
  setCurrentCell();
  
  ///@todo NV: It is possible to do as in SpaceTime_RD_SplitStrategy :
    
  ///      Then it is stored and we do not need to recompute it
  if (SubSystemStatusStack::getActive()->isFirstStep()){
    DataHandle<State*> pastStatesStorage = socket_pastStates.getDataHandle();
    
    // get the paststates in this cell
    for (CFuint i = 0; i < _nbStatesInCell; ++i) {
      const CFuint stateID = (*ddata.states)[i]->getLocalID();
      State & temp =  *_pastStates[i];
      //We need the local ID when we do Linearized Euler to access the meanflow
      temp.clone(*(*ddata.states)[i]);
      *_pastStates[i] = *pastStatesStorage[stateID];
    }
    
    //We point the past_residual from the DistributeData to the past residual of cellID
    //ddata.past_residuals = &_pastResiduals[ddata.cellID];
    
    // linearize the states in the cell
    // includes the transformation from update to linearization
    // variables to evaluate the jacobians in the average state
    // and then do the transformation from update to consistent variables
    tPastStates = computeConsistentStates(&_pastStates);
    // m_splitter->setConsStates(_pastStates);
    computeSTFluxIntegral(true);
    for (CFuint i = 0; i < nbEqs; ++i) {
      m_pastResiduals[ddata.cellID*nbEqs+i] = m_flux[i];
    }
  }
  else{
    // First load the past residual
    for (CFuint i = 0; i < nbEqs; ++i) {
      m_flux[i] = m_pastResiduals[ddata.cellID*nbEqs+i];
    }
  }
  
  //We do the same for the present states
  ddata.tStates = computeConsistentStates(ddata.states);
  // set the conservative states !! Not necessary here because we use CRD
  // m_splitter->setConsStates(*ddata.states);
  
  computeSTFluxIntegral(false);
    
  ddata.phi = m_flux;
  m_splitter->computeK(*ddata.states, normals[ddata.cellID]);
  
  ddata.phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux);  
  
  // Distribute the residual to the states of the present layer
  m_splitter->distribute(residual); 
}
    
//////////////////////////////////////////////////////////////////////////////

void STU_CRD_SplitStrategy::unsetup()
{
  // deallocating data for the temporary local residual
  for (CFuint i = 0; i < _pastStates.size(); ++i) {
    deletePtr(_pastStates[i]);
  }

  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    deletePtr(m_qdExtraVars[i]);
  }

  FluctuationSplitStrategy::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void STU_CRD_SplitStrategy::setup()
{
  CFAUTOTRACE;

  FluctuationSplitStrategy::setup();

  _null.resize(PhysicalModelStack::getActive()->getDim());
  _null = 0.;

  _dim = PhysicalModelStack::getActive()->getDim();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // get the splitter
  m_splitter = getMethodData().getSplitter().d_castTo<SpaceTime_Splitter>();
  m_splitter->setUpdateCoeff(updateCoeff);
  
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  // We assume that we have only one type of element
  // Then maxNbStatesInCell is also the number of states per cell
  _nbStatesInCell = maxNbStatesInCell;
  
  _pastStates.resize(maxNbStatesInCell);
  // Resizing pastStates
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _pastStates[i] = new State();
   }

  // Resizing m_flux_past, m_flux, m_phi
  nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_flux.resize(nbEqs);
  m_phi.resize(nbEqs);
 (getMethodData().getDistributionData().phi_time).resize(nbEqs);
 // AL: check if this is effective
  getMethodData().getUpdateVar()->setExtraData(true);

  // now complement it
  m_unitFaceNormals.resize(MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell());
  for (CFuint  i = 0; i < m_unitFaceNormals.size(); ++i) {
    m_unitFaceNormals[i].resize(PhysicalModelStack::getActive()->getDim());
  }

  // set up the contour integrator
  m_contourIntegrator = getMethodData().getContourIntegrator();
  m_contourIntegrator->setFaceNormals(&m_unitFaceNormals);

  const CFuint maxNbQPts = m_contourIntegrator->getMaxIntegratorPattern().totalNbPts();
  
  // physical data evaluated in the quadrature points
  m_pdata.resize(maxNbQPts);
  for (CFuint  i = 0; i < maxNbQPts; ++i) {
    PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
      resizePhysicalData(m_pdata[i]);
  }  
  
  m_contourIntegrator->getValues(m_qdstates);

  const CFuint extra_var_size = getMethodData().getUpdateVar()->getExtraPhysicalVarsSize();
  m_qdExtraVars.resize(m_qdstates.size());
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    m_qdExtraVars[i] = new RealVector(extra_var_size);
  }
  SafePtr<TopologicalRegionSet> innerCells = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbGeoEnts = innerCells->getLocalNbGeoEnts();
  m_nbQPointsInCell.resize(nbGeoEnts);
 
  // back up cell state pointers
  m_statesBkp.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());

  _stdTrsGeoBuilder.getDataGE().trs = innerCells;
 
  // loop over all the cells and set the number of quadrature
  // points in the m_nbQPointsInCell vector
  for (CFuint iGeoEnt = 0; iGeoEnt < nbGeoEnts; ++iGeoEnt) {
    
    // build the GeometricEntity
    _stdTrsGeoBuilder.getDataGE().idx = iGeoEnt;
    GeometricEntity& cell = *_stdTrsGeoBuilder.buildGE();
    // CFout << "m_nbQPointsInCell["<<iGeoEnt<<"]: "<< m_nbQPointsInCell[iGeoEnt]<<"\n";
    // CFout << "cell["<<iGeoEnt<<"].nbNodes: "<< cell.nbNodes()<<"\n";
    m_contourIntegrator->setNbSolQuadraturePoints(&cell, m_nbQPointsInCell[iGeoEnt]);
    
    // release the GeometricEntity
    _stdTrsGeoBuilder.releaseGE();
  }
  
  getMethodData().getDistributionData().isHO = false;
  m_pastResiduals.resize(nbGeoEnts*nbEqs);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
STU_CRD_SplitStrategy::needsSockets()
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

void STU_CRD_SplitStrategy::setCurrentCell()
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

void STU_CRD_SplitStrategy::computeSTFluxIntegral(bool isPast)
{ 
  const CFreal Area = CellVolume/(_dim+1.);
  (isPast) ? pastFlux(Area) : presentFlux(Area);
  
  // set the unit face normals
  DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();
  DistributionData& ddata = getMethodData().getDistributionData();
  const InwardNormalsData& normalsData = *normals[ddata.cellID];
  const CFuint nbFaces = normalsData.nbFaces();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    for (CFuint iDim = 0; iDim < _dim; ++iDim) {
      m_unitFaceNormals[iFace][iDim] = normalsData.getFaceUnitNormComp(iFace,iDim);
    }
  }
  
  // cell states substituted by the corresponding linearization
  // variables in which the interpolation will be performed
  // before the contour integration
  const CFuint nbCellStates = ddata.states->size();
  for (CFuint iState = 0; iState < nbCellStates; ++iState) {
    ddata.cell->setState(iState, (*m_linearStates)[iState]);
  }
  
  // first do interpolation then compute the physical related data for
  // the interpolated states
  m_contourIntegrator->computeAllAtQuadraturePointsOnGeo(ddata.cell);
  
  // compute the physical data for each interpolated state
  const CFuint nbQ = m_nbQPointsInCell[getMethodData().getDistributionData().cellID];
  computeStatesData(nbQ, getMethodData().getLinearVar(), m_qdstates,  m_pdata, m_qdExtraVars);

  // compute the contour integration of the fluxes
  m_contourIntegrator->integrateSolutionFunctorOnGeo<ConvectiveVarSet::Flux>
    (ddata.cell, getMethodData().getLinearVar()->getFlux(), m_pdata, m_phi);
  
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  const CFreal timeStep = dt/2.0;
  m_phi *= timeStep;
  
  // add to the space residual m_phi the contribution of the source term in the past n-level or 
  // in the present n+1-level
  ddata.time = (isPast) ? SubSystemStatusStack::getActive()->getCurrentTime() :
    SubSystemStatusStack::getActive()->getCurrentTime() + dt;
  
  // here there is a bug: dt/2*phi_S is missing! it has to be put inside the source term !!
  const CFuint nbST = getMethodData().getSourceTermSplitter()->size();
  for (CFuint i = 0; i < nbST; ++i) {
    ddata.sourceTermID = i;
    // computes the source term and places it in ddata.phiS
    getMethodData().getSourceTermSplitter(i)->computeSourceTerm(*normals[ddata.cellID]);
    m_phi -= ddata.phiS;
  }
  
  m_flux += m_phi;

  for (CFuint iState = 0; iState < nbCellStates; ++iState) {
    ddata.cell->setState(iState, m_statesBkp[iState]);
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit

} // namespace COOLFluiD
