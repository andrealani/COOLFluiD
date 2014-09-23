#include "FluctSplit/STM_CRD_SplitStrategy.hh"

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

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STM_CRD_SplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitSpaceTimeModule>
spaceTimemcrdFluctSplitStrategyProvider("STM_CRD");

//////////////////////////////////////////////////////////////////////////////

STM_CRD_SplitStrategy::STM_CRD_SplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  socket_pastNormals("pastNormals",false),
  socket_pastStates("pastStates"),
  socket_pastCellVolume("pastCellVolume",false),
  socket_cellSpeed("cellSpeed",false),
  m_splitter(CFNULL),
  _pastStates(0),
  _null(),
  m_flux_past_time(),
  m_flux_time(),
  m_flux_past_space(),
  m_flux_space(),
  CellVolume(),
  m_qdstates(0),
  m_qdExtraVars(0),
  m_contourIntegrator(CFNULL),
  m_nbQPointsInCell(),
  m_unitFaceNormals(),
  m_statesBkp(0),
  m_phi(),
  _pastResiduals(0),
  _pastResiduals_order1(0),
  temp_residual(0)
{
}

//////////////////////////////////////////////////////////////////////////////

STM_CRD_SplitStrategy::~STM_CRD_SplitStrategy()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void STM_CRD_SplitStrategy::computeFluctuation(vector<RealVector>& residual)
{
  ///@todo Nv: There is no moving mesh implemented

  DistributionData& ddata = getMethodData().getDistributionData();
  const CFuint nbStatesInCell = ddata.states->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
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

  // If it is the first step, then first
  // compute the residuals due to the past
  if (SubSystemStatusStack::getActive()->isFirstStep()){
    DataHandle<State*> pastStatesStorage = socket_pastStates.getDataHandle();

    // get the paststates in this cell
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
    State & temp =  *_pastStates[i];
    //We need the local ID when we do Linearized Euler
    //to access the meanflow
  temp.clone(*(*ddata.states)[i]);
     const CFuint stateID = (*ddata.states)[i]->getLocalID();
      *_pastStates[i] = *pastStatesStorage[stateID];

    }

// back up the flag telling if you are perturbing and set it to true
    // so that the update coefficient is not computed
    bool backUpPerturb = ddata.isPerturb;
    ddata.isPerturb = true;


    //We point the past_residual from the DistributeData to the past residual of cellID
    ddata.past_residuals = &_pastResiduals[ddata.cellID];

    //We point the past_residual of order1 from the DistributeData to the past residual of order 1 of cellID
    // If we are not unsing blending scheme these are vector of dimension 0
    ddata.past_residuals_order1 = &_pastResiduals_order1[ddata.cellID];

    m_splitter->setConsStates(_pastStates);

    //We do the same for the present states
    ddata.tStates = computeConsistentStates(&_pastStates);

    computeSTFluxIntegral_past();

    getMethodData().getDistributionData().phi = m_flux_space;
   // CF_DEBUG_OBJ(_pastStates[0]->getLocalID());
    m_splitter->computeK(_pastStates, normals[ddata.cellID]);

// restore the backed up flag
    ddata.isPerturb = backUpPerturb;


//We point the past_residual from the DistributeData to the past residual of cellID
    ddata.past_residuals = &_pastResiduals[ddata.cellID];

    const CFreal dt = SubSystemStatusStack::getActive()->getDT();
    ddata.time = SubSystemStatusStack::getActive()->getCurrentTime()-dt;

    const CFuint nbST = getMethodData().getSourceTermSplitter()->size();



    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
       for (CFuint iEq = 0; iEq <  nbEqs; ++iEq){
        temp_residual[iState][iEq]=0.0;
      }
    }


   if (nbST == 1) {
      ddata.sourceTermID = 0;

      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        ddata.cell->setState(iState, _pastStates[iState]);
      }

      getMethodData().getSourceTermSplitter(0)->computeSourceTerm(*normals[ddata.cellID]);

      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        ddata.cell->setState(iState, m_statesBkp[iState]);
      }

      m_splitter->distributePast(_pastStates);

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
 if (!SubSystemStatusStack::getActive()->isFirstStep())
  {
    ddata.past_residuals = &_pastResiduals[ddata.cellID];
    ddata.past_residuals_order1 = &_pastResiduals_order1[ddata.cellID];
  }

  computeSTFluxIntegral();
ddata.tStates = computeConsistentStates(ddata.states);
   m_splitter->computeK(*ddata.states, normals[ddata.cellID]);
   ddata.phi = m_flux_space;
    // Distribute the residual to the states of the present layer
ddata.time = SubSystemStatusStack::getActive()->getCurrentTime();

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

void STM_CRD_SplitStrategy::computeSTFluxIntegral(){

  DistributionData& ddata = getMethodData().getDistributionData();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();
  const InwardNormalsData& normalsData = *normals[ddata.cellID];

  // set the face normals
  const CFuint nbFaces = normalsData.nbFaces();

 // Compute the contour integral of the flux
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      m_unitFaceNormals[iFace][iDim] =
        normalsData.getFaceUnitNormComp(iFace,iDim);
    }
  }

  // first do interpolation then compute the physical related data for
  // the interpolated states
  m_contourIntegrator->computeAllAtQuadraturePointsOnGeo(ddata.cell);

  // compute the physical data for each interpolated state
  const CFuint nbQ = m_nbQPointsInCell[getMethodData().getDistributionData().cellID];
  computeStatesData(nbQ, getMethodData().getLinearVar(), m_qdstates,  m_pdata, m_qdExtraVars);
    
  // compute the contour integration of the fluxes
  m_contourIntegrator->
    integrateSolutionFunctorOnGeo<ConvectiveVarSet::Flux>
    (ddata.cell, getMethodData().getLinearVar()->getFlux(), m_pdata, m_phi);

  m_flux_space = m_phi;

}

//////////////////////////////////////////////////////////////////////////////

void STM_CRD_SplitStrategy::computeSTFluxIntegral_past(){

  DistributionData& ddata = getMethodData().getDistributionData();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();
  const InwardNormalsData& normalsData = *normals[ddata.cellID];

  // set the face normals
  const CFuint nbFaces = normalsData.nbFaces();

  // Compute the contour integral of the flux
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      m_unitFaceNormals[iFace][iDim] =
        normalsData.getFaceUnitNormComp(iFace,iDim);
    }
  }

  // cell states substituted by the corresponding linearization
  // variables in which the interpolation will be performed
  // before the contour integration
  const CFuint nbCellStates = ddata.states->size();
  for (CFuint iState = 0; iState < nbCellStates; ++iState) {
    ddata.cell->setState(iState, _pastStates[iState]);
  }

  // first do interpolation then compute the physical related data for
  // the interpolated states
  m_contourIntegrator->computeAllAtQuadraturePointsOnGeo(ddata.cell);

  // compute the physical data for each interpolated state
  const CFuint nbQ = m_nbQPointsInCell[getMethodData().getDistributionData().cellID];
  computeStatesData(nbQ, getMethodData().getLinearVar(), m_qdstates,  m_pdata, m_qdExtraVars);
  
  // compute the contour integration of the fluxes
  m_contourIntegrator->
    integrateSolutionFunctorOnGeo<ConvectiveVarSet::Flux>
    (ddata.cell, getMethodData().getLinearVar()->getFlux(), m_pdata, m_phi);

  m_flux_space = m_phi;

  for (CFuint iState = 0; iState < nbCellStates; ++iState) {
    ddata.cell->setState(iState, m_statesBkp[iState]);
  }

}

//////////////////////////////////////////////////////////////////////////////

void STM_CRD_SplitStrategy::unsetup()
{
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

void STM_CRD_SplitStrategy::setup()
{
  CFAUTOTRACE;

  FluctuationSplitStrategy::setup();

  _null.resize(PhysicalModelStack::getActive()->getDim());
  _null = 0.;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // get the splitter
  m_splitter = getMethodData().getSplitter().d_castTo<SpaceTime_Splitter>();
  m_splitter->setUpdateCoeff(updateCoeff);

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  _pastStates.resize(maxNbStatesInCell);
  // Resizing pastStates
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _pastStates[i] = new State();
   }

  // Resizing m_flux_past, m_flux, m_phi
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_flux_past_space.resize(nbEqs);
  m_flux_space.resize(nbEqs);
  m_flux_past_time.resize(nbEqs);
  m_flux_time.resize(nbEqs);
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
   CFuint nbStatesInCell;
  // back up cell state pointers
  m_statesBkp.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());

  // Resizing the storage for the full past Residuals
  _pastResiduals_order1.resize(nbGeoEnts);

  // Resizing the storage for the full past Residuals
  _pastResiduals.resize(nbGeoEnts);
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

    nbStatesInCell = cell.nbStates();
    _pastResiduals_order1[iGeoEnt].resize(nbStatesInCell*nbEqs);
   _pastResiduals[iGeoEnt].resize(nbStatesInCell*nbEqs);
    // release the GeometricEntity
    _stdTrsGeoBuilder.releaseGE();
  }

  temp_residual.resize(maxNbStatesInCell);

  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState)
    temp_residual[iState].resize(nbEqs);

  getMethodData().getDistributionData().isHO = false;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
STM_CRD_SplitStrategy::needsSockets()
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

void STM_CRD_SplitStrategy::setCurrentCell()
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
