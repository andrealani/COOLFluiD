#include "FluctSplit/CRD_SplitStrategy.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<CRD_SplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitModule>
                       crdFluctSplitStrategyProvider("CRD");

//////////////////////////////////////////////////////////////////////////////

CRD_SplitStrategy::CRD_SplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  m_phiT(),
  m_splitter(CFNULL),
  m_unitFaceNormals(),
  m_contourIntegrator(CFNULL),
  m_qdstates(),
  m_qExtraVars(),
  m_nbQPointsInCell(),
  m_statesBkp()
{
}

//////////////////////////////////////////////////////////////////////////////

CRD_SplitStrategy::~CRD_SplitStrategy()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void CRD_SplitStrategy::setup()
{
  // first call parent method
  FluctuationSplitStrategy::setup();

  m_phiT.resize(PhysicalModelStack::getActive()->getNbEq());
  
//   cout << "CRD " << PhysicalModelStack::getActive()->getNbEq() << "\n";

  if (!getMethodData().isMultipleSplitter()) {
    m_splitter = getMethodData().getSplitter();
  }

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
  
  m_qExtraVars.resize(m_qdstates.size());
  
  // this could just be avoided: supposing that the RealVetor* is set from an existing data
  for (CFuint i = 0; i < m_qExtraVars.size(); ++i) {
    m_qExtraVars[i] = new RealVector(getMethodData().getLinearVar()->getExtraPhysicalVarsSize());
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
}

//////////////////////////////////////////////////////////////////////////////

void CRD_SplitStrategy::computeFluctuation(vector<RealVector>& residual)
{
  cf_assert(!getMethodData().isMultipleSplitter());

  setCurrentCell();

  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DistributionData& ddata = getMethodData().getDistributionData();
  
  // compute the residual and the upwind parameters k in this cell
  m_splitter->computeK(*ddata.states, normals[ddata.cellID]);

  computeFluxIntegral();


  const CFuint nbST = getMethodData().getSourceTermSplitter()->size();
  if (nbST == 1) {
    ddata.sourceTermID = 0;
    getMethodData().getSourceTermSplitter(0)->computeSourceTerm(*normals[ddata.cellID]);

    if (getMethodData().includeSourceInFlux()) {
      // in this case converctive and source fluctuations will be distributed together
      m_phiT -= ddata.phiS;
    }

    // transform fluxes + source term to distribute variables
    ddata.phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_phiT);
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

void CRD_SplitStrategy::computeFluxIntegral()
{
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DistributionData& ddata = getMethodData().getDistributionData();
  const InwardNormalsData& normalsData = *normals[ddata.cellID];

  // set the face normals
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbFaces = normalsData.nbFaces();
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
    ddata.cell->setState(iState, (*m_linearStates)[iState]);
  }

  // first do interpolation then compute the physical related data for
  // the interpolated states
  m_contourIntegrator->computeAllAtQuadraturePointsOnGeo(ddata.cell);

  // compute the physical data for each interpolated state
  const CFuint nbQ = m_nbQPointsInCell[getMethodData().getDistributionData().cellID];
  computeStatesData(nbQ, getMethodData().getLinearVar(), m_qdstates,  m_pdata, m_qExtraVars);
 
  // compute the contour integration of the fluxes
  m_contourIntegrator->integrateSolutionFunctorOnGeo<ConvectiveVarSet::Flux>
    (ddata.cell, getMethodData().getLinearVar()->getFlux(), m_pdata, m_phiT);

 
  for (CFuint iState = 0; iState < nbCellStates; ++iState) {
    ddata.cell->setState(iState, m_statesBkp[iState]);
  }

}

//////////////////////////////////////////////////////////////////////////////

void CRD_SplitStrategy::setCurrentCell()
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

void CRD_SplitStrategy::unsetup()
{
  // this could just be avoided: supposing that the RelaVetor* is set from an existing data
  for (CFuint i = 0; i < m_qExtraVars.size(); ++i) {
    deletePtr(m_qExtraVars[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
