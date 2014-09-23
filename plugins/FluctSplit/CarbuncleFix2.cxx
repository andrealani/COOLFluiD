#include "Framework/GeometricEntity.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"

#include "FluctSplit/CarbuncleFix2.hh"
#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<CarbuncleFix2,
                       FluctuationSplitData,
                       ComputeJacobianFix,
                       FluctSplitNavierStokesModule>
carbuncleFix2Provider("Carbuncle2");

//////////////////////////////////////////////////////////////////////////////

CarbuncleFix2::CarbuncleFix2(const std::string& name) :
  ComputeJacobianFix(name),
  m_normals(CFNULL),
  m_delta(CFNULL),
  m_pGradTimesVolDim(),
  m_nShock(),
  m_pData(),
  m_lambda(),
  m_avSpeed(),
  m_shockFaceNodes(2),
  m_cellBuilder(),
  m_mapNode2CellID()
{
}

//////////////////////////////////////////////////////////////////////////////

CarbuncleFix2::~CarbuncleFix2()
{
}

//////////////////////////////////////////////////////////////////////////////

void CarbuncleFix2::setup()
{
  ComputeJacobianFix::setup();
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_pGradTimesVolDim.resize(dim);
  m_nShock.resize(dim);
  
  const CFuint nbCellStates = MeshDataStack::getActive()->
    Statistics().getMaxNbStatesInCell();
  m_pData.resize(nbCellStates);
  for (CFuint i = 0; i < nbCellStates; ++i) {
    PhysicalModelStack::getActive()->getImplementor()->
      getConvectiveTerm()->resizePhysicalData(m_pData[i]);
  }  
  
  m_lambda.resize(nbCellStates);
  for (CFuint i = 0; i < nbCellStates; ++i) {
    m_lambda[i].resize(nbEqs);
  }
  m_avSpeed.resize(dim);

  // set up of the local cell builder
  m_cellBuilder.setup();

  // compute all the node-to-cell connectivity
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder.getDataGE();
  geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbGeos = geoData.trs->getLocalNbGeoEnts();
  
  m_mapNode2CellID.reserve
    (nbGeos*MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());

  for (CFuint cellID = 0; cellID < nbGeos; ++cellID) {
    // build the GeometricEntity
    geoData.idx = cellID;
    GeometricEntity *const cell = m_cellBuilder.buildGE();
    vector<State*> *const states = cell->getStates();
    for (CFuint is = 0; is < states->size(); ++is) {
      m_mapNode2CellID.insert((*states)[is]->getLocalID(), cellID);
    }
    m_cellBuilder.releaseGE();
  }
  m_mapNode2CellID.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void CarbuncleFix2::computeFix(const InwardNormalsData& normalsData,
			       RealVector& delta)
{
  DistributionData& data = getMethodData().getDistributionData();
  GeometricEntity* const currCell = data.cell;
  m_delta = &delta;
  m_normals = const_cast<InwardNormalsData*>(&normalsData);
  
  computeCellDelta(currCell,-1);
  m_shockFaceNodes[0] = 0; m_shockFaceNodes[1] = 1;
  updateDelta(2);
  
  m_shockFaceNodes[0] = 1; m_shockFaceNodes[1] = 2;
  updateDelta(0);
  
  m_shockFaceNodes[0] = 2; m_shockFaceNodes[1] = 0;
  updateDelta(1);
  
  // take the half of the computed delta
  *m_delta *= 0.5;
  
  // try this
  const CFreal avV = sqrt(m_avSpeed[XX]*m_avSpeed[XX] + 
			    m_avSpeed[YY]*m_avSpeed[YY]);
  const CFreal cosD = m_avSpeed[XX]/avV;
  const CFreal sinD = m_avSpeed[YY]/avV;
  const CFuint nbStates = data.cell->nbStates();
  for (CFuint is = 0; is < nbStates; ++is) {
    // only 2D case for now
    const CFreal nEtaScalarNormal = -normalsData.getNodalUnitNormComp(is,XX)*sinD + 
      normalsData.getNodalUnitNormComp(is,YY)*cosD;
    (*m_delta)[is] = std::abs(nEtaScalarNormal)*(*m_delta)[is];
  }  
}

//////////////////////////////////////////////////////////////////////////////

void CarbuncleFix2::computeCellDelta(GeometricEntity *const cell, 
				     CFint iStateCurr)
{
  // here we are assuming that the eigenvalues of the K matrices 
  // related to update variables are the same as for distribution
  // it is not true in case characteristic variables for HE splitting
  // are used to distribute the residual

// unused //  DistributionData& data = getMethodData().getDistributionData();
  //  if (!data.isPerturb) {
  const CFuint nbStates = cell->nbStates();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  if (iStateCurr == -1) {
    m_pGradTimesVolDim = 0.0;
    m_avSpeed = 0.0;
    for (CFuint is = 0; is < nbStates; ++is) {
      getMethodData().getUpdateVar()->computePhysicalData(*cell->getState(is), m_pData[is]);
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
	m_pGradTimesVolDim[iDim] += m_pData[is][EulerTerm::P]*m_normals->
	  getNodalNormComp(is,iDim); 
      }
      
      m_avSpeed[XX] += m_pData[is][EulerTerm::VX];
      m_avSpeed[YY] += m_pData[is][EulerTerm::VY];
      if (dim == DIM_3D) {
	m_avSpeed[ZZ] += m_pData[is][EulerTerm::VZ];
      }
    }
    m_avSpeed /= nbStates;
  }
  
  for (CFuint is = 0; is < nbStates; ++is) {
    const RealVector& pv = m_pData[is];
    CFreal un =  pv[EulerTerm::VX]*m_normals->getNodalUnitNormComp(is,XX) +
      pv[EulerTerm::VY]*m_normals->getNodalUnitNormComp(is,YY);
    if (dim == DIM_3D) {
      un += pv[EulerTerm::VZ]*m_normals->getNodalUnitNormComp(is,ZZ);
    }
    
    m_lambda[is] = un;
    m_lambda[is][dim]   += pv[EulerTerm::A];
    m_lambda[is][dim+1] -= pv[EulerTerm::A];
  }
  
  CFreal tmpDelta = 0.0;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      CFreal lambdaUpDen = 0.0;
      CFreal lambdaUpNum = 0.0;
      CFreal lambdaDownDen = 0.0;
      CFreal lambdaDownNum = 0.0;
      for (CFuint is = 0; is < nbStates; ++is) {
	CFreal kj = 0.0;
	for (CFuint iDim = 0; iDim < dim; ++iDim) {
	  // should we use here the unit inward normal ?? 
	  kj += m_normals->getNodalUnitNormComp(iState,iDim)*
	    m_normals->getNodalNormComp(is,iDim);
	}
	const CFreal kmin = min(0.,kj);
	const CFreal kplus = max(0.,kj);
	lambdaUpDen += kmin;
	lambdaUpNum += kmin*m_lambda[is][iEq];
	lambdaDownDen += kplus;
	lambdaDownNum += kplus*m_lambda[is][iEq];
      }
      const CFreal lambdaUp = lambdaUpNum/lambdaUpDen;
      const CFreal lambdaDown = lambdaDownNum/lambdaDownDen;
      tmpDelta = max(tmpDelta, std::abs(lambdaDown - lambdaUp));
    }
  }
  
  /*  const CFreal avV = sqrt(m_avSpeed[XX]*m_avSpeed[XX] + 
      m_avSpeed[YY]*m_avSpeed[YY]);
      const CFreal cosD = m_avSpeed[XX]/avV;
      const CFreal sinD = m_avSpeed[YY]/avV;
      cf_assert(dim == DIM_2D);
      
      for (CFuint is = 0; is < nbStates; ++is) {
      // only 2D case for now
      const CFreal nEtaScalarNormal = -m_normals->getNodalUnitNormComp(is,XX)*sinD + 
      m_normals->getNodalUnitNormComp(is,YY)*cosD;
      delta[is] = std::abs(nEtaScalarNormal)*tmpDelta;
      }*/
  
  if (iStateCurr > -1) {
    // just update the iState component of the delta
    (*m_delta)[iStateCurr] = max((*m_delta)[iStateCurr], tmpDelta);
  }
  else {
    // compute the reference delta for each state 
    *m_delta = tmpDelta;
  }
  
  //  }
}

//////////////////////////////////////////////////////////////////////////////

void CarbuncleFix2::updateDelta(CFuint iState)
{
  CFint neighborCellID = -1;
  
  typedef CFMultiMap<CFuint,CFuint> MapNode2Cell;
  typedef MapNode2Cell::MapIterator MapIt;
  
  GeometricEntity* const currCell = 
    getMethodData().getDistributionData().cell;
  
  vector<CFuint> c1Vec;
  vector<CFuint> c2Vec;
  typedef CFMultiMap<CFuint,CFuint> MapNode2Cell;
  typedef MapNode2Cell::MapIterator MapIt;
  
  bool fo = false;
  pair<MapIt, MapIt> cells1 = m_mapNode2CellID.
    find(currCell->getState(m_shockFaceNodes[0])->getLocalID(), fo);
  cf_assert(fo);
  for (MapIt c1Ptr = cells1.first; c1Ptr != cells1.second; ++c1Ptr) {
    c1Vec.push_back(c1Ptr->second);
  }
  
  fo = false;
  pair<MapIt, MapIt> cells2 = m_mapNode2CellID.
    find(currCell->getState(m_shockFaceNodes[1])->getLocalID(), fo);
  cf_assert(fo);
  for (MapIt c2Ptr = cells2.first; c2Ptr != cells2.second; ++c2Ptr) {
    c2Vec.push_back(c2Ptr->second);
  }
  
  for (CFuint i = 0; i < c1Vec.size(); ++i) {
    for (CFuint j = 0; j < c2Vec.size(); ++j) {
      if (c1Vec[i] == c2Vec[j] && c1Vec[i] != currCell->getID()) {
	neighborCellID = c1Vec[i];
       	break;
      }
    }
    if (neighborCellID > -1) break;
  }
  
  cf_assert(neighborCellID != static_cast<CFint>(currCell->getID()));
  
  if (neighborCellID > -1) {
    StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder.getDataGE();
    geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");
    geoData.idx = static_cast<CFuint>(neighborCellID);
    GeometricEntity* const neighborCell = m_cellBuilder.buildGE();
    computeCellDelta(neighborCell,iState);
    //release the GeometricEntity
    m_cellBuilder.releaseGE();
  }
}
//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit



} // namespace COOLFluiD
