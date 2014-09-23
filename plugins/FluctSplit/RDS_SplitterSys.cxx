#include "RDS_SplitterSys.hh"
#include "Common/CFLog.hh"
#include "MathTools/MatrixInverter.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

RDS_SplitterSys::RDS_SplitterSys(const std::string& name) :
  Splitter(name),
  socket_updateCoeff("updateCoeff"),
  m_eValues(0),
  _tempKp(),
  _tempKm(),
  _inverter(CFNULL),
  _kPlus(0),
  _kMin(0),
  m_kCoeff(0.),
  m_nodeArea(0.),
  m_maxNbStatesInCell(0),  
  m_shockFaceNodes(2),
  m_delta(),
  m_cellBuilder(),
  m_mapNode2CellID(),
  m_nxShock(0.),
  m_nyShock(0.),
  m_lambda(),
  m_nShock(),
  m_eValuesP(0),
  m_eValuesM(0),
  m_rightEv(),
  m_leftEv()
{
}

//////////////////////////////////////////////////////////////////////////////

RDS_SplitterSys::~RDS_SplitterSys()
{
  for (CFuint i = 0; i < m_eValues.size(); ++i) {
    deletePtr(m_eValues[i]);
  }

  for (CFuint i = 0; i < _kPlus.size(); ++i) {
    deletePtr(_kPlus[i]);
  }

  for (CFuint i = 0; i < _kMin.size(); ++i) {
    deletePtr(_kMin[i]);
  }

  deletePtr(_inverter);
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterSys::setBlockData()
{
  CFLogDebugMin( "RDS_SplitterSys::setBlockData() => " <<
		 "blockSeparator = " << _blockSeparator << "\n" <<
		 "totalNbEq = " <<  PhysicalModelStack::getActive()->getNbEq() << "\n");

  cf_assert(_blockSeparator <=  PhysicalModelStack::getActive()->getNbEq() );
  
  if (!getMethodData().isScalarFirst()) {
    _nbEquations = _blockSeparator;
    _firstVarID = 0;
    _lastVarID = _blockSeparator;
    
  }
  else {
    _firstVarID = _blockSeparator;
    _lastVarID = PhysicalModelStack::getActive()->getNbEq();
    _nbEquations = _lastVarID - _firstVarID;
  }
  
  cf_assert(_nbEquations == _lastVarID - _firstVarID);
  cf_assert(_lastVarID <= PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterSys::setup()
{
  CFAUTOTRACE;

  Splitter::setup();
  
  m_maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_kCoeff = 1.0 / PhysicalModelStack::getActive()->getDim();
  
  CFLogDebugMin( "RDS_SplitterSys::_nbEquations: "
  << _nbEquations
  << ", _firstVar = " << _firstVarID
  << ", _lastVar = " << _lastVarID
  << "\n");

  _tempKp.resize(_nbEquations,_nbEquations);
  _tempKm.resize(_nbEquations,_nbEquations);
  _inverter = MatrixInverter::create(_nbEquations, false);

  const CFuint maxNbStatesInCell =
    MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  _kPlus.resize(maxNbStatesInCell);
  _kMin.resize(maxNbStatesInCell);
  m_eValues.resize(maxNbStatesInCell);

  m_eValuesP.resize(_nbEquations);
  m_eValuesM.resize(_nbEquations);
  m_rightEv.resize(_nbEquations,_nbEquations);
  m_leftEv.resize(_nbEquations,_nbEquations);

  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    m_eValues[i] = new RealVector(PhysicalModelStack::getActive()->getNbEq());

    _kPlus[i] = new RealMatrix(_nbEquations, _nbEquations);
    _kMin[i]  = new RealMatrix(_nbEquations, _nbEquations);
  }

  m_delta.resize(maxNbStatesInCell);
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

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  // carbuncle fix
  m_lambda.resize(nbEqs);
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    m_lambda[iEq].resize(3,3);
  }

  m_nShock.resize(3);
  for (CFuint i = 0; i < 3; ++i) {
    m_nShock[i].resize(dim);
  }
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterSys::doComputeK(CFuint iState)
{
  Splitter::setAdimensionalNormal(iState);

  getMethodData().getDistribVar()->splitJacobian(*_kPlus[iState],
					      *_kMin[iState],
					      *m_eValues[iState],
					      _adimNormal);
    
  CFuint nbEqs =  (*m_eValues[iState]).size();
//   First we check if we are at a stagnation point which coorespond 
//   to some eigen value that are null
  bool istagnpoint = false;
  const CFreal Eps = 1.0e-13;//1.0e-13; //1.0e-8
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq){
    if ( std::abs( (*m_eValues[iState])[iEq] ) < Eps )
     istagnpoint = true;
  }
  

   m_nodeArea = m_normals->getAreaNode(iState);
// 
  if (!istagnpoint){
    *_kPlus[iState] *= m_kCoeff * m_nodeArea;
    *_kMin[iState]  *= m_kCoeff * m_nodeArea;
  }
  else {
    getMethodData().getDistribVar()->computeEigenValuesVectors(m_rightEv, m_leftEv, *m_eValues[iState],_adimNormal);
    
    // If we are at a stagnation point we add an epsilon to all the eigenvalue to be sure that
    // the K will keep invertible
    for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      if ((*m_eValues[iState])[iEq] >= 0.)
        (*m_eValues[iState])[iEq] += Eps;
      else
        (*m_eValues[iState])[iEq] -= Eps;
        
      m_eValuesP[iEq] = m_kCoeff * m_nodeArea*max(0.,(*m_eValues[iState])[iEq]);
      m_eValuesM[iEq] = m_kCoeff * m_nodeArea*min(0.,(*m_eValues[iState])[iEq]);
    }

    // compute jacobian + and -
    *_kPlus[iState] = m_rightEv*(m_eValuesP*m_leftEv);
    *_kMin[iState]  = m_rightEv*(m_eValuesM*m_leftEv);
  }
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterSys::computeK(const vector<State*>& states,
			       const InwardNormalsData* const normalsData)
{
  m_normals = normalsData;
  _nbStatesInCell = states.size();

  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  // apply the entropy or carbuncle fix
  getMethodData().getJacobianFixComputer()->computeFix(*normalsData, m_delta);

  // applyCarbuncleFix(states, normalsData);
  //  applyCarbuncleFix();

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    getMethodData().getDistribVar()->setDelta(m_delta[iState]);
    doComputeK(iState);

    if (!getMethodData().getDistributionData().isPerturb) {
      //    const CFreal maxEigenValue = std::max(0.0, m_eValues[iState]->max());
      //    CFreal maxEigenValue = 0.0;
      //for (CFuint i = 0; i < _nbEquations; ++i) {
      //maxEigenValue = max(maxEigenValue, 0.5*(*m_eValues[iState])[i] + std::abs((*m_eValues[iState])[i]));
      //}

      //// remove this
      //    const CFreal u = (*states[iState])[1];
      //const CFreal v = (*states[iState])[2];
      // const CFreal V2 = u*u + v*v;
      //const CFreal T = (*states[iState])[3];
      // const CFreal a = sqrt(1.4*287.046*T);
      // const CFreal Vn = u*_adimNormal[XX] + v*_adimNormal[YY];
      //const CFreal maxEigenValue = std::max(1e-8, a + std::abs(Vn));
      // const CFreal maxEigenValue = 0.5*(a + Vn + std::abs(a + Vn));

      //// remove this
      //    const CFreal rho = (*states[iState])[0];
      // const CFreal u = (*states[iState])[1]/rho;
      //const CFreal v = (*states[iState])[2]/rho;
      // const CFreal V2 = u*u + v*v;
      //const CFreal p = 0.4*((*states[iState])[3] - 0.5*rho*V2);
      // const CFreal a = sqrt(1.4*p/rho);
      //const CFreal Vn = u*_adimNormal[XX] + v*_adimNormal[YY];
      //const CFreal maxEigenValue = std::max(0.0, a + std::abs(Vn));
      const CFreal maxEigenValue = std::max(0.0, m_eValues[iState]->max());
      updateCoeff[states[iState]->getLocalID()] += m_kCoeff*m_nodeArea*maxEigenValue;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterSys::applyCarbuncleFix(const vector<State*>& states,
					const InwardNormalsData* const normalsData)
{
  std::cout << "Deprecated method: void RDS_SplitterSys::applyCarbuncleFix()\n"; abort();
  
//   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
// 
//   // carbuncle fix
//   vector<RealVector> lambda(nbEqs);
//   for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
//     lambda[iEq].resize(3);
//   }
// 
//   const CFreal rho0 =(*states[0])[0];
//   const CFreal u0 =(*states[0])[1]/rho0;
//   const CFreal v0 = (*states[0])[2]/rho0;
//   const CFreal V20 = u0*u0 + v0*v0;
//   const CFreal p0 = 0.4*((*states[0])[3]-rho0*0.5*V20);
//   const CFreal a0 = sqrt(1.4*p0/rho0);
//   const CFreal M0 = sqrt(V20)/a0;
// 
//   CFreal maxMach = M0;
//   CFreal minMach = M0;
// 
//   const CFreal rho1 =(*states[1])[0];
//   const CFreal u1 =(*states[1])[1]/rho1;
//   const CFreal v1 = (*states[1])[2]/rho1;
//   const CFreal V21 = u1*u1 + v1*v1;
//   const CFreal p1 = 0.4*((*states[1])[3]-rho1*0.5*V21);
//   const CFreal a1 = sqrt(1.4*p1/rho1);
//   const CFreal M1 = sqrt(V21)/a1;
// 
//   maxMach = max(maxMach, M1);
//   minMach = min(maxMach, M1);
// 
//   const CFreal rho2 =(*states[2])[0];
//   const CFreal u2 =(*states[2])[1]/rho2;
//   const CFreal v2 = (*states[2])[2]/rho2;
//   const CFreal V22 = u2*u2 + v2*v2;
//   const CFreal p2 = 0.4*((*states[2])[3]-rho2*0.5*V22);
//   const CFreal a2 = sqrt(1.4*p2/rho2);
//   const CFreal M2 = sqrt(V22)/a2;
// 
//   maxMach = max(maxMach, M2);
//   minMach = min(maxMach, M2);
// 
//   const CFreal pGradTimesVolX2 =
//     p0*m_normals->getNodalNormComp(0,XX) +
//     p1*m_normals->getNodalNormComp(1,XX) +
//     p2*m_normals->getNodalNormComp(2,XX);
// 
//   const CFreal pGradTimesVolY2 =
//     p0*m_normals->getNodalNormComp(0,YY) +
//     p1*m_normals->getNodalNormComp(1,YY) +
//     p2*m_normals->getNodalNormComp(2,YY);
// 
//   CFreal maxNormalGradP = 0.0;
//   CFuint nMaxGradPID = 0;
//   for (CFuint is = 0; is < states.size(); ++is) {
//     const CFreal gradP =
//       std::abs(pGradTimesVolX2*m_normals->getNodalUnitNormComp(is,XX) +
// 	       pGradTimesVolY2*m_normals->getNodalUnitNormComp(is,YY));
//     if (gradP > maxNormalGradP) {
//       maxNormalGradP = gradP;
//       nMaxGradPID = is;
//     }
//   }
// 
//   // n_shock corresponds to the unit normal along which
//   // the maximum pressure gradient is detected
//   CFreal nx = m_normals->getNodalUnitNormComp(0,XX);
//   CFreal ny = m_normals->getNodalUnitNormComp(0,YY);
//   lambda[0][0] = lambda[1][0] = u0*nx + v0*ny;
//   lambda[2][0] = u0*nx + v0*ny + a0;
//   lambda[3][0] = u0*nx + v0*ny - a0;
// 
//   nx = m_normals->getNodalUnitNormComp(1,XX);
//   ny = m_normals->getNodalUnitNormComp(1,YY);
//   lambda[0][1] = lambda[1][1] = u1*nx + v1*ny;
//   lambda[2][1] = u1*nx + v1*ny + a1;
//   lambda[3][1] = u1*nx + v1*ny - a1;
// 
//   nx = m_normals->getNodalUnitNormComp(2,XX);
//   ny = m_normals->getNodalUnitNormComp(2,YY);
//   lambda[0][2] = lambda[1][2] = u2*nx + v2*ny;
//   lambda[2][2] = u2*nx + v2*ny + a2;
//   lambda[3][2] = u2*nx + v2*ny - a2;
// 
//   //  const CFreal nxParallelShock = -m_normals->getNodalUnitNormComp(nMaxGradPID,YY);
//   //  const CFreal nyParallelShock = m_normals->getNodalUnitNormComp(nMaxGradPID,XX);
// 
//   CFreal delta = 0.0;
//   for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
//     for (CFuint is = 0; is < 3; ++is) {
//       CFreal lambdaUpDen = 0.0;
//       CFreal lambdaUpNum = 0.0;
//       CFreal lambdaDownDen = 0.0;
//       CFreal lambdaDownNum = 0.0;
// 
//       for (CFuint i = 0; i < 3; ++i) {
// 	const CFreal kj =
// 	  m_normals->getNodalUnitNormComp(is,XX)*m_normals->getNodalNormComp(i,XX) +
// 	  m_normals->getNodalUnitNormComp(is,YY)*m_normals->getNodalNormComp(i,YY);
// 	const CFreal kmin = min(0.,kj);
// 	const CFreal kplus = max(0.,kj);
// 	lambdaUpDen += kmin;
// 	lambdaUpNum += kmin*lambda[iEq][i];
// 	lambdaDownDen += kplus;
// 	lambdaDownNum += kplus*lambda[iEq][i];
//       }
//       const CFreal lambdaUp = lambdaUpNum/lambdaUpDen;
//       const CFreal lambdaDown = lambdaDownNum/lambdaDownDen;
//       delta = max(delta, std::abs(lambdaDown - lambdaUp));
//     }
//   }
//   delta *= 0.5;
// 
//   const CFreal nxShock = m_normals->getNodalUnitNormComp(nMaxGradPID,XX);
//   const CFreal nyShock = m_normals->getNodalUnitNormComp(nMaxGradPID,YY);
//   for (CFuint is = 0; is < 3; ++is) {
//     const CFreal etaCorr = std::abs(-m_normals->getNodalUnitNormComp(is,XX)*nyShock +
// 				    m_normals->getNodalUnitNormComp(is,YY)*nxShock);
// 
//     m_delta[is] = delta*etaCorr;
//   }
}

//////////////////////////////////////////////////////////////////////////////

// void RDS_SplitterSys::updateDelta(GeometricEntity *const cell, bool reuseNShock)
// {
//   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
//   // carbuncle fix
//   vector<RealMatrix> lambda(nbEqs);
//   for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
//     lambda[iEq].resize(3,3);
//   }

//   const vector<State*>& states = *cell->getStates();
//   const CFreal rho0 =(*states[0])[0];
//   const CFreal u0 =(*states[0])[1]/rho0;
//   const CFreal v0 = (*states[0])[2]/rho0;
//   const CFreal V20 = u0*u0 + v0*v0;
//   const CFreal p0 = 0.4*((*states[0])[3]-rho0*0.5*V20);
//   const CFreal a0 = sqrt(1.4*p0/rho0);
//   const CFreal M0 = sqrt(V20)/a0;

//   CFreal maxMach = M0;
//   CFreal minMach = M0;

//   const CFreal rho1 =(*states[1])[0];
//   const CFreal u1 =(*states[1])[1]/rho1;
//   const CFreal v1 = (*states[1])[2]/rho1;
//   const CFreal V21 = u1*u1 + v1*v1;
//   const CFreal p1 = 0.4*((*states[1])[3]-rho1*0.5*V21);
//   const CFreal a1 = sqrt(1.4*p1/rho1);
//   const CFreal M1 = sqrt(V21)/a1;

//   maxMach = max(maxMach, M1);
//   minMach = min(maxMach, M1);

//   const CFreal rho2 =(*states[2])[0];
//   const CFreal u2 =(*states[2])[1]/rho2;
//   const CFreal v2 = (*states[2])[2]/rho2;
//   const CFreal V22 = u2*u2 + v2*v2;
//   const CFreal p2 = 0.4*((*states[2])[3]-rho2*0.5*V22);
//   const CFreal a2 = sqrt(1.4*p2/rho2);
//   const CFreal M2 = sqrt(V22)/a2;

//   maxMach = max(maxMach, M2);
//   minMach = min(maxMach, M2);

//   CFuint nMaxGradPID = 0;
//   if (!reuseNShock) {
//     const CFreal pGradTimesVolX2 =
//       p0*m_normals->getNodalUnitNormComp(0,XX) +
//       p1*m_normals->getNodalUnitNormComp(1,XX) +
//       p2*m_normals->getNodalUnitNormComp(2,XX);

//     const CFreal pGradTimesVolY2 =
//       p0*m_normals->getNodalUnitNormComp(0,YY) +
//       p1*m_normals->getNodalUnitNormComp(1,YY) +
//       p2*m_normals->getNodalUnitNormComp(2,YY);

//     CFreal maxNormalGradP = 0.0;
//     CFuint nMaxGradPID = 0;
//     for (CFuint is = 0; is < states.size(); ++is) {
//       const CFreal gradP = pGradTimesVolX2*m_normals->getNodalUnitNormComp(is,XX) +
// 	pGradTimesVolY2*m_normals->getNodalUnitNormComp(is,YY);
//       if (gradP > maxNormalGradP) {
// 	maxNormalGradP = gradP;
// 	nMaxGradPID = is;
//       }
//     }

//     if (nMaxGradPID == 0) {
//       m_shockFaceNodes[0] = 1;
//       m_shockFaceNodes[1] = 2;
//     }
//     if (nMaxGradPID == 1) {
//       m_shockFaceNodes[0] = 2;
//       m_shockFaceNodes[1] = 0;
//     }
//     if (nMaxGradPID == 2) {
//       m_shockFaceNodes[0] = 0;
//       m_shockFaceNodes[1] = 1;
//     }
//   }

//   // n_shock corresponds to the unit normal along which
//   // the maximum pressure gradient is detected
//   for (CFuint is = 0; is < 3; ++is) {
//     const CFreal nx = m_normals->getNodalUnitNormComp(is,XX);
//     const CFreal ny = m_normals->getNodalUnitNormComp(is,YY);

//     lambda[0](is,0) = lambda[1](is,0) = u0*nx + v0*ny;
//     lambda[0](is,1) = lambda[1](is,1) = u1*nx + v1*ny;
//     lambda[0](is,2) = lambda[1](is,2) = u2*nx + v2*ny;

//     lambda[2](is,0) = u0*nx + v0*ny + a0;
//     lambda[2](is,1) = u1*nx + v1*ny + a1;
//     lambda[2](is,2) = u2*nx + v2*ny + a2;

//     lambda[3](is,0) = u0*nx + v0*ny - a0;
//     lambda[3](is,1) = u1*nx + v1*ny - a1;
//     lambda[3](is,2) = u2*nx + v2*ny - a2;
//   }

//   CFreal delta = 0.0;
//   for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
//     for (CFuint is = 0; is < 3; ++is) {
//       CFreal lambdaUpDen = 0.0;
//       CFreal lambdaUpNum = 0.0;
//       CFreal lambdaDownDen = 0.0;
//       CFreal lambdaDownNum = 0.0;

//       for (CFuint i = 0; i < 3; ++i) {
// 	const CFreal kj =
// 	  m_normals->getNodalNormComp(is,XX)*
// 	  m_normals->getNodalNormComp(i,XX) +
// 	  m_normals->getNodalNormComp(is,YY)*
// 	  m_normals->getNodalNormComp(i,YY);
// 	const CFreal kmin = min(0.,kj);
// 	const CFreal kplus = max(0.,kj);
// 	lambdaUpDen += kmin;
// 	lambdaUpNum += kmin*lambda[iEq](is,i);
// 	lambdaDownDen += kplus;
// 	lambdaDownNum += kplus*lambda[iEq](is,i);
//       }
//       const CFreal lambdaUp = lambdaUpNum/lambdaUpDen;
//       const CFreal lambdaDown = lambdaDownNum/lambdaDownDen;
//       delta = max(delta, std::abs(lambdaDown - lambdaUp));
//     }
//   }
//   delta *= 0.5;

//   // (inward) normal to shock in the other cell has opposite sign
//   m_nxShock = (!reuseNShock) ? m_normals->getNodalUnitNormComp(nMaxGradPID,XX) : -m_nxShock;
//   m_nyShock = (!reuseNShock) ? m_normals->getNodalUnitNormComp(nMaxGradPID,YY) : -m_nyShock;
//   for (CFuint is = 0; is < 3; ++is) {
//     // const CFreal etaCorr = std::abs(-m_normals->getNodalNormComp(is,XX)*m_nxShock +
//     //   m_normals->getNodalNormComp(is,YY)*m_nyShock);

//     // m_delta[is] = max(m_delta[is], delta*etaCorr);
//     m_delta[is] = max(m_delta[is], delta);
//   }
// }

//////////////////////////////////////////////////////////////////////////////

// void RDS_SplitterSys::applyCarbuncleFix()
// {
//   DistributionData& ddata = getMethodData().getDistributionData();
//   GeometricEntity* const currCell = ddata.cell;
//   updateDelta(currCell, false);

//   vector<CFuint> c1Vec;
//   vector<CFuint> c2Vec;
//   typedef CFMultiMap<CFuint,CFuint> MapNode2Cell;
//   typedef  MapNode2Cell::MapIterator MapIt;

//   pair<MapIt, MapIt> cells1 = m_mapNode2CellID.
//     find(currCell->getState(m_shockFaceNodes[0])->getLocalID());
//   for (MapIt c1Ptr = cells1.first; c1Ptr != cells1.second; ++c1Ptr) {
//     c1Vec.push_back(c1Ptr->second);
//   }

//   pair<MapIt, MapIt> cells2 = m_mapNode2CellID.
//     find(currCell->getState(m_shockFaceNodes[1])->getLocalID());
//   for (MapIt c2Ptr = cells2.first; c2Ptr != cells2.second; ++c2Ptr) {
//     c2Vec.push_back(c2Ptr->second);
//   }

//   CFint neighborCellID = -1;
//   for (CFuint i = 0; i < c1Vec.size(); ++i) {
//     for (CFuint j = 0; j < c2Vec.size(); ++j) {
//       if (c1Vec[i] == c2Vec[j] && c1Vec[i] != currCell->getID()) {
// 	neighborCellID = c1Vec[i];
//        	break;
//       }
//     }
//     if (neighborCellID > -1) break;
//   }

//   cf_assert(neighborCellID != currCell->getID());

//   if (neighborCellID > -1) {
//     StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder.getDataGE();
//     geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");
//     geoData.idx = static_cast<CFuint>(neighborCellID);
//     GeometricEntity* const neighborCell = m_cellBuilder.buildGE();
//     updateDelta(neighborCell, true);
//     //release the GeometricEntity
//     m_cellBuilder.releaseGE();
//   }
// }

// //////////////////////////////////////////////////////////////////////////////

void RDS_SplitterSys::applyCarbuncleFix()
{
  std::cout << "Deprecated method: void RDS_SplitterSys::applyCarbuncleFix()\n"; abort();
//   DistributionData& ddata = getMethodData().getDistributionData();
//   GeometricEntity* const currCell = ddata.cell;
//   computeCellDelta(currCell,-1);
//   m_shockFaceNodes[0] = 0; m_shockFaceNodes[1] = 1;
//   updateDelta(2);
// 
//   m_shockFaceNodes[0] = 1; m_shockFaceNodes[1] = 2;
//   updateDelta(0);
// 
//   m_shockFaceNodes[0] = 2; m_shockFaceNodes[1] = 0;
//   updateDelta(1);
// 
//   // take the half of the computed delta
//   m_delta *= 0.5;
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterSys::computeCellDelta(GeometricEntity *const cell, CFint iState)
{
  std::cout << "Deprecated method: void RDS_SplitterSys::computeCellDelta(GeometricEntity *const cell, CFint iState)\n"; abort();
//   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
// 
//   const vector<State*>& states = *cell->getStates();
//   const CFreal rho0 =(*states[0])[0];
//   const CFreal u0 =(*states[0])[1]/rho0;
//   const CFreal v0 = (*states[0])[2]/rho0;
//   const CFreal V20 = u0*u0 + v0*v0;
//   const CFreal p0 = 0.4*((*states[0])[3]-rho0*0.5*V20);
//   const CFreal a0 = sqrt(1.4*p0/rho0);
//   // const CFreal M0 = sqrt(V20)/a0;
// 
//   // CFreal maxMach = M0;
//   // CFreal minMach = M0;
// 
//   const CFreal rho1 =(*states[1])[0];
//   const CFreal u1 =(*states[1])[1]/rho1;
//   const CFreal v1 = (*states[1])[2]/rho1;
//   const CFreal V21 = u1*u1 + v1*v1;
//   const CFreal p1 = 0.4*((*states[1])[3]-rho1*0.5*V21);
//   const CFreal a1 = sqrt(1.4*p1/rho1);
//   // const CFreal M1 = sqrt(V21)/a1;
// 
//   // maxMach = max(maxMach, M1);
//   // minMach = min(maxMach, M1);
// 
//   const CFreal rho2 =(*states[2])[0];
//   const CFreal u2 =(*states[2])[1]/rho2;
//   const CFreal v2 = (*states[2])[2]/rho2;
//   const CFreal V22 = u2*u2 + v2*v2;
//   const CFreal p2 = 0.4*((*states[2])[3]-rho2*0.5*V22);
//   const CFreal a2 = sqrt(1.4*p2/rho2);
//   // const CFreal M2 = sqrt(V22)/a2;
// 
//   // maxMach = max(maxMach, M2);
//   // minMach = min(maxMach, M2);
// 
//   // compute the lambda_ik (n_j) corresponding to all
//   // each eigen value component in each state projected in the
//   // direction of each inward normal
//   for (CFuint is = 0; is < 3; ++is) {
//     const CFreal nx = m_normals->getNodalUnitNormComp(is,XX);
//     const CFreal ny = m_normals->getNodalUnitNormComp(is,YY);
// 
//     m_lambda[0](is,0) = m_lambda[1](is,0) = u0*nx + v0*ny;
//     m_lambda[0](is,1) = m_lambda[1](is,1) = u1*nx + v1*ny;
//     m_lambda[0](is,2) = m_lambda[1](is,2) = u2*nx + v2*ny;
// 
//     m_lambda[2](is,0) = u0*nx + v0*ny + a0;
//     m_lambda[2](is,1) = u1*nx + v1*ny + a1;
//     m_lambda[2](is,2) = u2*nx + v2*ny + a2;
// 
//     m_lambda[3](is,0) = u0*nx + v0*ny - a0;
//     m_lambda[3](is,1) = u1*nx + v1*ny - a1;
//     m_lambda[3](is,2) = u2*nx + v2*ny - a2;
//   }
// 
//   CFreal delta = 0.0;
//   for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
//     for (CFuint is = 0; is < 3; ++is) {
//       CFreal lambdaUpDen = 0.0;
//       CFreal lambdaUpNum = 0.0;
//       CFreal lambdaDownDen = 0.0;
//       CFreal lambdaDownNum = 0.0;
// 
//       for (CFuint i = 0; i < 3; ++i) {
// 	const CFreal kj =
// 	  m_normals->getNodalNormComp(is,XX)*
// 	  m_normals->getNodalNormComp(i,XX) +
// 	  m_normals->getNodalNormComp(is,YY)*
// 	  m_normals->getNodalNormComp(i,YY);
// 	const CFreal kmin = min(0.,kj);
// 	const CFreal kplus = max(0.,kj);
// 	lambdaUpDen += kmin;
// 	lambdaUpNum += kmin*m_lambda[iEq](is,i);
// 	lambdaDownDen += kplus;
// 	lambdaDownNum += kplus*m_lambda[iEq](is,i);
//       }
//       const CFreal lambdaUp = lambdaUpNum/lambdaUpDen;
//       const CFreal lambdaDown = lambdaDownNum/lambdaDownDen;
//       delta = max(delta, std::abs(lambdaDown - lambdaUp));
//     }
//   }
// 
//   if (iState > -1) {
//     // just update the iState component of the delta
//     m_delta[iState] = max(m_delta[iState], delta);
//   }
//   else {
//     // compute the reference delta for each state
//     m_delta = delta;
//   }
}

//////////////////////////////////////////////////////////////////////////////

void RDS_SplitterSys::updateDelta(CFuint iState)
{
  std::cout << "Deprecated method: void RDS_SplitterSys::updateDelta(CFuint iState)\n"; abort();
//   CFint neighborCellID = -1;
// 
//   typedef CFMultiMap<CFuint,CFuint> MapNode2Cell;
//   typedef  MapNode2Cell::MapIterator MapIt;
// 
//   GeometricEntity* const currCell =
//     getMethodData().getDistributionData().cell;
// 
//   vector<CFuint> c1Vec;
//   vector<CFuint> c2Vec;
//   typedef CFMultiMap<CFuint,CFuint> MapNode2Cell;
//   typedef  MapNode2Cell::MapIterator MapIt;
//   
//   bool fo = false;
//   pair<MapIt, MapIt> cells1 = m_mapNode2CellID.
//     find(currCell->getState(m_shockFaceNodes[0])->getLocalID(), fo);
//   cf_assert(fo);
//   
//   for (MapIt c1Ptr = cells1.first; c1Ptr != cells1.second; ++c1Ptr) {
//     c1Vec.push_back(c1Ptr->second);
//   }
//   
//   fo = false;
//   pair<MapIt, MapIt> cells2 = m_mapNode2CellID.
//     find(currCell->getState(m_shockFaceNodes[1])->getLocalID(), fo);
//   cf_assert(fo);
//   
//   for (MapIt c2Ptr = cells2.first; c2Ptr != cells2.second; ++c2Ptr) {
//     c2Vec.push_back(c2Ptr->second);
//   }
//   
//   for (CFuint i = 0; i < c1Vec.size(); ++i) {
//     for (CFuint j = 0; j < c2Vec.size(); ++j) {
//       if (c1Vec[i] == c2Vec[j] && c1Vec[i] != currCell->getID()) {
// 	neighborCellID = c1Vec[i];
//        	break;
//       }
//     }
//     if (neighborCellID > -1) break;
//   }
// 
//   cf_assert(neighborCellID != static_cast<CFint>(currCell->getID()));
// 
//   if (neighborCellID > -1) {
//     StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder.getDataGE();
//     geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");
//     geoData.idx = static_cast<CFuint>(neighborCellID);
//     GeometricEntity* const neighborCell = m_cellBuilder.buildGE();
//     computeCellDelta(neighborCell,iState);
//     //release the GeometricEntity
//     m_cellBuilder.releaseGE();
//   }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
RDS_SplitterSys::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result =
    Splitter::needsSockets();

  result.push_back(&socket_updateCoeff);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

