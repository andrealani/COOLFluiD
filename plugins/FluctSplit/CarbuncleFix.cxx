#include "Framework/GeometricEntity.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"

#include "FluctSplit/CarbuncleFix.hh"
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

MethodStrategyProvider<CarbuncleFix,
                       FluctuationSplitData,
                       ComputeJacobianFix,
                       FluctSplitNavierStokesModule>
carbuncleFixProvider("Carbuncle");

//////////////////////////////////////////////////////////////////////////////

CarbuncleFix::CarbuncleFix(const std::string& name) :
  ComputeJacobianFix(name),
  m_pGradTimesVolDim(),
  m_nShock(),
  m_pData(),
  m_lambda(),
  m_avSpeed()
{
}

//////////////////////////////////////////////////////////////////////////////

CarbuncleFix::~CarbuncleFix()
{
}

//////////////////////////////////////////////////////////////////////////////

void CarbuncleFix::setup()
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
}

//////////////////////////////////////////////////////////////////////////////

void CarbuncleFix::computeFix(const InwardNormalsData& normalsData,
			      RealVector& delta)
{
  // here we are assuming that the eigenvalues of the K matrices 
  // related to update variables are the same as for distribution
  // it is not true in case characteristic variables for HE splitting
  // are used to distribute the residual
  DistributionData& data = getMethodData().getDistributionData();
  //  if (!data.isPerturb) {
  const CFuint nbStates = data.states->size();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  m_pGradTimesVolDim = 0.0;
  m_avSpeed = 0.0;
  for (CFuint is = 0; is < nbStates; ++is) {
    getMethodData().getUpdateVar()->computePhysicalData(*(*data.states)[is], m_pData[is]);
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      m_pGradTimesVolDim[iDim] += m_pData[is][EulerTerm::P]*normalsData.
	getNodalNormComp(is,iDim); 
    }
    
    m_avSpeed[XX] += m_pData[is][EulerTerm::VX];
    m_avSpeed[YY] += m_pData[is][EulerTerm::VY];
    if (dim == DIM_3D) {
      m_avSpeed[ZZ] += m_pData[is][EulerTerm::VZ];
    }
  }
  m_avSpeed /= nbStates;
    
  /*    CFreal maxNormalGradP = 0.0;
	CFuint nMaxGradPID = 0;
	for (CFuint is = 0; is < nbStates; ++is) {
	CFreal gradP = 0.0;
	for (CFuint iDim = 0; iDim < dim; ++iDim) {
	gradP += m_pGradTimesVolDim[iDim]*normalsData.getNodalUnitNormComp(is,iDim);
	}
	// take the absolute value of the pressure gradient
	gradP = std::abs(gradP);
	if (gradP > maxNormalGradP) {
	maxNormalGradP = gradP;
	nMaxGradPID = is;
	}
	}
	
	// n_shock corresponds to the unit normal along which
	// the maximum pressure gradient is detected
	for (CFuint iDim = 0; iDim < dim; ++iDim) {
	m_nShock[iDim] = normalsData.getNodalUnitNormComp(nMaxGradPID,iDim);
	}
	
	// the normal shock direction has to be downwind with respect 
	// to the cell average velocity
	if (MathFunctions::innerProd(m_avSpeed, m_nShock) < 0.0) {
	m_nShock *= -1.0;
	}
  */
  
  for (CFuint is = 0; is < nbStates; ++is) {
    const RealVector& pv = m_pData[is];
    CFreal un =  pv[EulerTerm::VX]*normalsData.getNodalUnitNormComp(is,XX) +
      pv[EulerTerm::VY]*normalsData.getNodalUnitNormComp(is,YY);
    if (dim == DIM_3D) {
      un += pv[EulerTerm::VZ]*normalsData.getNodalUnitNormComp(is,ZZ);
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
	  kj += normalsData.getNodalUnitNormComp(iState,iDim)*
	    normalsData.getNodalNormComp(is,iDim);
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
  tmpDelta *= 0.5;
  
  const CFreal avV = sqrt(m_avSpeed[XX]*m_avSpeed[XX] + 
			  m_avSpeed[YY]*m_avSpeed[YY]);
  const CFreal cosD = m_avSpeed[XX]/avV;
  const CFreal sinD = m_avSpeed[YY]/avV;
  cf_assert(dim == DIM_2D);
  
  for (CFuint is = 0; is < nbStates; ++is) {
    // only 2D case for now
    const CFreal nEtaScalarNormal = -normalsData.getNodalUnitNormComp(is,XX)*sinD + 
      normalsData.getNodalUnitNormComp(is,YY)*cosD;
    delta[is] = std::abs(nEtaScalarNormal)*tmpDelta;
  }
  //  }
}
      
// //////////////////////////////////////////////////////////////////////////////

// void CarbuncleFix::computeFix(const InwardNormalsData& normalsData,
// 			      RealVector& delta)
// {
//   // here we are assuming that the eigenvalues of the K matrices 
//   // related to update variables are the same as for distribution
//   // it is not true in case characteristic variables for HE splitting
//   // are used to distribute the residual
//   SafePtr<DistributionData> data = getMethodData().getDistributionData();
//   //  if (!data.isPerturb) {
//   const CFuint nbStates = data.states->size();
//   const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
//   vector<CFuint> bigger1;
//   vector<CFuint> smaller1;  
  
//   m_pGradTimesVolDim = 0.0;
//   m_avSpeed = 0.0;
//   for (CFuint is = 0; is < nbStates; ++is) {
//     getMethodData().getUpdateVar()->computePhysicalData(*(*data.states)[is], m_pData[is]);
//     for (CFuint iDim = 0; iDim < dim; ++iDim) {
//       m_pGradTimesVolDim[iDim] += m_pData[is][EulerTerm::P]*normalsData.
// 	getNodalNormComp(is,iDim); 
//     }
    
//     m_avSpeed[XX] += m_pData[is][EulerTerm::VX];
//     m_avSpeed[YY] += m_pData[is][EulerTerm::VY];
//     if (dim == DIM_3D) {
//       m_avSpeed[ZZ] += m_pData[is][EulerTerm::VZ];
//     }
    
//     const CFreal mach = m_pData[is][EulerTerm::V]/m_pData[is][EulerTerm::A];
//     if (mach < 1.0) {
//       smaller1.push_back(is);
//     }
//     else {
//       bigger1.push_back(is);
//     }    
//   }
//   m_avSpeed /= nbStates;
  
//   const bool applyFix = (bigger1.size() > 0 && smaller1.size() > 0) ? true : false;
  
//   if (applyFix) {
//     /*    CFreal maxNormalGradP = 0.0;
// 	CFuint nMaxGradPID = 0;
// 	for (CFuint is = 0; is < nbStates; ++is) {
// 	CFreal gradP = 0.0;
// 	for (CFuint iDim = 0; iDim < dim; ++iDim) {
// 	gradP += m_pGradTimesVolDim[iDim]*normalsData.getNodalUnitNormComp(is,iDim);
// 	}
// 	// take the absolute value of the pressure gradient
// 	gradP = std::abs(gradP);
// 	if (gradP > maxNormalGradP) {
// 	maxNormalGradP = gradP;
// 	nMaxGradPID = is;
// 	}
// 	}
	
// 	// n_shock corresponds to the unit normal along which
// 	// the maximum pressure gradient is detected
// 	for (CFuint iDim = 0; iDim < dim; ++iDim) {
// 	m_nShock[iDim] = normalsData.getNodalUnitNormComp(nMaxGradPID,iDim);
// 	}
	
// 	// the normal shock direction has to be downwind with respect 
// 	// to the cell average velocity
// 	if (MathFunctions::innerProd(m_avSpeed, m_nShock) < 0.0) {
// 	m_nShock *= -1.0;
// 	}
//     */
    
//     for (CFuint is = 0; is < nbStates; ++is) {
//       const RealVector& pv = m_pData[is];
//       CFreal un =  pv[EulerTerm::VX]*normalsData.getNodalUnitNormComp(is,XX) +
// 	pv[EulerTerm::VY]*normalsData.getNodalUnitNormComp(is,YY);
//       if (dim == DIM_3D) {
// 	un += pv[EulerTerm::VZ]*normalsData.getNodalUnitNormComp(is,ZZ);
//       }
      
//       m_lambda[is] = un;
//       m_lambda[is][dim]   += pv[EulerTerm::A];
//       m_lambda[is][dim+1] -= pv[EulerTerm::A];
//     }
    
//     CFreal tmpDelta = 0.0;
//     const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
//     for (CFuint iState = 0; iState < nbStates; ++iState) {
//       for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
// 	CFreal lambdaUpDen = 0.0;
// 	CFreal lambdaUpNum = 0.0;
// 	CFreal lambdaDownDen = 0.0;
// 	CFreal lambdaDownNum = 0.0;
// 	for (CFuint is = 0; is < nbStates; ++is) {
// 	  CFreal kj = 0.0;
// 	  for (CFuint iDim = 0; iDim < dim; ++iDim) {
// 	    // should we use here the unit inward normal ?? 
// 	    kj += normalsData.getNodalUnitNormComp(iState,iDim)*
// 	      normalsData.getNodalNormComp(is,iDim);
// 	  }
// 	  const CFreal kmin = min(0.,kj);
// 	  const CFreal kplus = max(0.,kj);
// 	  lambdaUpDen += kmin;
// 	  lambdaUpNum += kmin*m_lambda[is][iEq];
// 	  lambdaDownDen += kplus;
// 	  lambdaDownNum += kplus*m_lambda[is][iEq];
//       }
// 	const CFreal lambdaUp = lambdaUpNum/lambdaUpDen;
// 	const CFreal lambdaDown = lambdaDownNum/lambdaDownDen;
// 	tmpDelta = max(tmpDelta, std::abs(lambdaDown - lambdaUp));
//       }
//     }
//     tmpDelta *= 0.5;
    
//     const CFreal avV = sqrt(m_avSpeed[XX]*m_avSpeed[XX] + 
// 			    m_avSpeed[YY]*m_avSpeed[YY]);
//     const CFreal cosD = m_avSpeed[XX]/avV;
//     const CFreal sinD = m_avSpeed[YY]/avV;
//     cf_assert(dim == DIM_2D);
    
//     for (CFuint is = 0; is < nbStates; ++is) {
//       // only 2D case for now
//       const CFreal nEtaScalarNormal = -normalsData.getNodalUnitNormComp(is,XX)*sinD + 
// 	normalsData.getNodalUnitNormComp(is,YY)*cosD;
//       delta[is] = std::abs(nEtaScalarNormal)*tmpDelta;
//     }
//   }
//   //  }
// }
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
