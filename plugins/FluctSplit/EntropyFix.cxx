#include "Framework/GeometricEntity.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"

#include "FluctSplit/EntropyFix.hh"
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

MethodStrategyProvider<EntropyFix,
                       FluctuationSplitData,
                       ComputeJacobianFix,
                       FluctSplitNavierStokesModule>
entropyFixProvider("Entropy");

//////////////////////////////////////////////////////////////////////////////

EntropyFix::EntropyFix(const std::string& name) :
  ComputeJacobianFix(name),
  m_pGradTimesVolDim(),
  m_nShock(),
  m_pData(),
  m_lambda(),
  m_avSpeed()
{
}

//////////////////////////////////////////////////////////////////////////////

EntropyFix::~EntropyFix()
{
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFix::setup()
{
  ComputeJacobianFix::setup();
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
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
  m_avSpeed.resize(dim);
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFix::computeFix(const InwardNormalsData& normalsData,
			    RealVector& delta)
{ 
  // here we are assuming that the eigenvalues of the K matrices 
  // related to update variables are the same as for distribution
  // it is not true in case characteristic variables for HE splitting
  // are used to distribute the residual
  DistributionData& ddata = getMethodData().getDistributionData();
  if (!ddata.isPerturb) {
    const CFuint nbStates = ddata.states->size();
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    
    m_pGradTimesVolDim = 0.0;
    m_avSpeed          = 0.0;
    for (CFuint is = 0; is < nbStates; ++is) {
      getMethodData().getUpdateVar()->computePhysicalData(*(*ddata.states)[is], m_pData[is]);
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
        m_pGradTimesVolDim[iDim] += m_pData[is][EulerTerm::P]*normalsData.getNodalNormComp(is,iDim); 
      }
      
      m_avSpeed[XX] += m_pData[is][EulerTerm::VX];
      m_avSpeed[YY] += m_pData[is][EulerTerm::VY];
      if (dim == DIM_3D) {
        m_avSpeed[ZZ] += m_pData[is][EulerTerm::VZ];
      }
    }
    m_avSpeed /= nbStates;
    
    CFreal maxNormalDeltaP = 0.0;
    CFuint nMaxGradPID     = 0;
    for (CFuint is = 0; is < nbStates; ++is) {
      CFreal deltaP = 0.0;
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
        deltaP += m_pGradTimesVolDim[iDim]*normalsData.getNodalUnitNormComp(is,iDim);
      }
      // take the absolute value of the pressure gradient
      deltaP = std::abs(deltaP);
      
      if (deltaP > maxNormalDeltaP) {
        maxNormalDeltaP = deltaP;
        nMaxGradPID    = is;
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

	// JGM: Do not remove this assertion till extensive experience with the entropy fix has been collected.
    assert( MathFunctions::innerProd(m_pGradTimesVolDim, m_nShock) < 0. );

    /*  for (CFuint is = 0; is < nbStates; ++is) {
	const RealVector& pv = m_pData[is];
	m_lambda[is] = pv[EulerTerm::VX]*m_nShock[XX] + 
	pv[EulerTerm::VY]*m_nShock[YY] - pv[EulerTerm::A];
	if (dim == DIM_3D) { 
	m_lambda[is] += pv[EulerTerm::VZ]*m_nShock[ZZ];
	}
	}
    */
    
    for (CFuint is = 0; is < nbStates; ++is) {
      const RealVector& pv = m_pData[is];
      m_lambda[is] = pv[EulerTerm::VX]*m_nShock[XX] + pv[EulerTerm::VY]*m_nShock[YY] - pv[EulerTerm::A];

      if (dim == DIM_3D) {
        m_lambda[is] += pv[EulerTerm::VZ]*m_nShock[ZZ];
      }
    }
    
    CFreal lambdaUpDen   = 0.0;
    CFreal lambdaUpNum   = 0.0;
    CFreal lambdaDownDen = 0.0;
    CFreal lambdaDownNum = 0.0;
    
    for (CFuint is = 0; is < nbStates; ++is) {
      const CFreal proj = -m_pGradTimesVolDim[0]*m_nShock[0] -m_pGradTimesVolDim[1]*m_nShock[1];assert(dim == 2);
      const CFuint sig  = proj >= 0 ? proj > 0 : -1;

      CFreal kj = normalsData.getNodalUnitNormComp(is, XX)*m_nShock[XX] + normalsData.getNodalUnitNormComp(is, YY)*m_nShock[YY];
      kj *= sig;
      
      const CFreal kmin  = std::min(0.,kj);
      const CFreal kplus = std::max(0.,kj);
      
      lambdaUpDen += kmin;
      lambdaUpNum += kmin*m_lambda[is];
      
      lambdaDownDen += kplus;
      lambdaDownNum += kplus*m_lambda[is];
    }
    
    const CFreal lambdaUp   = lambdaUpNum/lambdaUpDen;
    const CFreal lambdaDown = lambdaDownNum/lambdaDownDen;

    delta = 0.;
     for (CFuint is = 0; is < nbStates; ++is) {
       if (is == nMaxGradPID)
         delta[is] = 0.5*std::max(0., lambdaDown - lambdaUp);
     }
//     delta = 0.5*std::max(0., lambdaDown - lambdaUp);
// //    delta = std::max(0., lambdaDown - lambdaUp);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
} // namespace COOLFluiD
