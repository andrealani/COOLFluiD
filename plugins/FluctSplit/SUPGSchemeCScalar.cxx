#include "FluctSplit/SUPGSchemeCScalar.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitScalar.hh"
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

MethodStrategyProvider<SUPGSchemeCScalar,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitScalarModule>
SUPGSchemeScalarCProvider("ScalarSUPGC");

//////////////////////////////////////////////////////////////////////////////

SUPGSchemeCScalar::SUPGSchemeCScalar(const std::string& name) :
  RDS_SplitterScalar(name),
  m_C1(0.),
  m_C2(0.),
  m_h(0.),
  m_tau(),
  m_lambda(),
  m_avgSpeed(),
  m_normGrad(),
  m_kappa_hat(),
  m_grad(),
  m_stateNormal(),
  m_p()
{
}

//////////////////////////////////////////////////////////////////////////////

SUPGSchemeCScalar::~SUPGSchemeCScalar()
{
}

//////////////////////////////////////////////////////////////////////////////

void SUPGSchemeCScalar::setup()
{
   RDS_SplitterScalar::setup();

   m_dim = Framework::PhysicalModelStack::getActive()->getDim();
   m_invDim = 1.0 / m_dim;

   m_C1 = 0.5;
   m_C2 = 0.5;

   m_tau.resize(_nbEquations);
   m_lambda.resize(_nbEquations);
   m_normGrad.resize(_nbEquations);
   m_kappa_hat.resize(_nbEquations);

   m_grad.resize(_nbEquations);
   m_avgSpeed.resize(_nbEquations);
   for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
      m_grad[iEq].resize(m_dim);
      m_avgSpeed[iEq].resize(m_dim);
   }

   m_stateNormal.resize(m_maxNbStatesInCell);
   m_p.resize(m_maxNbStatesInCell);
   for (CFuint iState = 0; iState < m_maxNbStatesInCell; ++iState) {
      m_stateNormal[iState].resize(m_dim);
      m_p[iState].resize(_nbEquations);
   }
}

//////////////////////////////////////////////////////////////////////////////

void SUPGSchemeCScalar::distribute(vector<RealVector>& residual)
{
  DistributionData& ddata = getMethodData().getDistributionData();
  const vector<State*>& tStates = *ddata.tStates;
  const RealVector& phiT = ddata.phi;
  
  const CFreal invnstates = 1.0 / _nbStatesInCell;

   // compute h
   m_h = m_normals->getAreaFace(0);
   for (CFuint iFace = 1; iFace < m_normals->nbFaces(); ++iFace) {
      m_h = std::max(m_normals->getAreaFace(iFace),m_h);
   }

CF_DEBUG_OBJ(m_h);

 m_invCellArea = 1.0 / ddata.cell->computeVolume();

CF_DEBUG_OBJ(m_invCellArea);

   for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
     
     // get average speed vector per variable
     getMethodData().getDistribVar()->getAverageAdvectionVector(m_avgSpeed[iEq],iEq);
     CF_DEBUG_OBJ(m_avgSpeed[iEq]);
     
     // compute lambda and tau
     m_lambda[iEq] = std::sqrt(m_avgSpeed[iEq].norm2()) + MathTools::MathConsts::CFrealEps();
     m_tau[iEq] = m_C1 * m_h / m_lambda[iEq];
     
     // reset the gradients
     m_grad[iEq] = 0.;
   }
   
   CF_DEBUG_OBJ(m_lambda);
   CF_DEBUG_OBJ(m_tau);

   // compute the gradient of the variables
   for (CFuint iDim = 0; iDim < m_dim; ++iDim) {
     for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
       
       // get nodal normal component
       m_stateNormal[iState][iDim] = m_normals->getNodalNormComp(iState,iDim);
       
       for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
	 m_grad[iEq][iDim] += (*tStates[iState])[iEq] * m_stateNormal[iState][iDim];
       }
     }
   }
   for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
     
      m_grad[iEq] *= m_invDim;
      m_grad[iEq] *= m_invCellArea;
      
      m_normGrad[iEq] = m_grad[iEq].norm2();
      
      // compute kappa hat
      m_kappa_hat[iEq] = m_C2 * m_h / (m_normGrad[iEq] + m_h);
      if( phiT[iEq] < 0.0 ) {
	m_kappa_hat[iEq] *= - 1.0;
      }
      CF_DEBUG_OBJ(m_grad[iEq]);
      CF_DEBUG_OBJ(m_normGrad[iEq]);
      CF_DEBUG_OBJ(m_kappa_hat[iEq]);
   }
   
   
   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
     CF_DEBUG_OBJ(m_stateNormal[iState]);
     CF_DEBUG_OBJ(*tStates[iState]);
     for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
       
       // compute p parameter
       m_p[iState][iEq] = m_invDim * MathFunctions::innerProd(m_grad[iEq],m_stateNormal[iState]);
       
       // distribute to nodes
       //          residual[iState][iEq] =
//             (invnstates + (m_tau[iEq]*_k[iState][iEq] + m_kappa_hat[iEq] * m_p[iState][iEq]) * m_invCellArea) * phiT[iEq] ;
       
       residual[iState][iEq] =
	 (invnstates + (m_tau[iEq]*_k[iState][iEq]) * m_invCellArea) * phiT[iEq] ;
     }
     CF_DEBUG_OBJ(m_p[iState]);
     CF_DEBUG_OBJ(residual[iState]);
   }
   CF_DEBUG_OBJ(phiT);
   CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////////////

void SUPGSchemeCScalar::distributePart(vector<RealVector>& residual)
{
  throw Common::NotImplementedException (FromHere(),"SUPGSchemeCScalar::distributePart()");
}

//////////////////////////////////////////////////////////////////////////////

void SUPGSchemeCScalar::computePicardJacob(vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"SUPGSchemeCScalar::computePicardJacob()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
