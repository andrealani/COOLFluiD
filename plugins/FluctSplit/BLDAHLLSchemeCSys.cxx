#include "MathTools/MatrixInverter.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FluctSplit/FluctSplitSystem.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/BLDAHLLSchemeCSys.hh"

// LW
#include "NavierStokes/EulerTerm.hh"
// end - LW
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<BLDAHLLSchemeCSys,			
FluctuationSplitData,
                        Splitter,
                        FluctSplitSystemModule>
aBLDAHLLSchemeCSysProvider("SysBLDAHLLC");


//////////////////////////////////////////////////////////////////////////////

void BLDAHLLSchemeCSys::defineConfigOptions(Config::OptionList& options)
{ 
  options.addConfigOption< CFuint, Config::DynamicOption<> >("FirstOrder","Run first order.");
}

//////////////////////////////////////////////////////////////////////////////

BLDAHLLSchemeCSys::BLDAHLLSchemeCSys(const std::string& name) :
  BSchemeBase<RDHLLSchemeCSys>(name),
  m_phiLDA(),
  m_phiLW(),
  m_phiLLxF(),
  m_sumKmin(),
  m_inv_sumAbsK()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
  
  m_firstOrder = 0;
  setParameter("FirstOrder",&m_firstOrder); 
}
      
//////////////////////////////////////////////////////////////////////////////

BLDAHLLSchemeCSys::~BLDAHLLSchemeCSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void BLDAHLLSchemeCSys::configure ( Config::ConfigArgs& args )
{
  BSchemeBase<RDHLLSchemeCSys>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BLDAHLLSchemeCSys::setup()
{
  CFAUTOTRACE;
  
  BSchemeBase<RDHLLSchemeCSys>::setup();
  
  m_phiLDA.resize(_nbEquations);
//   m_phiLW.resize(_nbEquations);
//   m_dissLLxF.resize(_nbEquations); 
//   m_phiLLxF.resize(_nbEquations);
//   m_sumKmin.resize(_nbEquations, _nbEquations);
//   m_sumAbsK.resize(_nbEquations, _nbEquations);
//   m_inv_sumAbsK.resize(_nbEquations, _nbEquations);
  _Uavg.resize(_nbEquations);  
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->resizePhysicalData(_pData);  
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  n_IJ.resize(dim);
  n_IK.resize(dim);
    
  F_I.resize(_nbEquations);
  G_I.resize(_nbEquations);
  F_J.resize(_nbEquations);
  G_J.resize(_nbEquations);
  F_K.resize(_nbEquations);
  G_K.resize(_nbEquations);
  Fn_IJ.resize(_nbEquations);
  Fn_IK.resize(_nbEquations);
}
      
//////////////////////////////////////////////////////////////////////////////

void BLDAHLLSchemeCSys::distribute(vector<RealVector>& residual)
{
//std::cout << "void BLDAHLLSchemeCSys::distribute(vector<RealVector>& residual)\n";
  
  DistributionData& distdata = getMethodData().getDistributionData();
  const bool isPerturb = distdata.isPerturb;
  
  if ((m_firstOrder == 0 && !m_firstOrderJacob) || 
      (m_firstOrder == 0 && m_firstOrderJacob && !isPerturb)) {
    
    const RealVector& phiT = distdata.phi;
    const vector<State*>& tStates = *distdata.tStates;
    const CFuint nbEqs = _nbEquations;
    const CFuint nbStates = _nbStatesInCell;
    
    _sumKplusU = (*_kPlus[0])*(*tStates[0]);
    _sumKplus  = *_kPlus[0];
// 	m_sumKmin  = *_kMin[0];
    for (CFuint iState = 1; iState < nbStates; ++iState) {
      _sumKplusU += (*_kPlus[iState])*(*tStates[iState]);
      _sumKplus  += *_kPlus[iState];
// 	  // LW
// 	  m_sumKmin  += *_kMin[iState];
// 	  // end - LW
    }
    
    _inverter->invert(_sumKplus, _invK);
    _sumKplusU -= phiT;
    _uInflow = _invK * _sumKplusU;
    m_uTemp = _invK*phiT;
    
    // computate the residual of the P scheme, and its sum
    m_sumPhiN = 0.0;
    for (CFuint iState = 0; iState < nbStates; ++iState) {
	  const CFuint I = iState       % nbStates;
	  const CFuint J = (iState + 1) % nbStates ;
	  const CFuint K = (iState + 2) % nbStates ;

	  m_phiN[iState] = computeDimSplittedFVContribution(I, J, K); //computeDimSplittedFVContribution(m_phiN[iState], I, J, K);
	  
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  	m_sumPhiN[iEq] += std::abs(m_phiN[iState][iEq]);
      }
    }    
    
    computeBlendingCoeff();
    
    if ( m_store_thetas ) storeThetas();
    
    if ( m_addExtraDiss ) addExtraDissipation(residual);
	
	////////////////////////////////////////////////////////////////	
	// cell wise blending function
	  const CFreal mach_number = getMach();

	  CFreal sigma = 10.*mach_number;// chosen to be ~ 1. until Ma ~ 0.2 => Ma^2 ~ 0.04 subsonic limit

	  const CFreal fraction_LDA =std::tanh(sigma);
	  const CFreal fraction_LW =  1. - fraction_LDA;	  
	  // end - cell wise blending function
	  
// 	  m_sumAbsK = _sumKplus - m_sumKmin;
// 	  _inverter->invert(m_sumAbsK, m_inv_sumAbsK);
	  // end - LW
	  ////////////////////////////////////////////////////////////////
    
    
//     // variant #1 computation of LDA residual and the blending
//     for (CFuint iState = 0; iState < nbStates; ++iState)
//     {	  
// 	  m_phiLDA = (*_kPlus[iState])*m_uTemp;
// 	  // LW
// 	  m_phiLW  = 0.75*( (*_kPlus[iState]) + (*_kMin[iState]) ) * m_inv_sumAbsK*phiT; // hard coded cell CFL number
// 	  // end - LW
// 	  
//   	  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
// 		m_phiLW[iEq] += 1./3.*phiT[iEq];
// 		
// 	  	residual[iState][iEq] = m_theta[iEq]*m_phiN[iState][iEq] + 
// 	  	(1. - m_theta[iEq]) * ( fraction_LDA*m_phiLDA[iEq] + fraction_LW*m_phiLW[iEq] );	     
// // 	     m_theta[iEq]*m_phiN[iState][iEq] + (1. - m_theta[iEq])*m_phiLDA[iEq];
// 	  }
// 
// 	  if (distdata.computeBetas)
//       {
// 	    (*distdata.currBetaMat)[iState] = (*_kPlus[iState])*_invK;
// 	  }
//     }
//   }
//   else {
//     assert(m_firstOrder == 1 || (m_firstOrderJacob && isPerturb));
//     BSchemeBase<RDHLLSchemeCSys>::distribute(residual);
//   }
	// variant #2 computation of LDA residual and the blending	
    for (CFuint iState = 0; iState < nbStates; ++iState)
    {  
	  m_phiLDA = (*_kPlus[iState])*m_uTemp;

  	  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  	residual[iState][iEq] = m_theta[iEq]*m_phiN[iState][iEq] + (1. - m_theta[iEq]) * m_phiLDA[iEq];
	  }

	  if (distdata.computeBetas)
      {
	    (*distdata.currBetaMat)[iState] = (*_kPlus[iState])*_invK;
		std::cout << "Computing betas...\n"; abort();
	  }
    }
  }
  else {
	std::cout << "Aborting from else branch\n"; abort();
    assert(m_firstOrder == 1 || (m_firstOrderJacob && isPerturb));
    BSchemeBase<RDHLLSchemeCSys>::distribute(residual);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void BLDAHLLSchemeCSys::distributePart(vector<RealVector>& residual)
{  
 //std::cout << "void BLDAHLLSchemeCSys::distributePart(vector<RealVector>& residual)\n";
 throw Common::NotImplementedException (FromHere(),"void BLDAHLLSchemeCSys::distributePart()");
//  DistributionData& distdata = getMethodData().getDistributionData();
//   const bool isPerturb = distdata.isPerturb;
//   
//   if ((m_firstOrder == 0 && !m_firstOrderJacob) || 
//       (m_firstOrder == 0 && m_firstOrderJacob && !isPerturb)) {
//     
//     const CFuint nbEqs = _nbEquations;
//     RealVector& phi = distdata.phi;
//     const vector<State*>& tStates = *distdata.tStates;
//     
//     _sumKplusU.slice(0, nbEqs) = (*_kPlus[0]) * tStates[0]->slice(_firstVarID, nbEqs);
//     _sumKplus  =*_kPlus[0];
//     for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
//       _sumKplusU.slice(0, nbEqs) +=
// 	(*_kPlus[iState]) *  tStates[iState]->slice(_firstVarID, nbEqs);
//       _sumKplus  += *_kPlus[iState];
//     }
//     
//     _inverter->invert(_sumKplus, _invK);
//     
//     _sumKplusU.slice(0, nbEqs) -= phi.slice(_firstVarID, nbEqs);
//     _uInflow = _invK*_sumKplusU;
//     m_uTemp.slice(0, nbEqs) = _invK * phi.slice(_firstVarID, nbEqs);
//     
//     m_sumPhiN = 0.0;
// 	for (CFuint iState = 0; iState < _nbStatesInCell; ++iState){
// 	  m_phiN[iState].slice(0, nbEqs) = (*_kPlus[iState]) * ( tStates[iState]->slice(_firstVarID, nbEqs) -  _uInflow.slice(0, nbEqs) );
// 	  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
// 		m_sumPhiN[iEq] += std::abs(m_phiN[iState][iEq]);
// 	  }
// 	}
//     
//     computeBlendingCoeff();
//     
//     if ( m_store_thetas ) storeThetas();
//     
//     for (CFuint iState = 0; iState < _nbStatesInCell; ++iState){
// 		m_phiLDA = (*_kPlus[iState])*m_uTemp;
//     
// 	  for (CFuint iEq = _firstVarID, jEq = 0; iEq < _lastVarID; ++iEq, ++jEq){
// 		  residual[iState][iEq] = 0.;abort();
// 		  // 	      m_theta[jEq]*m_phiN[iState][jEq] +  (1. - m_theta[jEq])*
// 	  }
// 	}
//   }
//   else {
//     assert(m_firstOrder == 1 || (m_firstOrderJacob && isPerturb));
//     BSchemeBase<RDHLLSchemeCSys>::distributePart(residual);
//   }
}

//////////////////////////////////////////////////////////////////////////////

void BLDAHLLSchemeCSys::computeBlendingCoeff()
{
  std::cout << "void BLDAHLLSchemeCSys::computeBlendingCoeff()\n";
  const RealVector& phi = getMethodData().getDistributionData().phi;
  const CFuint nbEqs = _nbEquations;
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
  {
    m_theta[iEq] = std::max ( std::abs(phi[iEq])/std::max(MathTools::MathConsts::CFrealEps(), m_sumPhiN[iEq]) , m_min_theta );
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal BLDAHLLSchemeCSys::getMach(){
  
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Physics::NavierStokes;

  const CFuint N0 = 0;
  const CFuint N1 = 1;
  const CFuint N2 = 2;

  const vector<State*>& tStates = *(this->getMethodData().getDistributionData().tStates);
  
  const CFuint nbStates = tStates.size();

  _Uavg = (*tStates[N0] + *tStates[N1] + *tStates[N2])/this->_nbStatesInCell;
  
  this->getMethodData().getLinearVar()->computePhysicalData(_Uavg, _pData);
  
  const CFreal a = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::A];
  const CFreal V = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::V];
  
  return V/a;  
}

//////////////////////////////////////////////////////////////////////////////

CFreal BLDAHLLSchemeCSys::getMaxLambda(){
  
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Physics::NavierStokes;

  const vector<State*>& tStates = *(this->getMethodData().getDistributionData().tStates);
  
  const CFuint nbStates = tStates.size();  
  
  this->getMethodData().getLinearVar()->computePhysicalData(*tStates[0], _pData);
   
  CFreal v_plus_a  =_pData[COOLFluiD::Physics::NavierStokes::EulerTerm::V] + _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::A];  
  CFreal lambdaMax = v_plus_a ;
   
  for (CFuint iState = 1; iState < nbStates; iState++){
	this->getMethodData().getLinearVar()->computePhysicalData(*tStates[iState], _pData);
	v_plus_a = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::V] + _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::A];
	
	lambdaMax = std::max( v_plus_a, lambdaMax);
  }
  
  return lambdaMax;  
}

//////////////////////////////////////////////////////////////////////////////

CFreal BLDAHLLSchemeCSys::getFactor(){
  
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Physics::NavierStokes;

  const vector<State*>& tStates = *(this->getMethodData().getDistributionData().tStates);
  
  const CFuint nbStates = tStates.size(); 

  RealVector avgState =  1./3.*(*tStates[0])  + 1./3.*(*tStates[1]) + 1./3.*(*tStates[2]); 
  
  this->getMethodData().getLinearVar()->computePhysicalData( avgState, _pData);
   
  CFreal V =_pData[COOLFluiD::Physics::NavierStokes::EulerTerm::V] ;  
  
  return V;  
}
//////////////////////////////////////////////////////////////////////////////

RealVector BLDAHLLSchemeCSys::computeDimSplittedFVContribution(const CFuint I, const CFuint J, const CFuint K ){
  //std::cout << "RealVector& BLDAHLLSchemeCSys::computeDimSplittedFVContribution(const CFuint I, const CFuint J, const CFuint K)\n";

  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Physics::NavierStokes;

  const CFuint XX = 0; const CFuint YY = 1;

  const CFuint rhoID  = 0;
  const CFuint rhoUID = 1;
  const CFuint rhoVID = 2;
  const CFuint rhoEID = 3;

  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  DataHandle< InwardNormalsData*> normals = this->socket_normals.getDataHandle();
  FluctSplit::DistributionData& ddata = getMethodData().getDistributionData();

  GeometricEntity *const cell = this->getMethodData().getDistributionData().cell;    
  const CFuint cellID = cell->getID();

  n_IJ[XX] = (1./6.) *( normals[cellID]->getNodalNormComp(J,XX) - normals[cellID]->getNodalNormComp(I,XX) );
  n_IJ[YY] = (1./6.) *( normals[cellID]->getNodalNormComp(J,YY) - normals[cellID]->getNodalNormComp(I,YY) );

  n_IK[XX] = (1./6.) *( normals[cellID]->getNodalNormComp(K,XX) - normals[cellID]->getNodalNormComp(I,XX) );
  n_IK[YY] = (1./6.) *( normals[cellID]->getNodalNormComp(K,YY) - normals[cellID]->getNodalNormComp(I,YY) );


  const vector<State*>& tStates = *(this->getMethodData().getDistributionData().tStates);
  
  const CFuint nbStates = tStates.size();  
  
  this->getMethodData().getLinearVar()->computePhysicalData(*tStates[I], _pData);
  const CFreal rho_I = _pData[EulerTerm::RHO];
  const CFreal u_I   = _pData[EulerTerm::VX];
  const CFreal v_I   = _pData[EulerTerm::VY];
  const CFreal p_I   = _pData[EulerTerm::P];
  const CFreal E_I   = _pData[EulerTerm::E];
  const CFreal V_I   = _pData[EulerTerm::V];
  const CFreal A_I   = _pData[EulerTerm::A];

  F_I = 0.;                               G_I = 0.;
  F_I[rhoID]  = rho_I*u_I;                G_I[rhoID]  = rho_I*v_I;
  F_I[rhoUID] = rho_I*u_I*u_I + p_I;      G_I[rhoUID] = rho_I*u_I*v_I;
  F_I[rhoVID] = rho_I*v_I*u_I;            G_I[rhoVID] = rho_I*v_I*v_I + p_I;
  F_I[rhoEID] = (rho_I*E_I + p_I)*u_I;    G_I[rhoEID] = (rho_I*E_I + p_I)*v_I;

  this->getMethodData().getLinearVar()->computePhysicalData(*tStates[J], _pData);
  const CFreal rho_J = _pData[EulerTerm::RHO];
  const CFreal u_J   = _pData[EulerTerm::VX];
  const CFreal v_J   = _pData[EulerTerm::VY];
  const CFreal p_J   = _pData[EulerTerm::P];
  const CFreal E_J   = _pData[EulerTerm::E];
  const CFreal V_J   = _pData[EulerTerm::V];
  const CFreal A_J   = _pData[EulerTerm::A];

  F_J = 0.;                            G_J = 0.;
  F_J[rhoID]  = rho_J*u_J;             G_J[rhoID]  = rho_J*v_J;
  F_J[rhoUID] = rho_J*u_J*u_J + p_J;   G_J[rhoUID] = rho_J*u_J*v_J ;
  F_J[rhoVID] = rho_J*v_J*u_J;         G_J[rhoVID] = rho_J*v_J*v_J + p_J;
  F_J[rhoEID] = (rho_J*E_J + p_J)*u_J; G_J[rhoEID] = (rho_J*E_J + p_J)*v_J;

  this->getMethodData().getLinearVar()->computePhysicalData(*tStates[K], _pData);
  const CFreal rho_K = _pData[EulerTerm::RHO];
  const CFreal u_K   = _pData[EulerTerm::VX];
  const CFreal v_K   = _pData[EulerTerm::VY];
  const CFreal p_K   = _pData[EulerTerm::P];
  const CFreal E_K   = _pData[EulerTerm::E];
  const CFreal V_K   = _pData[EulerTerm::V];
  const CFreal A_K   = _pData[EulerTerm::A];

  F_K = 0.;                              G_K = 0.;
  F_K[rhoID]  = rho_K*u_K;               G_K[rhoID]  = rho_K*v_K;
  F_K[rhoUID] = rho_K*u_K*u_K + p_K;     G_K[rhoUID] = rho_K*u_K*v_K;
  F_K[rhoVID] = rho_K*v_K*u_K;           G_K[rhoVID] = rho_K*v_K*v_K + p_K;
  F_K[rhoEID] = (rho_K*E_K + p_K)*u_K;   G_K[rhoEID] = (rho_K*E_K + p_K)*v_K;

///////////////////////////////////////////////////////////////////////////////
//  HLL 

 const CFreal vI_n_IJ = u_I*n_IJ[XX] + v_I*n_IJ[YY];
 const CFreal vJ_n_IJ = u_J*n_IJ[XX] + v_J*n_IJ[YY];

 const CFreal SLeft_IJ  = std::min(vI_n_IJ - A_I, vJ_n_IJ - A_J); 

 const CFreal SRight_IJ = std::max(vI_n_IJ + A_I, vJ_n_IJ + A_J);

 Fn_IJ = 0.;

 if ( 0. <= SLeft_IJ ){
   Fn_IJ = F_I*n_IJ[XX] +  G_I*n_IJ[YY];
 }
 else if ( SRight_IJ <= 0. ){
   Fn_IJ = F_J*n_IJ[XX] +  G_J*n_IJ[YY];
 }
 else if ( (SLeft_IJ < 0.) && (0. < SRight_IJ) ){
  Fn_IJ = SRight_IJ*(F_I*n_IJ[XX] + G_I*n_IJ[YY]) - SLeft_IJ*(F_J*n_IJ[XX] + G_J*n_IJ[YY]) 
         + SLeft_IJ*SRight_IJ*( (*tStates[J]) -  (*tStates[I]) );

  Fn_IJ /= SRight_IJ - SLeft_IJ;
  if (std::abs(SRight_IJ - SLeft_IJ) < 1.e-9)
	std::cout << "SRight_IJ - SLeft_IJ is too small: " << (SRight_IJ - SLeft_IJ) << "\n";
 }
 else{
   std::cout << "SLeft_IJ = " << SLeft_IJ << "\tSRight_IJ = " << SRight_IJ << "\n";
   std::cout << "IJ Problem!\n"; abort();
 }

  const CFreal vI_n_IK =  u_I*n_IK[XX] + v_I*n_IK[YY];
  const CFreal vK_n_IK =  u_K*n_IK[XX] + v_K*n_IK[YY];

  const CFreal SLeft_IK  = std::min(vI_n_IK - A_I, vK_n_IK - A_K);
  const CFreal SRight_IK = std::max(vI_n_IK + A_I, vK_n_IK + A_K);

  Fn_IK = 0.;
  if ( 0. <= SLeft_IK ){
	Fn_IK = F_I*n_IK[XX] +  G_I*n_IK[YY];
  }
  else if ( SRight_IK <= 0. ){
	Fn_IK = F_K*n_IK[XX] +  G_K*n_IK[YY];
  }
  else if ( (SLeft_IK < 0.) && (0. < SRight_IK) ){
	Fn_IK = SRight_IK*(F_I*n_IK[XX] + G_I*n_IK[YY]) - SLeft_IK*(F_K*n_IK[XX] + G_K*n_IK[YY]) 
			+ SLeft_IK*SRight_IK*( (*tStates[K]) -  (*tStates[I]) );

	Fn_IK /= SRight_IK - SLeft_IK;
  	if (std::abs(SRight_IK - SLeft_IK) < 1.e-9)
	  std::cout << "SRight_IK - SLeft_IK is too small: " << (SRight_IK - SLeft_IK) << "\n";

  }
 else{
   std::cout << "SLeft_IK = " << SLeft_IK << "\tSRight_IK = " << SRight_IK << "\n";
   std::cout << "IK Problem!\n"; abort();
 }


  RealVector phiP_nodal = Fn_IJ + Fn_IK;

// end -  HLL 
/////////////////////////////////////////////////////////////////////////////////

  return phiP_nodal;  
}

//////////////////////////////////////////////////////////////////////////////

void BLDAHLLSchemeCSys::addExtraDissipation(vector<RealVector>& residual)
{
 // std::cout << "\tnothing is done\n";
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

   } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
