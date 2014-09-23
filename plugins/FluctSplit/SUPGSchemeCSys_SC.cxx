#include "SUPGSchemeCSys_SC.hh"
#include "MathTools/MatrixInverter.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitSystem.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/MeshData.hh"

#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<SUPGSchemeCSys_SC,
		       FluctuationSplitData,
		       Splitter,
               FluctSplitSystemModule>
supg_scSchemeCSysProvider("SysSUPG_SC");

//////////////////////////////////////////////////////////////////////////////

void SUPGSchemeCSys_SC::defineConfigOptions(Config::OptionList& options)
{ 
  options.addConfigOption< CFuint, Config::DynamicOption<> >("FirstOrder","Run first order.");
}

//////////////////////////////////////////////////////////////////////////////

SUPGSchemeCSys_SC::SUPGSchemeCSys_SC(const std::string& name) :
  BSchemeBase<NSchemeCSys>(name)
  // _sumKplusU(),
 //  _sumKplus(),
//   _invK(),
//   _uInflow(),
//   _uDiff(),
//   _temp(),
//   _tempBkp(),
//   _tempMat(),
//   _tmp(),
//   _sumKU(),
//   _stab_Adv(),
//   _stab_Shock(),
//   _invTau(),
//   _Tau(),
//   _Ax(),
//   _Ay(),
//   _Uavg(),
//   _dUHdU(),
//   m_phiLDA()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);
  
  m_firstOrder = 0;
  setParameter("FirstOrder",&m_firstOrder); 
}
      
//////////////////////////////////////////////////////////////////////////////

SUPGSchemeCSys_SC::~SUPGSchemeCSys_SC()
{
}

//////////////////////////////////////////////////////////////////////////////

void SUPGSchemeCSys_SC::configure ( Config::ConfigArgs& args )
{
  BSchemeBase<NSchemeCSys>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

void SUPGSchemeCSys_SC::setup()
{

  CFAUTOTRACE;
  
  BSchemeBase<NSchemeCSys>::setup();
  
  m_phiLDA.resize(_nbEquations);
  
  _sumKplusU.resize(_nbEquations);
  _sumKplus.resize(_nbEquations, _nbEquations);
  _invK.resize(_nbEquations, _nbEquations);
  _uInflow.resize(_nbEquations);
  _uDiff.resize(_nbEquations);
  _temp.resize(_nbEquations);
  _tempBkp.resize(_nbEquations);
  _tempMat.resize(_nbEquations, _nbEquations);
  _tmp.resize(_nbEquations,_nbEquations);
  _sumKU.resize(_nbEquations);

  // Upwind stabilization:
  const CFuint dim      = PhysicalModelStack::getActive()->getDim();
  // Valid only for P1 elements
  const CFuint nbStates = dim + 1;
  
  _stab_Adv.resize(_nbEquations);
  _stab_Shock.resize(_nbEquations);
  _invTau.resize(_nbEquations, _nbEquations);
  _Tau.resize(_nbEquations, _nbEquations);
  _Ax.resize(_nbEquations, _nbEquations);
  _Ay.resize(_nbEquations, _nbEquations);
  
  _speedNormal.resize( dim );
  _Uavg.resize(_nbEquations);
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->resizePhysicalData(_pData);
  
  // end -- Upwind stabilization
  
  // -- Shock capturing term:  
  _gradXi.resize(dim);  
  _gradEta.resize(dim);
  
  _dUdx.resize(_nbEquations);
  _dUdy.resize(_nbEquations);
  _dUdz.resize(_nbEquations);
  
  _dVdU.resize(_nbEquations, _nbEquations);
  _dUdV.resize(_nbEquations, _nbEquations);
  _M.resize(dim, nbStates);

  _dUHdU.resize(_nbEquations, _nbEquations);
  
  // end -- Shock capturing term  
}
      
//////////////////////////////////////////////////////////////////////////////

void SUPGSchemeCSys_SC::distribute(vector<RealVector>& residual)
{

  const CFuint rhoID  = 0;
  const CFuint rhoUID = 1;
  const CFuint rhoVID = 2;
  const CFuint rhoEID = 3; const CFuint rhoHID = rhoEID;

  const CFuint N0 = 0; const CFuint N1 = 1; const CFuint N2 = 2;

  DistributionData& ddata = getMethodData().getDistributionData();
  
  const CFuint cellID = ddata.cellID;
  
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;

  _Uavg = 0.;
  
  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
    _Uavg += *tStates[iState];
  }

  _Uavg /= static_cast<CFreal>(_nbStatesInCell);

  this->getMethodData().getLinearVar()->computePhysicalData( _Uavg, _pData);

  const CFreal rhoAvg = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::RHO];
  const CFreal VAvg   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::V];
  const CFreal aAvg   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::A];
  const CFreal TAvg   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::T];

  // Element length-scale
  CFreal h = 0.;

  _speedNormal[XX] = -1.;
  _speedNormal[YY] = 0.;

  if ( 1.e-5 < VAvg ){
    _speedNormal[XX] = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::VX]/VAvg;
    _speedNormal[YY] = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::VY]/VAvg;
  }


  DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<CFreal> volumes             = socket_volumes.getDataHandle();

  const CFreal vol = volumes[cellID];

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

    CFreal nIx = normals[cellID]->getNodalNormComp(iState, XX);
    CFreal nIy = normals[cellID]->getNodalNormComp(iState, YY);
    
    const CFreal unitU_DOT_nI = (_speedNormal[XX]*nIx + _speedNormal[YY]*nIy)/(2.*vol);
    h += std::abs(unitU_DOT_nI) ;
  }

  h = 2. / h;

  // Anderson's Hypersonics book (new edition), p 292
  const CFreal muRef = 1.789e-5;
  const CFreal TRef  = 288.;
  const CFreal S     = 110.;

  const CFreal muAvg = muRef*std::pow(TAvg/TRef, 1.5) * (TRef + S)/(TAvg + S);

  const CFreal nuAvg = muAvg / rhoAvg;

  const CFreal Pr = 0.71;
  const CFreal alphaAvg = nuAvg / Pr;

  // Tau matrix
  const CFreal inv_tauCont   =  2. * (VAvg + aAvg) / h;
  const CFreal tauCont   = 1./inv_tauCont;
  const CFreal tauMom    = std::pow( std::pow(inv_tauCont, 2) + std::pow(4.*nuAvg   /(h*h) , 2), -0.5 );
  const CFreal tauEnergy = std::pow( std::pow(inv_tauCont, 2) + std::pow(4.*alphaAvg/(h*h) , 2), -0.5 );

  _Tau = 0.;

//   _Tau(rhoID, rhoID)   = tauCont;
//   _Tau(rhoUID, rhoUID) = _Tau(rhoVID, rhoVID) = tauMom;
//   _Tau(rhoEID, rhoEID) = tauEnergy;

   // Kirk et al., AIAA 2013, Eq. 78
  _invTau = (*_kPlus[N0])  - (*_kMin[N0]) + (*_kPlus[N1])  - (*_kMin[N1]) + (*_kPlus[N2])  - (*_kMin[N2]);
  _inverter->invert(_invTau, _Tau); 

  _Tau *= vol;

   // end - AIAA 2013

  const CFreal one_OVER_dimPlusOne = 1. / static_cast<CFreal>(_nbStatesInCell);

//   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
//   // residual[iState] = one_OVER_dimPlusOne * phiT + (1./vol) * _Tau *  ( (*_kPlus[iState]) + (*_kMin[iState]) ) * phiT; // Kirk ordering
//     residual[iState] = one_OVER_dimPlusOne * phiT + (1./vol) * ( (*_kPlus[iState]) + (*_kMin[iState]) ) * _Tau *  phiT; // vdW ordering
// 
//   }

  // -----------------------------------------------------------------------------------
  
  const bool isPerturb = ddata.isPerturb;
  
  if ((m_firstOrder == 0 && !m_firstOrderJacob) ||
      (m_firstOrder == 0 && m_firstOrderJacob && !isPerturb)) {

    /////////////////////////////////////////////
    const CFreal cellVolume = getMethodData().getDistributionData().cell->computeVolume();

    const CFreal PI    = MathTools::MathConsts::CFrealPi();
    //    const CFreal h     = std::sqrt(4*cellVolume/PI);

    const CFreal ALPHA = (VAvg + aAvg)*h;
    /////////////////////////////////////////////


    // ********************************************************************** //
    // Dissipative part of LxF scheme, and its sum
    // store it in m_sumPhiN
    m_sumPhiN = 0.0;
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {

      const CFuint I = iState       % _nbStatesInCell;
      const CFuint J = (iState + 1) % _nbStatesInCell ;
      const CFuint K = (iState + 2) % _nbStatesInCell ;

      m_phiN[iState] = ALPHA*( *tStates[I] - *tStates[J] + *tStates[I] - *tStates[K] ) / static_cast<CFreal>(_nbStatesInCell);

      for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
        m_sumPhiN[iEq] += std::abs(m_phiN[iState][iEq]);
      }
    }


    /// Compute cell gradient, and ///
    /// Compute cell averaged state from nodal (conservative) variables ///

    const CFreal coeffGrad = 1./(2.*vol);
    _dUdx = 0.;
    _dUdy = 0.;
    _dUdz = 0.;

    for (CFuint kStates = 0; kStates < _nbStatesInCell; ++kStates){
      const CFreal nK_x = normals[cellID]->getNodalNormComp(kStates,XX);
      const CFreal nK_y = normals[cellID]->getNodalNormComp(kStates,YY);

      _dUdx += nK_x * (*tStates[kStates]);
      _dUdy += nK_y * (*tStates[kStates]);

    }
    _dUdx *= coeffGrad;
    _dUdy *= coeffGrad;

    /// end - Compute cell gradient and average conservative state///
  // ********************************************************************** //
  
  
  // ********************************************************************** //
  // * Copied from SUPG_ShockCapturingStrategy.cxx ************************ //
    /// Compute gradients of xi, eta (and dseta) variables - Rao book on FEM, eq. 3.76 and eq. 3.75///
  /// Eq. 3.75 : dfdx_j = sum_{k=1}^{3} { dLkdxj * dfdLk}

  /// Eq. 3.76 for a P1 triangle ///

  GeometricEntity *const m_cell = m_cellBuilder.buildGE();

  const Node& node0 = *m_cell->getNode(N0);
  const Node& node1 = *m_cell->getNode(N1);
  const Node& node2 = *m_cell->getNode(N2);
  
  const CFreal x0 = node0[XX]; const CFreal y0 = node0[YY];
  const CFreal x1 = node1[XX]; const CFreal y1 = node1[YY];
  const CFreal x2 = node2[XX]; const CFreal y2 = node2[YY];


  /***************************************************************************/

  // Modify _dUdxj, to have drhoHdxj - IJNMF 2010 Kirk
  // Matrix dU^HdU is needed, where U^H = [rho rhoU rhoV rhoH ]'
  // rhoH = rhoE + p = rhoE;
  // p = Rg*rho*T = (Rg/cv)*rho*e = gamma_MINUS_one*(rhoE - 0.5*(rhoU**2 + rhoV**2)*rho**-1 );
  // Finally rhoH = gamma*rhoE - 0.5*gamma_MINUS_one*(rhoU**2 + rhoV**2)*rho**-1;
	
  const CFreal gamma    = 1.4;
  const CFreal gammaBar = gamma - 1.;
  const CFreal u   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::VX];
  const CFreal v   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::VY];

  _dUHdU = 0.;
  
  _dUHdU(rhoID , rhoID ) = 1.;

  _dUHdU(rhoUID, rhoUID) = 1.;

  _dUHdU(rhoVID, rhoVID) = 1.;

  _dUHdU(rhoHID, rhoID ) = 0.5 * gammaBar * (u*u + v*v);
  _dUHdU(rhoHID, rhoUID) = -gammaBar * u;
  _dUHdU(rhoHID, rhoVID) = -gammaBar * v;
  _dUHdU(rhoHID, rhoHID) = gamma;
 // end - IJNMF 2010 Kirk
    
    computeBlendingCoeff();
    
    if ( m_store_thetas ) storeThetas();
    
    if ( m_addExtraDiss ) addExtraDissipation(residual);
    
    
    // computation of LDA residual and the blending
    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
    {
      _stab_Adv   = (1./vol) * ( (*_kPlus[iState]) + (*_kMin[iState]) ) * _Tau *  phiT;

      const CFreal nK_x = normals[cellID]->getNodalNormComp(iState, XX);
      const CFreal nK_y = normals[cellID]->getNodalNormComp(iState, YY);
      
      _stab_Shock = ALPHA*(0.5*nK_x*_dUHdU*_dUdx + 0.5*nK_y*_dUHdU*_dUdy);  // WORKS
// ALPHA*(0.5*nK_x*_dUdx + 0.5*nK_y*_dUdy);  // WORKS
 
  	  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
        residual[iState][iEq] = one_OVER_dimPlusOne*phiT[iEq] + _stab_Adv[iEq] + m_theta[iEq]*_stab_Shock[iEq];
      }

      if (ddata.computeBetas)
      {
        (*ddata.currBetaMat)[iState] = (1./vol) * ( (*_kPlus[iState]) + (*_kMin[iState]) ) * _Tau; // vdW, still missing 1./3. in the diagonal
        for(CFuint iEq = 0; iEq < _nbEquations; ++iEq)
          (*ddata.currBetaMat)[iState](iEq, iEq) += one_OVER_dimPlusOne;
      }
    }
  }
  else {
	std::cout << "You shouldn't be here! --- aborting\n";abort();
    assert(m_firstOrder == 1 || (m_firstOrderJacob && isPerturb));
    BSchemeBase<NSchemeCSys>::distribute(residual);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void SUPGSchemeCSys_SC::distributePart(vector<RealVector>& residual)
{
//   const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
//   const RealVector& phiT = getMethodData().getDistributionData().phi;
//
//   _sumKplusU.slice(0, _nbEquations) = (*_kPlus[0]) * tStates[0]->slice(_firstVarID, _nbEquations);
//   _sumKplus  = *_kPlus[0];
//   for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
//     _sumKplusU.slice(0, _nbEquations) +=
//       (*_kPlus[iState]) * tStates[iState]->slice(_firstVarID, _nbEquations);
//     _sumKplus  += *_kPlus[iState];
//   }
//
//   _inverter->invert(_sumKplus, _invK);
//
//   RealVector& phi = const_cast<RealVector&>(phiT);
//
//   _sumKplusU.slice(0, _nbEquations) -= phi.slice(_firstVarID, _nbEquations);
//   _uInflow = _invK * _sumKplusU;
//
//   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
//     residual[iState].slice(_firstVarID, _nbEquations) = (*_kPlus[iState]) *
//       (tStates[iState]->slice(_firstVarID, _nbEquations) - _uInflow.slice(0, _nbEquations));
//   }
}

//////////////////////////////////////////////////////////////////////////////

void SUPGSchemeCSys_SC::computeBlendingCoeff()
{
  const RealVector& phi = getMethodData().getDistributionData().phi;
  const CFuint nbEqs = _nbEquations;
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
  {
    m_theta[iEq] = max ( std::abs(phi[iEq])/max(MathTools::MathConsts::CFrealEps(), m_sumPhiN[iEq]) , m_min_theta );
  }
}


//////////////////////////////////////////////////////////////////////////////

void SUPGSchemeCSys_SC::addExtraDissipation(vector<RealVector>& residual)
{
 // std::cout << "\tnothing is done\n";
}

//////////////////////////////////////////////////////////////////////////////

void SUPGSchemeCSys_SC::computePicardJacob(vector<RealMatrix*>& jacob)
{
   std::cout << "SUPGSchemeCSys_SC::computePicardJacob(vector<RealMatrix*>& jacob)\n";
   std::cout << "Aborting!\n"; abort();
//   // careful with the signs !!!
// //   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
// //     const CFuint nStart = iState*_nbStatesInCell;
// //     for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
// //       RealMatrix *const block = jacob[nStart + jState];
// //
// //       _tempMat = _invK*(*_kMin[jState]);
// //       if (iState == jState) {
// //  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
// //    _tempMat(iEq,iEq) -= 1.0;
// //  }
// //       }
// //
// //       _tempMat *= -1.0;
// //
// //       (*block) = (*_kPlus[iState])*_tempMat;
// //     }
// //   }
}

//////////////////////////////////////////////////////////////////////////////

CFreal SUPGSchemeCSys_SC::getNormSq(const RealVector vectorW){
  
  CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  assert(vectorW.size() == nbEqs);

//   CFreal normSqBydQdU = 0.;
//   for (CFuint iEq = 0; iEq < nbEqs; iEq++){
//       normSqBydQdU += vectorW[iEq]*vectorW[iEq];
//   }
// 
//   assert(normSqBydQdU >= 0.);

  RealVector intermediateVector(nbEqs);
//   intermediateVector = _dVdU*vectorW;

  for (CFuint iEq = 0; iEq < nbEqs; iEq++){
    intermediateVector[iEq] = 0.;
    for (CFuint jEq = 0; jEq < nbEqs; jEq++){
      intermediateVector[iEq] += _dVdU(iEq,jEq) * vectorW[jEq];
//       intermediateVector[iEq] += 1.*vectorW[jEq];
    }
  }
  CFreal normSqBydQdU = 0.;

  for (CFuint iEq = 0; iEq < nbEqs; iEq++){
      normSqBydQdU += vectorW[iEq] * intermediateVector[iEq];
  }

  return normSqBydQdU;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////



  // Advective Jacobian matrices

//  const CFreal Rg    = 8314.4/28.84;
//  const CFreal gamma = 1.4;
//  const CFreal gamma_MINUS_one = gamma - 1.;
//
//  const CFreal cv    = Rg/gamma_MINUS_one;
//
//  const CFreal VxAvg   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::VX];
//  const CFreal VyAvg   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::VY];
//  const CFreal H       = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::H];
//
//  _Ax = 0.;
//
//  _Ax(rhoID, rhoUID) = 1.0;
//
//  _Ax(rhoUID, rhoID)  = 0.5*gamma_MINUS_one * VAvg * VAvg - VxAvg*VxAvg ;
//  _Ax(rhoUID, rhoUID) = (3. - gamma ) * VxAvg;
//  _Ax(rhoUID, rhoVID) = -gamma_MINUS_one * VyAvg;
//  _Ax(rhoUID, rhoEID) = gamma_MINUS_one;
//
//  _Ax(rhoVID, rhoID)  = -VxAvg*VyAvg;
//  _Ax(rhoVID, rhoUID) = VyAvg;
//  _Ax(rhoVID, rhoVID) = VxAvg;
////   _Ax(rhoVID, rhoEID) = 0.;
//
//  _Ax(rhoEID, rhoID)  = VxAvg * ( 0.5*gamma_MINUS_one*VAvg*VAvg - H ) ;
//  _Ax(rhoEID, rhoUID) = H - gamma_MINUS_one*VxAvg*VxAvg;
//  _Ax(rhoEID, rhoVID) = -gamma_MINUS_one*VxAvg*VyAvg;
//  _Ax(rhoEID, rhoEID) = gamma*VxAvg;
//
//  _Ay = 0.;
//
//  _Ay(rhoID, rhoVID) = 1.0;
//
//  _Ay(rhoUID, rhoID)  = -VxAvg*VyAvg;
//  _Ay(rhoUID, rhoUID) = VyAvg;
//  _Ay(rhoUID, rhoVID) = VxAvg;
////   _Ay(rhoUID, rhoEID) = 0.;
//
//  _Ay(rhoVID, rhoID)  = 0.5*gamma_MINUS_one*VAvg*VAvg - VyAvg*VyAvg ;
//  _Ay(rhoVID, rhoUID) = -gamma_MINUS_one*VxAvg;
//  _Ay(rhoVID, rhoVID) = (3. - gamma) * VyAvg;
//  _Ay(rhoVID, rhoEID) = gamma_MINUS_one;
//
//  _Ay(rhoEID, rhoID)  = VyAvg* ( 0.5*gamma_MINUS_one*VAvg*VAvg - H );
//  _Ay(rhoEID, rhoUID) = -gamma_MINUS_one*VxAvg*VyAvg;
//  _Ay(rhoEID, rhoVID) = H - gamma_MINUS_one*VyAvg*VyAvg ;
//  _Ay(rhoEID, rhoEID) = gamma*VyAvg;

// ****************************************************************************
// ****************************************************************************
// ****************************************************************************
// ****************************************************************************
// ****************************************************************************
// ****************************************************************************
// ****************************************************************************
// ****************************************************************************



// method distribute, version 2014/02/06
//void SUPGSchemeCSys_SC::distribute(vector<RealVector>& residual)
// {
// 
//   const CFuint rhoID  = 0;
//   const CFuint rhoUID = 1;
//   const CFuint rhoVID = 2;
//   const CFuint rhoEID = 3; const CFuint rhoHID = rhoEID;
// 
//   const CFuint N0 = 0; const CFuint N1 = 1; const CFuint N2 = 2;
// 
//   DistributionData& ddata = getMethodData().getDistributionData();
//   
//   const CFuint cellID = ddata.cellID;
//   
//   const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
//   const RealVector& phiT = getMethodData().getDistributionData().phi;
// 
//   _Uavg = 0.;
//   
//   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
//     _Uavg += *tStates[iState];
//   }
// 
//   _Uavg /= static_cast<CFreal>(_nbStatesInCell);
// 
//   this->getMethodData().getLinearVar()->computePhysicalData( _Uavg, _pData);
// 
//   const CFreal rhoAvg = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::RHO];
//   const CFreal VAvg   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::V];
//   const CFreal aAvg   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::A];
//   const CFreal TAvg   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::T];
// 
//   // Element length-scale
//   CFreal h = 0.;
// 
//   _speedNormal[XX] = -1.;
//   _speedNormal[YY] = 0.;
// 
//   if ( 1.e-5 < VAvg ){
//     _speedNormal[XX] = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::VX]/VAvg;
//     _speedNormal[YY] = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::VY]/VAvg;
//   }
// 
// 
//   DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();
//   DataHandle<CFreal> volumes             = socket_volumes.getDataHandle();
// 
//   const CFreal vol = volumes[cellID];
// 
//   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
// 
//     CFreal nIx = normals[cellID]->getNodalNormComp(iState, XX);
//     CFreal nIy = normals[cellID]->getNodalNormComp(iState, YY);
//     
//     const CFreal unitU_DOT_nI = (_speedNormal[XX]*nIx + _speedNormal[YY]*nIy)/(2.*vol);
//     h += std::abs(unitU_DOT_nI) ;
//   }
// 
//   h = 2. / h;
// 
//   // Anderson's Hypersonics book (new edition), p 292
//   const CFreal muRef = 1.789e-5;
//   const CFreal TRef  = 288.;
//   const CFreal S     = 110.;
// 
//   const CFreal muAvg = muRef*std::pow(TAvg/TRef, 1.5) * (TRef + S)/(TAvg + S);
// 
//   const CFreal nuAvg = muAvg / rhoAvg;
// 
//   const CFreal Pr = 0.71;
//   const CFreal alphaAvg = nuAvg / Pr;
// 
//   // Tau matrix
//   const CFreal inv_tauCont   =  2. * (VAvg + aAvg) / h;
//   const CFreal tauCont   = 1./inv_tauCont;
//   const CFreal tauMom    = std::pow( std::pow(inv_tauCont, 2) + std::pow(4.*nuAvg   /(h*h) , 2), -0.5 );
//   const CFreal tauEnergy = std::pow( std::pow(inv_tauCont, 2) + std::pow(4.*alphaAvg/(h*h) , 2), -0.5 );
// 
//   _Tau = 0.;
// 
// //   _Tau(rhoID, rhoID)   = tauCont;
// //   _Tau(rhoUID, rhoUID) = _Tau(rhoVID, rhoVID) = tauMom;
// //   _Tau(rhoEID, rhoEID) = tauEnergy;
// 
//    // Kirk et al., AIAA 2013, Eq. 78
//   _invTau = (*_kPlus[N0])  - (*_kMin[N0]) + (*_kPlus[N1])  - (*_kMin[N1]) + (*_kPlus[N2])  - (*_kMin[N2]);
//   _inverter->invert(_invTau, _Tau); 
// 
//   _Tau *= vol;
// 
//    // end - AIAA 2013
// 
//   const CFreal one_OVER_dimPlusOne = 1. / static_cast<CFreal>(_nbStatesInCell);
// 
// //   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
// //   // residual[iState] = one_OVER_dimPlusOne * phiT + (1./vol) * _Tau *  ( (*_kPlus[iState]) + (*_kMin[iState]) ) * phiT; // Kirk ordering
// //     residual[iState] = one_OVER_dimPlusOne * phiT + (1./vol) * ( (*_kPlus[iState]) + (*_kMin[iState]) ) * _Tau *  phiT; // vdW ordering
// // 
// //   }
// 
//   // -----------------------------------------------------------------------------------
//   
//   const bool isPerturb = ddata.isPerturb;
//   
//   if ((m_firstOrder == 0 && !m_firstOrderJacob) ||
//       (m_firstOrder == 0 && m_firstOrderJacob && !isPerturb)) {
// 
//     /////////////////////////////////////////////
//     const CFreal cellVolume = getMethodData().getDistributionData().cell->computeVolume();
// 
//     const CFreal PI    = MathTools::MathConsts::CFrealPi();
//     //    const CFreal h     = std::sqrt(4*cellVolume/PI);
// 
//     const CFreal ALPHA = (VAvg + aAvg)*h;
//     /////////////////////////////////////////////
// 
// 
//     // ********************************************************************** //
//     // Dissipative part of LxF scheme, and its sum
//     // store it in m_sumPhiN
//     m_sumPhiN = 0.0;
//     for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
// 
//       const CFuint I = iState       % _nbStatesInCell;
//       const CFuint J = (iState + 1) % _nbStatesInCell ;
//       const CFuint K = (iState + 2) % _nbStatesInCell ;
// 
//       m_phiN[iState] = ALPHA*( *tStates[I] - *tStates[J] + *tStates[I] - *tStates[K] ) / static_cast<CFreal>(_nbStatesInCell);
// 
//       for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
//         m_sumPhiN[iEq] += std::abs(m_phiN[iState][iEq]);
//       }
//     }
// 
// 
//     /// Compute cell gradient, and ///
//     /// Compute cell averaged state from nodal (conservative) variables ///
// 
//     const CFreal coeffGrad = 1./(2.*vol);
//     _dUdx = 0.;
//     _dUdy = 0.;
//     _dUdz = 0.;
// 
//     for (CFuint kStates = 0; kStates < _nbStatesInCell; ++kStates){
//       const CFreal nK_x = normals[cellID]->getNodalNormComp(kStates,XX);
//       const CFreal nK_y = normals[cellID]->getNodalNormComp(kStates,YY);
// 
//       _dUdx += nK_x * (*tStates[kStates]);
//       _dUdy += nK_y * (*tStates[kStates]);
// 
//     }
//     _dUdx *= coeffGrad;
//     _dUdy *= coeffGrad;
// 
//     /// end - Compute cell gradient and average conservative state///
//   // ********************************************************************** //
//   
//   
//   // ********************************************************************** //
//   // * Copied from SUPG_ShockCapturingStrategy.cxx ************************ //
//     /// Compute gradients of xi, eta (and dseta) variables - Rao book on FEM, eq. 3.76 and eq. 3.75///
//   /// Eq. 3.75 : dfdx_j = sum_{k=1}^{3} { dLkdxj * dfdLk}
// 
//   /// Eq. 3.76 for a P1 triangle ///
// 
//   GeometricEntity *const m_cell = m_cellBuilder.buildGE();
// 
//   const Node& node0 = *m_cell->getNode(N0);
//   const Node& node1 = *m_cell->getNode(N1);
//   const Node& node2 = *m_cell->getNode(N2);
//   
//   const CFreal x0 = node0[XX]; const CFreal y0 = node0[YY];
//   const CFreal x1 = node1[XX]; const CFreal y1 = node1[YY];
//   const CFreal x2 = node2[XX]; const CFreal y2 = node2[YY];
// 
// //   std::cout << "volume: " << volume << "\n";
// //   std::cout << "x0 y0: " << x0 << " " << y0 << "\n";
// //   std::cout << "x1 y1: " << x1 << " " << y1 << "\n";
// //   std::cout << "x2 y2: " << x2 << " " << y2 << "\n";
// //   std::cout << 0.5*(x1*y2 +x2*y0 + x0*y1 - y0*x1 - y1*x2 - y2*x0) << "\n";
//   
//   
//   const CFreal dL0dx = (y1 - y2)*coeffGrad;
//   const CFreal dL0dy = (x2 - x1)*coeffGrad;
// 
//   const CFreal dL1dx = (y2 - y0)*coeffGrad;
//   const CFreal dL1dy = (x0 - x2)*coeffGrad;
// 
//   const CFreal dL2dx = (y0 - y1)*coeffGrad;
//   const CFreal dL2dy = (x1 - x0)*coeffGrad;
// 
//   _M(XX,0) = dL0dx;  _M(XX,1) = dL1dx;  _M(XX,2) = dL2dx;
//   _M(YY,0) = dL0dy;  _M(YY,1) = dL1dy;  _M(YY,2) = dL2dy;
//   /// end - Eq. 3.76 ///
// 
//   /// Vectors containing dXidLk, dEtadLk, dDsetadLk
// //   RealVector dXi(3) ;  dXi[0] = -1.; dXi[1]  = 1.;  dXi[2]  = -1.;
// //   RealVector dEta(3); dEta[0] = -1.; dEta[1] = -1.; dEta[2] = 1.;
//   RealVector dXi(3);  dXi[0]  = -1.; dXi[1]  = 1.; dXi[2]  = 0.;
//   RealVector dEta(3); dEta[0] = -1.; dEta[1] = 0.; dEta[2] = 1.;
//   /// end - Vectors containing dXidLk, dEtadLk, dDsetadLk
// 
//   _gradXi  = _M*dXi;
//   _gradEta = _M*dEta;
//   
//   /// end -  Compute gradients of xi, eta (and dseta) variables
// 
//   /// Matrix transforming Conservation to Entropy variables: _dVdU///
//   /// dV = _dVdU * dU
//   /// as defined in
//   /// CompMethAppMechEng, vol 89, pp 141-219, Shakib et al.
//   this->getMethodData().getLinearVar()->computePhysicalData( _Uavg, _pData);
// 
//   const CFreal rho = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::RHO];
//   const CFreal u   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::VX];
//   const CFreal v   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::VY];
//   const CFreal VSq = u*u + v*v;
//   const CFreal a   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::A];
//   const CFreal aSq = a*a;
// 
//   const CFreal E   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::E];
//   const CFreal H   = _pData[COOLFluiD::Physics::NavierStokes::EulerTerm::H];
// 
//   const CFreal gamma    = 1.4;
//   const CFreal gammaBar = gamma - 1.;
// 
//   const CFreal rhoU = rho*u;
//   const CFreal rhoV = rho*v;
//   const CFreal rhoE = rho*E;
// 
//   const CFreal rhoe = rhoE - 0.5*rho*VSq;
//   const CFreal s = std::log( (gamma-1.)*rhoe/std::pow(rho, gamma) );
// 
//   const CFreal rho0   = 1.4;
//   const CFreal rhoU0  = -12.6;
//   const CFreal rhoE0  = 59.2;
//   
//   const CFreal rhoe0 = rhoE0 - 0.5*rhoU0*rhoU0/rho0;
//   const CFreal s0 = std::log( (gamma-1.)*rhoe0/std::pow(rho0, gamma) );
// 
//   /// Entropy variable:
//   const CFreal V1 = ( -rhoE + rhoe*(gamma + 1. - s + s0) ) / rhoe;
//   const CFreal V2 = rhoU / rhoe;
//   const CFreal V3 = rhoV / rhoe;
//   const CFreal V5 = -rho / rhoe;
// 
//   // Auxiliary constants:
//   const CFreal k1 = 0.5*(V2*V2 + V3*V3)/V5;
//   const CFreal k2 = k1 - gamma;
//   const CFreal k3 = k1*k1 - 2.*gamma*k1 + gamma;
// 
//   const CFreal c1 = gammaBar*V5 - V2*V2;
//   const CFreal c2 = gammaBar*V5 - V3*V3;
// 
//   const CFreal d1 = -V2*V3;
// 
//   const CFreal e1 = V2*V5;
//   const CFreal e2 = V3*V5;
// 
//   _dVdU = 0.;
// 
//   _dVdU(rhoID, rhoID)  = k1*k1 + gamma;
//   _dVdU(rhoID, rhoUID) = k1*V2;
//   _dVdU(rhoID, rhoVID) = k1*V3;
//   _dVdU(rhoID, rhoEID) = (k1 + 1.)*V5;
// 
//   _dVdU(rhoUID, rhoID)  = k1*V2;
//   _dVdU(rhoUID, rhoUID) = V2*V2 - V5;
//   _dVdU(rhoUID, rhoVID) = -d1;
//   _dVdU(rhoUID, rhoEID) = e1;
// 
//   _dVdU(rhoVID, rhoID)  = k1*V3;
//   _dVdU(rhoVID, rhoUID) = -d1;
//   _dVdU(rhoVID, rhoVID) = V3*V3 - V5;
//   _dVdU(rhoVID, rhoEID) = e2;
// 
//   _dVdU(rhoEID, rhoID)  = (k1 + 1.)*V5;
//   _dVdU(rhoEID, rhoUID) = e1;
//   _dVdU(rhoEID, rhoVID) = e2;
//   _dVdU(rhoEID, rhoEID) = V5*V5;
// 
//   _dVdU *= -1./(rhoe*V5);
// 
//   /// Verifcation: _dVdU*_dUdV == Id
// //   _dUdV = 0.;
// // 
// //   _dUdV(rhoID, rhoID)  = -V5*V5;
// //   _dUdV(rhoID, rhoUID) = e1;
// //   _dUdV(rhoID, rhoVID) = e2;
// //   _dUdV(rhoID, rhoEID) = V5*(1. - k1);
// // 
// //   _dUdV(rhoUID, rhoID)  = e1;
// //   _dUdV(rhoUID, rhoUID) = c1;
// //   _dUdV(rhoUID, rhoVID) = d1;
// //   _dUdV(rhoUID, rhoEID) = k2*V2;
// // 
// //   _dUdV(rhoVID, rhoID)  = e2;
// //   _dUdV(rhoVID, rhoUID) = d1;
// //   _dUdV(rhoVID, rhoVID) = c2;
// //   _dUdV(rhoVID, rhoEID) = k2*V3;
// // 
// //   _dUdV(rhoEID, rhoID)  = V5*(1. - k1);
// //   _dUdV(rhoEID, rhoUID) = k2*V2;
// //   _dUdV(rhoEID, rhoVID) = k2*V3;
// //   _dUdV(rhoEID, rhoEID) = -k3;
// // 
// //   _dUdV *= rhoe/(gammaBar*V5);
// // 
// /*  RealMatrix tempProduct1 = _dVdU*_dUdV;
//   RealMatrix tempProduct2 = _dUdV*_dVdU;
// //   std::cout << "_dVdU*_dUdV\n " << (tempProduct1) << std::endl;
// //   std::cout << "_dUdV*_dVdU\n " << (tempProduct2) << std::endl;
//   for (CFuint ii = 0 ; ii < 4 ; ii++){
//     for (CFuint jj = 0 ; jj < 4 ; jj++){
//       if(ii != jj){
//         assert(std::abs(tempProduct1(ii,jj)) < 1.e-7 );
//         assert(std::abs(tempProduct2(ii,jj)) < 1.e-7 );
//       }
//       else{
//         assert( std::abs(tempProduct1(ii,ii) - 1.) < 1.e-7 );
//         assert( std::abs(tempProduct2(ii,ii) - 1.) < 1.e-7 );
// 
//       }
//     }
//   }*/
//   
//   /// end - Matrix _dVdU///
//   
//   /// Compute delta
// 
// //   const CFreal numerator = std::max(1.e-15, getNorm(cellResidual) );
// //   const CFreal numeratorSq = getNormSq(cellResidual);
// 
//   const RealVector& phiT = getMethodData().getDistributionData().phi;
//   
//   const CFreal numeratorSq = getNormSq( phiT/static_cast<CFreal>(vol) );/*BEWARE 1./vol! Is it correct? Yes, it is.*/
// 
//   /// Compute the denominator of delta like in Kirk - IJNMF 2010
// //   RealVector vXi  = _gradXi[XX] *_dUdx + _gradXi[YY] *_dUdy;
// //   RealVector vEta = _gradEta[XX]*_dUdx + _gradEta[YY]*_dUdy;
// // 
// //   const CFreal denXiSq  = getNormSq(vXi) ; 
// //   const CFreal denEtaSq = getNormSq(vEta); 
// //   
// //   const CFreal denominatorSq = denXiSq + denEtaSq;
// 
//   /// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   /// Compute the denominator of delta like in Kirk et al. - SUPG 2013
//   //_dUdx*invA0*dUdx + gXY*(_dUdx*invA0dUdy + _dUdy*temp_dVdU_DOT_dUdx) + gYY*_dUdy*temp_dVdU_DOT_dUdy
//   
//   // g^{xI,xJ} = dxIdcsiK * dyIdcsiK
//   // Use eq. 3.71, Rao book
//   // x = x0*L0 + x1*L1 + x2*L2; = x0*(1. - csi - eta ) + x1*csi + x2*eta = x0 + (x1-x0)*csi + (x2-x0)*eta
//   // y = y0*L0 + y1*L1 + y2*L2; = y0*(1. - csi - eta ) + y1*csi + y2*eta = y0 + (y1-y0)*csi + (y2-y0)*eta
//     
//   const CFreal gXX = (x1-x0)*(x1-x0) + (x2-x0)*(x2-x0);
//   const CFreal gXY = (x1-x0)*(y1-y0) + (x2-x0)*(y2-y0);
//   const CFreal gYY = (y1-y0)*(y1-y0) + (y2-y0)*(y2-y0);
// 
//   RealVector temp_dVdU_DOT_dUdx = _dVdU * _dUdx;
//   RealVector temp_dVdU_DOT_dUdy = _dVdU * _dUdy;
// 
//   CFreal pxx = 0.;
//   CFreal pxy = 0.;
//   CFreal pyx = 0.;
//   CFreal pyy = 0.;
//   
//   for (CFuint iEq = 0 ; iEq < _nbEquations; iEq++ ){
//     pxx += _dUdx[iEq]*temp_dVdU_DOT_dUdx[iEq];
//     pxy += _dUdx[iEq]*temp_dVdU_DOT_dUdy[iEq];
//     pyx += _dUdy[iEq]*temp_dVdU_DOT_dUdx[iEq];
//     pyy += _dUdy[iEq]*temp_dVdU_DOT_dUdy[iEq];
//   }
// 
//   const CFreal denominatorSq = ( gXX*pxx + gXY*(pxy + pyx) + gYY*pyy ); assert (0. < denominatorSq);
// 
//   //   std::cout << "denominator - denominatorSq2 = " << (denominatorSq - denominatorSq2) << std::endl;
// 
//   /// end -  Compute gradients of xi, eta (and dseta) variables
// /// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 
//   const CFreal delta = std::sqrt( numeratorSq/denominatorSq );
//   
//   // * end - from SUPG_ShockCapturingStrategy.cxx ************************ //
// 
//   // Modify _dUdxj, to have drhoHdxj - IJNMF 2010 Kirk
//   // Matrix dU^HdU is needed, where U^H = [rho rhoU rhoV rhoH ]'
//   // rhoH = rhoE + p = rhoE;
//   // p = Rg*rho*T = (Rg/cv)*rho*e = gamma_MINUS_one*(rhoE - 0.5*(rhoU**2 + rhoV**2)*rho**-1 );
//   // Finally rhoH = gamma*rhoE - 0.5*gamma_MINUS_one*(rhoU**2 + rhoV**2)*rho**-1;
// 
//   _dUHdU = 0.;
//   
//   _dUHdU(rhoID , rhoID ) = 1.;
// 
//   _dUHdU(rhoUID, rhoUID) = 1.;
// 
//   _dUHdU(rhoVID, rhoVID) = 1.;
// 
//   _dUHdU(rhoHID, rhoID ) = 0.5 * gammaBar * (u*u + v*v);
//   _dUHdU(rhoHID, rhoUID) = -gammaBar * u;
//   _dUHdU(rhoHID, rhoVID) = -gammaBar * v;
//   _dUHdU(rhoHID, rhoHID) = gamma;
//   
//  // end - IJNMF 2010 Kirk
//     
//     computeBlendingCoeff();
//     
//     if ( m_store_thetas ) storeThetas();
//     
//     if ( m_addExtraDiss ) addExtraDissipation(residual);
//     
//     
//     // computation of LDA residual and the blending
//     for (CFuint iState = 0; iState < _nbStatesInCell; ++iState)
//     {
//       _stab_Adv   = (1./vol) * ( (*_kPlus[iState]) + (*_kMin[iState]) ) * _Tau *  phiT;
// 
//       const CFreal nK_x = normals[cellID]->getNodalNormComp(iState, XX);
//       const CFreal nK_y = normals[cellID]->getNodalNormComp(iState, YY);
//       
//       _stab_Shock = ALPHA*(0.5*nK_x*_dUHdU*_dUdx + 0.5*nK_y*_dUHdU*_dUdy);  // WORKS 
// //       _dUHdU*m_phiN[iState]; // RECIRCULATION
//        // ALPHA * 0.5 * ( nK_x*(gXX*_dUHdU*_dUdx + gXY*_dUHdU*_dUdy) + nK_y*(gXY*_dUHdU*_dUdx + gYY*_dUHdU*_dUdy) );//BLOWS UP
// //       delta * 0.5 * ( nK_x*(gXX*_dUHdU*_dUdx + gXY*_dUHdU*_dUdy) + nK_y*(gXY*_dUHdU*_dUdx + gYY*_dUHdU*_dUdy) ); // DOESN'T WORK
//       // m_phiN[iState];
//       // ALPHA*(0.5*nK_x*_dUdx + 0.5*nK_y*_dUdy);
//       
//   	  for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
//         residual[iState][iEq] = one_OVER_dimPlusOne*phiT[iEq] + _stab_Adv[iEq] + m_theta[iEq]*_stab_Shock[iEq];
//       }
// 
//       if (ddata.computeBetas)
//       {
//         (*ddata.currBetaMat)[iState] = (1./vol) * ( (*_kPlus[iState]) + (*_kMin[iState]) ) * _Tau; // vdW, still missing 1./3. in the diagonal
//         for(CFuint iEq = 0; iEq < _nbEquations; ++iEq)
//           (*ddata.currBetaMat)[iState](iEq, iEq) += one_OVER_dimPlusOne;
//       }
//     }
//   }
//   else {
// 	std::cout << "Oh, why!\n";abort();
//     assert(m_firstOrder == 1 || (m_firstOrderJacob && isPerturb));
//     BSchemeBase<NSchemeCSys>::distribute(residual);
//   }
// }
//       
// //////////////////////////////////////////////////////////////////////////////
