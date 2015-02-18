#include "NSchemeCSysNEQ.hh"
#include "MathTools/MatrixInverter.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitSystem.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/MeshData.hh"

#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::FluctSplit;

typedef MultiScalarTerm<COOLFluiD::Physics::NavierStokes::EulerTerm> NEQTerm;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplitNEQ {

CFuint NSchemeCSysNEQ::_actual_iter = 0;

CFuint NSchemeCSysNEQ::_last_accessed_at_iter = 0;

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<NSchemeCSysNEQ,
                       FluctSplit::FluctuationSplitData,
                       FluctSplit::Splitter,
                       FluctSplit::FluctSplitSystemModule>
ncNEQSchemeSysProvider("SysNCNEQ");

//////////////////////////////////////////////////////////////////////////////

NSchemeCSysNEQ::NSchemeCSysNEQ(const std::string& name): NSchemeCSys(name),
socket_ExtraDiss_active("extra_dissipation_active"),
_doAct(0.),
_varID(0)
{
  this->NSchemeCSys::addConfigOptionsTo(this);

  _deltaVar = 0.0;
  this->setParameter("Delta",&_deltaVar);

  _length = 1.0;
  this->setParameter("Length",&_length);

  _speed = 1.0;
  this->setParameter("Speed",&_speed);

  _varName = "p";
  this->setParameter("VarName",&_varName);
  
 _freezeAdditionalDiss = 0;
 this->setParameter("FreezeAdditionalDissipation", &_freezeAdditionalDiss);
}

//////////////////////////////////////////////////////////////////////////////

NSchemeCSysNEQ::~NSchemeCSysNEQ()
{
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCSysNEQ::configure ( Config::ConfigArgs& args )
{
  NSchemeCSys::configure(args);  
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCSysNEQ::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal>("Delta","Delta of variable.");
  options.addConfigOption< CFreal>("Length","Reference Length.");
  options.addConfigOption< CFreal>("Speed","Reference Speed.");
  options.addConfigOption< std::string>("VarName","Variable name.");
  options.addConfigOption< CFuint, Config::DynamicOption<> >
    ("FreezeAdditionalDissipation", "The switch for the additional dissipation is not recomputed anymore");
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCSysNEQ::setup()
{
  NSchemeCSys::setup();

  _Uavg.resize(_nbEquations);

  _model = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm().d_castTo<NEQTerm>();
  _model->resizePhysicalData(_pData);

  _gradVar.resize(2);
  
  // set _choiceVar to pressure if default value is required
  if (_varName == "rho")
  {
    _varID = Physics::NavierStokes::EulerTerm::RHO;
  }

  if (_varName == "p")
  {
    _varID = Physics::NavierStokes::EulerTerm::P;
  }

  if (_varName == "T")
  {
    _varID = Physics::NavierStokes::EulerTerm::T;
  }

  if ( _varName != "p" && _varName != "rho" && _varName != "T" )
    throw Common::BadValueException (FromHere(),"Variable to base the shock detector must be either [rho], [T] or [p]"); 
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCSysNEQ::distribute(vector<RealVector>& residual)
{
  const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
  const RealVector& phiT = getMethodData().getDistributionData().phi;


  this->_tmp = (*_kPlus[0]) + (*_kMin[0]);

  this->_sumKU = this->_tmp*(*tStates[0]);
  this->_sumKplusU = (*_kPlus[0]) * (*tStates[0]);
  this->_sumKplus = *_kPlus[0];
  
  for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
    this->_sumKplusU += (*_kPlus[iState])*(*tStates[iState]);
    this->_sumKplus += *_kPlus[iState];

    this->_tmp = (*_kPlus[iState])  + (*_kMin[iState]);
    this->_sumKU += this->_tmp*(*tStates[iState]);
  }

  this->_inverter->invert(this->_sumKplus, this->_invK);


  CFLogDebugMax( "invK = " << "\n" << (this->_invK) << "\n");

  this->_uInflow = this->_invK * (this->_sumKplusU - this->_sumKU);

  CFLogDebugMax( "uInflow = " << "\n" << (this->_uInflow) << "\n");

  for (CFuint iState = 0; iState < this->_nbStatesInCell; ++iState) {
    residual[iState] = (*_kPlus[iState])*(*tStates[iState] - this->_uInflow);

    RealMatrix& betaLDA = (*getMethodData().getDistributionData().currBetaMat)[iState];
    betaLDA = (*_kPlus[iState])*this->_invK;
 
    residual[iState] -= betaLDA*(this->_sumKU - phiT);
  }

  _doAct = 0.;
  _actual_iter = SubSystemStatusStack::getActive()->getNbIter();
  
  const CFuint N0 = 0;
  const CFuint N1 = 1;
  const CFuint N2 = 2;
  
  // RealVector alpha(_nbEquations);
  
  // CFreal alpha_max = 0.;
  
  //       for (CFuint iEqs = 0; iEqs < _nbEquations; iEqs++){
  // 
  //         const CFreal refVal = std::max(std::abs(phiT[iEqs]), std::abs(this->_sumKU[iEqs]));//??
  // //         const CFreal refVal = std::abs(phiT[iEqs]);
  // 
  //         alpha[iEqs] = refVal > 1.e-12? std::abs(this->_sumKU[iEqs] - phiT[iEqs])/refVal: 0.;
  // 
  //         assert( !(alpha[iEqs] != alpha[iEqs]) );
  //          
  //         if( alpha[iEqs] > alpha_max) alpha_max = alpha[iEqs];
  //        
  //       }//std::cout << "\nalpha_max = " << alpha_max <<  std::endl;
  // 
  //       assert(alpha_max >= 0.);assert(alpha_max <= 1.);
  
  /*-----*/
  
  _Uavg = (*tStates[N0] + *tStates[N1] + *tStates[N2])/this->_nbStatesInCell;
  getMethodData().getLinearVar()->computePhysicalData(_Uavg, _pData);
  
  const CFuint cellID = getMethodData().getDistributionData().cell->getID();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  FluctSplit::DistributionData& ddata = getMethodData().getDistributionData();
  const CFreal cellVolume = ddata.cell->computeVolume();
  const CFreal char_size = std::sqrt(cellVolume);
  
  const CFreal a = _pData[EulerTerm::A];
  const CFreal V = _pData[EulerTerm::V];
  
  /* ------------------------------------------------------------------------- */
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFuint uID = nbSpecies;
  const CFuint vID = nbSpecies + 1;
  
  CFreal aux_density = 0.;
  
  for(CFuint iSpecies = 0 ; iSpecies < nbSpecies; iSpecies++){
    aux_density += _Uavg[iSpecies];
  }
  
  const CFreal rho = aux_density;
  
  
  DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();
  
  DataHandle<CFreal> volumes             = socket_volumes.getDataHandle();
  const CFreal vol = volumes[cellID];
  
  _gradVar[XX] = 0.;
  _gradVar[YY] = 0.;
  
  for(CFuint iState = 0; iState < _nbStatesInCell; iState++){
    getMethodData().getLinearVar()->computePhysicalData(*tStates[iState], _pData);
    
    const CFreal nodalVar = _pData[_varID];
    
    for(CFuint iDim = 0; iDim < dim; iDim++){
      const CFreal n_s__d = normals[cellID]->getNodalNormComp(iState,iDim);
      _gradVar[iDim] += nodalVar*n_s__d;
    }
  }
  _gradVar /= (dim*vol);
  
  const CFreal ux = _Uavg[uID]/rho;
  const CFreal uy = _Uavg[vID]/rho;
  const CFreal l_OVER_UdeltaVar_ref = _length/(_speed*_deltaVar);
  const CFreal shockDetector_aux = (_gradVar[XX]*ux + _gradVar[YY]*uy)*l_OVER_UdeltaVar_ref;
  const CFreal shockDetector = std::max(0.,shockDetector_aux);
  
  // Diameter of a circle with the same area of the element:
  const CFreal h = 2.0*std::sqrt(static_cast<CFreal>(vol/MathTools::MathConsts::CFrealPi()));
  
  //       const CFreal theta = std::min(1., h*shockDetector*shockDetector);
  //       const CFreal theta = std::min(1., shockDetector*shockDetector);
  CFreal theta = std::min(1., (h/_length)*shockDetector*shockDetector);
  
  if( _freezeAdditionalDiss == 0 ) {
    if (_last_accessed_at_iter != SubSystemStatusStack::getActive()->getNbIter()){
      socket_ExtraDiss_active.getDataHandle() = 0;
    }
    
    socket_ExtraDiss_active.getDataHandle()[cellID] = theta;
  }
  else {
    theta = socket_ExtraDiss_active.getDataHandle()[cellID];
  }
  
  //       store_ExtraDiss(theta);
  
  //         if (true) {//if (_store_cells_extraDiss == 1) {
  //                 store_ExtraDiss(theta);
  //               }
  //       Needed to detect when the next timestep has just  been achieved:
  _last_accessed_at_iter = SubSystemStatusStack::getActive()->getNbIter();
  
  /*------*/
  //        CFreal K_diss = 0.; // Doesn't work: T < Tinf, and then blows up!!!
  //       CFreal K_diss = 0.01*V*char_size/3.;// Doesn't work!
  //       CFreal K_diss = 0.05*V*char_size/3.;
  //         CFreal K_diss = 0.1*V*char_size/3.; // Works!!! But they need the Carbuncle Fix ...
  const CFreal K_diss = 0.2*V*char_size/3.; // Works!!! But they need the Carbuncle Fix ...
  //		K_diss = V*char_size;
  //         CFreal K_diss = V*char_size/3.; // Works!!! Very diffusive, not improving much w.r.t. convergence for carbuncle fix active
  //         CFreal K_diss = V*char_size; // Doesn't work!
  
  //       CFreal K_diss = a*char_size/3.;      // Doesn't work
  //       CFreal K_diss = 0.1*a*char_size/3.;  // Doesn't work
  //       CFreal K_diss = 0.05*a*char_size/3.; // Doesn't work
  
  const CFreal u0 = std::sqrt((*tStates[N0])[uID]*(*tStates[N0])[uID] + (*tStates[N0])[vID]*(*tStates[N0])[vID]);
  const CFreal u1 = std::sqrt((*tStates[N1])[uID]*(*tStates[N1])[uID] + (*tStates[N1])[vID]*(*tStates[N1])[vID]);
  const CFreal u2 = std::sqrt((*tStates[N2])[uID]*(*tStates[N2])[uID] + (*tStates[N2])[vID]*(*tStates[N2])[vID]);
  
  const CFreal uMin = std::min(std::min(u0,u1),u2);
  const CFreal uMax = std::max(std::max(u0,u1),u2);
  
  //         CFreal K_diss = (uMax - uMin)*h/3.; //Doesn't work
  //         CFreal K_diss = (uMax - uMin)*h;//Doesn't work
  //         CFreal K_diss = 0.5*(uMax + uMin)*h;//Doesn't work
  //         CFreal K_diss = 0.1*0.5*(uMax + uMin)*h;//Doesn't work
  //         CFreal K_diss = (u0 + u1 + u2)*h/3.;
  //       CFreal K_diss = (u0 + u1 + u2)*h;//Blows up!
  //       CFreal K_diss = (u0 + u1 + u2)*char_size/3.;
  /*------*/
  
  residual[0] += theta*K_diss*(2.*(*tStates[N0]) - *tStates[N1] - *tStates[N2] );
  residual[1] += theta*K_diss*(2.*(*tStates[N1]) - *tStates[N2] - *tStates[N0] );
  residual[2] += theta*K_diss*(2.*(*tStates[N2]) - *tStates[N0] - *tStates[N1] );
  //       
  //         const CFreal rhoE0 = (*tStates[N0])[4];
  //         const CFreal rhoE1 = (*tStates[N1])[4];
  //         const CFreal rhoE2 = (*tStates[N2])[4];
  // 
  //         const CFreal rhoEMin = std::min(std::min(rhoE0,rhoE1),rhoE2);
  //         const CFreal rhoEMax = std::max(std::max(rhoE0,rhoE1),rhoE2);
  // 
  //       residual[0][4] += theta*V*char_size*(rhoE0 - 0.5*(rhoE1 + rhoE2))/3.;
  //       residual[1][4] += theta*V*char_size*(rhoE1 - 0.5*(rhoE2 + rhoE0))/3.;
  //       residual[2][4] += theta*V*char_size*(rhoE2 - 0.5*(rhoE0 + rhoE1))/3.;
  
  /* ------------------------------------------------------------------------- */
  
  /* ------------------------------------------------------------------------- */
  
  //       _doAct = alpha_max > 1.e-11? 1.: 0.;
  // 
  //       if( _freezeAdditionalDiss == 0 ) {
  //         if (_last_accessed_at_iter != SubSystemStatusStack::getActive()->getNbIter()){
  //           socket_ExtraDiss_active.getDataHandle() = 0;
  //         }
  // 
  //         socket_ExtraDiss_active.getDataHandle()[cellID] = _doAct;
  //       }
  //       else {
  //         _doAct = socket_ExtraDiss_active.getDataHandle()[cellID];
  //       }
  // //       Needed to detect when the next timestep has just  been achieved:
  //       _last_accessed_at_iter = SubSystemStatusStack::getActive()->getNbIter();
  // 
  //       CFreal K_diss = 0.1*V*char_size/3.; // Works!!! But they need the Carbuncle Fix ...
  // 
  //       residual[0] += _doAct*K_diss*(2.*(*tStates[N0]) - *tStates[N1] - *tStates[N2] );
  //       residual[1] += _doAct*K_diss*(2.*(*tStates[N1]) - *tStates[N2] - *tStates[N0] );
  //       residual[2] += _doAct*K_diss*(2.*(*tStates[N2]) - *tStates[N0] - *tStates[N1] );
  //       
  /* ------------------------------------------------------------------------- */
}

//////////////////////////////////////////////////////////////////////////////

void NSchemeCSysNEQ::distributePart(vector<RealVector>& residual)
{
  throw Common::NotImplementedException (FromHere(),"NSchemeCSysNEQ::distributePart()");
}
  
//////////////////////////////////////////////////////////////////////////////

void NSchemeCSysNEQ::computePicardJacob(vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"NSchemeCSysNEQ::computePicardJacob");
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

void NSchemeCSysNEQ::store_ExtraDiss(const CFreal value)
{
  using namespace COOLFluiD::Framework;

//  cf_assert(_store_cells_extraDiss);
  FluctSplit::DistributionData& distdata = getMethodData().getDistributionData();

  if (distdata.isPerturb) return; // skip if is being perturbed
  DataHandle< CFreal > extra_dissipation_active = socket_ExtraDiss_active.getDataHandle();

  const CFuint cellID = getMethodData().getDistributionData().cell->getID();
//   extra_dissipation_active[cellID] = _doAct? alpha_max : 0.;
  extra_dissipation_active[cellID] = value > 0.? value : 0.;

  //Reset the socket for the next iteration:
  if( _last_accessed_at_iter != _actual_iter ){
    extra_dissipation_active = 0;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////


// 2012_06_06
// void NSchemeCSysNEQ::distribute(vector<RealVector>& residual)
// {
// //   std::cout << "void NSchemeCSysNEQ::distribute(vector<RealVector>& residual)\n";
//   // AL: OLD implementation left here for the moment (performance comparison needed)
//   //  _sumKplusU = (*_kPlus[0]) * (*tStates[0]);
//   //   _sumKplus = *_kPlus[0];
//   //   for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
//   //     _sumKplusU += (*_kPlus[iState])*(*tStates[iState]);
//   //     _sumKplus += *_kPlus[iState];
//   //   }
// 
//   //   CFLogDebugMax( "sumKplusU = " << _sumKplusU << "\n");
//   //   CFLogDebugMax( "sumKplus = " << _sumKplus << "\n");
// 
//   //   _inverter->invert(_sumKplus, _invK);
// 
//   //   CFLogDebugMax( "invK = " << "\n" <<_invK << "\n");
// 
//   //   _uInflow = _invK * (_sumKplusU - phiT);
// 
//   //   CFLogDebugMax( "uInflow = " << "\n" << _uInflow << "\n");
// 
//   //   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
//   //     residual[iState] =
//   //       (*_kPlus[iState])*(*tStates[iState] - _uInflow);
//   //   }
// 
//   // beta LW
//   //   RealMatrix sumAbsK(_nbEquations,_nbEquations);
//   //   RealMatrix invSumAbsK(_nbEquations,_nbEquations);
//   //   RealMatrix k(_nbEquations,_nbEquations);
//   //   const CFuint dimPlus1 = PhysicalModelStack::getActive()->getDim() + 1;
// 
//   const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
//   const RealVector& phiT = getMethodData().getDistributionData().phi;
// 
//   _tmp = (*_kPlus[0])  + (*_kMin[0]);
// 
//   //  for (CFuint i=0; i < tmp.size(); ++i) {
//   //     sumAbsK[i] = std::abs(tmp[i]);
//   //   }
// 
//   _sumKU = _tmp*(*tStates[0]);
//   _sumKplusU = (*_kPlus[0]) * (*tStates[0]);
//   _sumKplus = *_kPlus[0];
//   for (CFuint iState = 1; iState < _nbStatesInCell; ++iState) {
//     _sumKplusU += (*_kPlus[iState])*(*tStates[iState]);
//     _sumKplus += *_kPlus[iState];
// 
//     _tmp = (*_kPlus[iState])  + (*_kMin[iState]);
//     _sumKU += _tmp*(*tStates[iState]);
// 
//     //     for (CFuint i=0; i < tmp.size(); ++i) {
//     //       sumAbsK[i] += std::abs(tmp[i]);
//     //     }
//   }
// 
//   //  for (CFuint iState = 0; iState < _kPlus.size(); ++iState) {
//   // for (CFuint i = 0; i < _nbEquations; ++i) {
//   //   (*_kPlus[iState])(i,i) += 1e-8;
//   //  }
//   // }
// 
//   _inverter->invert(_sumKplus, _invK);
// 
//   // _inverter->invert(sumAbsK, invSumAbsK);
// 
//   CFLogDebugMax( "invK = " << "\n" <<_invK << "\n");
// 
//   _uInflow = _invK * (_sumKplusU - _sumKU);
// 
//   CFLogDebugMax( "uInflow = " << "\n" << _uInflow << "\n");
// 
//   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
//     residual[iState] = (*_kPlus[iState])*(*tStates[iState] - _uInflow);
// 
//     RealMatrix& betaLDA = (*getMethodData().getDistributionData().currBetaMat)[iState];
//     betaLDA = (*_kPlus[iState])*_invK;
// 
//     // L_W
//     // k = (*_kPlus[iState])  + (*_kMin[iState]);
//     //     tmp = k*invSumAbsK;
//     //     tmp *= CFL::getInstance().getCFL();
//     //     for (CFuint i=0; i < _nbEquations; ++i) {
//     //       tmp(i,i) += (1./dimPlus1);
//     //     }
// 
//     //    RealVector verr(_nbEquations);
//     //verr = sumKU - phiT;
// 
//     //  cout << "sumKU = " << sumKU << endl;
//     //     cout << "phiT  = " << phiT << endl;
//     //     RealVector t(_nbEquations);
//     //     t = tmp*(sumKU - phiT);
//     //     cout << "t  = " << t << endl << endl;
// 
//     //    if (_isPerturb) {
//     //  _temp = _tempBkp;
//     //}
//     //else {
//     //  _temp = _tmp*(_sumKU - phiT);
//     //  _tempBkp = _temp;
//     // }
//     residual[iState] -= betaLDA*(_sumKU - phiT);
//   }
// 
//       /* ----JGM-------*/
// 
//     _enhance_dissipation = true;
//     if(_enhance_dissipation){
//       //       std::cout << "Enhancing dissipation\n";
// 
//       RealVector alpha(_nbEquations);
// 
//       CFreal alpha_max = 0.;
// 
//       for (CFuint iEqs = 0; iEqs < _nbEquations; iEqs++){
// 
//         const CFreal refVal = std::max(std::abs(phiT[iEqs]), std::abs(_sumKU[iEqs]));
// 
//         alpha[iEqs] = refVal > 1.e-12? std::abs(_sumKU[iEqs] - phiT[iEqs])/refVal: 0.;
// 
//         assert( !(alpha[iEqs] != alpha[iEqs]) );
// 
// //         if ( !(alpha[iEqs] != alpha[iEqs])){// alpha[iEqs] is NOT NaN
// 
//           if( alpha[iEqs] > alpha_max) alpha_max = alpha[iEqs];
// 
// //         }
// //         const CFreal ref = max(std::abs(phiT[iEqs]), 1.e-12);
// //
// //         alpha[iEqs] = std::abs(_sumKU[iEqs] - phiT[iEqs])/ref;
// 
//       }//std::cout << "\nalpha_max = " << alpha_max <<  std::endl;
// 
//       assert(alpha_max >= 0.);assert(alpha_max <= 1.);
// 
//       /*-----*/
// 
//       RealVector _Uavg = (*tStates[0] + *tStates[1] + *tStates[2])/_nbStatesInCell;
// 
//       RealVector _pData;
// 
//       typedef Framework::MultiScalarTerm<COOLFluiD::Physics::NavierStokes::EulerTerm> NEQTerm;
//       Common::SafePtr<NEQTerm> _model;
//       _model = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm().d_castTo<NEQTerm>();
// 
//       _model->resizePhysicalData(_pData);
// 
//       getMethodData().getLinearVar()->computePhysicalData(_Uavg, _pData);
// 
//       DistributionData& ddata = getMethodData().getDistributionData();
//       const CFreal cellVolume = ddata.cell->computeVolume();
//       const CFreal char_size = std::sqrt(cellVolume);
// 
//       //       std::cout << "_pData = "<< _pData << std::endl;
//       // a = _pData[4]; V = _pData[6];
//       const CFreal a = _pData[4];
//       const CFreal V = _pData[6];
// 
// //       CFreal K_diss = 0.1*(V + a)*char_size/3.;
// 
// //       CFreal K_diss = 0.1*(V + a)*char_size/3.;
//       // 0.001 is too small
// //       CFreal K_diss = 0.01*a*char_size/3.;
// 
// //       CFreal K_diss = 0.1*a*char_size/3.;
// 
// //       CFreal K_diss = 0.1*(a + V)*char_size/3.;
// 
// //       CFreal K_diss = std::abs(a - V)*char_size/3.;
// 
// //       CFreal K_diss = alpha_max*V*char_size/3.; Blows up!
// 
//       /*------*/
// //       CFreal K_diss = 0.;
//       CFreal K_diss = 0.1*V*char_size/3.; // Works!!! But they need the Carbuncle Fix ...
// //       CFreal K_diss = 0.01*V*char_size/3.;// Doesn't work!
// // //       CFreal K_diss = 0.001*V*char_size/3.;// Doesn't work!
// 
//       /*------*/
// 
//       CFreal doAct = alpha_max > 1.e-12? 1.: 0.;
// 
//       /*-----*/
// 
//       /*Nop*/
// //       residual[0] += alpha_max*K_diss*(2.*(*tStates[0]) - *tStates[1] - *tStates[2] );
// //       residual[1] += alpha_max*K_diss*(2.*(*tStates[1]) - *tStates[2] - *tStates[0] );
// //       residual[2] += alpha_max*K_diss*(2.*(*tStates[2]) - *tStates[0] - *tStates[1] );
// 
//       /*Nop*/
// //       residual[0] += alpha_max*_sumKplus*(2.*(*tStates[0]) - *tStates[1] - *tStates[2] )/3.;
// //       residual[1] += alpha_max*_sumKplus*(2.*(*tStates[1]) - *tStates[2] - *tStates[0] )/3.;
// //       residual[2] += alpha_max*_sumKplus*(2.*(*tStates[2]) - *tStates[0] - *tStates[1] )/3.;
// 
//     /*Psche...*/
// //       CFuint doAct = alpha_max > 1.e-12? 1.: 0.;
// //
// //       residual[0] += doAct*K_diss*(2.*(*tStates[0]) - *tStates[1] - *tStates[2] );
// //       residual[1] += doAct*K_diss*(2.*(*tStates[1]) - *tStates[2] - *tStates[0] );
// //       residual[2] += doAct*K_diss*(2.*(*tStates[2]) - *tStates[0] - *tStates[1] );
// 
// 
//     /* ---------------------------------------------------------------------*/
// //         RealMatrix D(_nbEquations,_nbEquations);
// //         RealMatrix S(_nbEquations,_nbEquations);
// //
// //         for (CFuint iEqs = 0; iEqs < _nbEquations; iEqs++){
// //           D(iEqs,iEqs) = alpha[iEqs];
// //         }
// 
//         //       // Upstream reference values
//         // //       CFreal rhoInf  = 0.0005952;
//         // //       CFreal rhoUInf = rhoInf*5590.;
//         // //       CFreal rhoEInf = rhoInf*rhoUInf;
// 
//         // Local reference values
// //         CFreal rhoInf  = _pData[0];
// //         CFreal rhoUInf = rhoInf*V;//~rho*A
// //         CFreal rhoEInf = rhoInf*rhoUInf;
// 
//         // //       D(0,0) = alpha[0]/rhoInf;
//         // //       D(1,1) = alpha[1]/rhoInf;
//         // //
//         // //       D(2,2) = alpha[2]/rhoUInf;
//         // //       D(3,3) = alpha[3]/rhoUInf;
//         // //
//         // //       D(4,4) = alpha[4]/rhoEInf;
//         // //       D(5,5) = alpha[5]/rhoEInf;
// 
//         // S is just a scaling matrix
// //         S(0,0) = 1./rhoInf;
// //         S(1,1) = 1./rhoInf;
// //
// //         S(2,2) = 1./rhoUInf;
// //         S(3,3) = 1./rhoUInf;
// //
// //         S(4,4) = 1./rhoEInf;
// //         S(5,5) = 1./rhoEInf;
// 
// //               D(0,0) = alpha[0]/rhoInf;
// //               D(1,1) = alpha[1]/rhoInf;
// //
// //               D(2,2) = alpha[2]/rhoUInf;
// //               D(3,3) = alpha[3]/rhoUInf;
// //
// //               D(4,4) = alpha[4]/rhoEInf;
// //               D(5,5) = alpha[5]/rhoEInf;
//         //
//         //
//         // //       D(0,0) = alpha[0];
//         // //       D(1,1) = alpha[1];
//         // //
//         // //       D(2,2) = alpha[2];
//         // //       D(3,3) = alpha[3];
//         // //
//         // //       D(4,4) = alpha[4];
//         // //       D(5,5) = alpha[5];
//         /* --- */
// //               residual[0] += doAct*D*(2.*(*tStates[0]) - *tStates[1] - *tStates[2] );
// //               residual[1] += doAct*D*(2.*(*tStates[1]) - *tStates[2] - *tStates[0] );
// //               residual[2] += doAct*D*(2.*(*tStates[2]) - *tStates[0] - *tStates[1] );
//         /* --- */
// 
//       /* --- */
// //       residual[0] += doAct*alpha_max*S*(2.*(*tStates[0]) - *tStates[1] - *tStates[2] );// D*
// //       residual[1] += doAct*alpha_max*S*(2.*(*tStates[1]) - *tStates[2] - *tStates[0] );// D*
// //       residual[2] += doAct*alpha_max*S*(2.*(*tStates[2]) - *tStates[0] - *tStates[1] );// D*
//       /* --- */
//       /* ---------------------------------------------------------------------*/
// 
// 
// //       residual[0] += doAct*D*(2.*(*tStates[0]) - *tStates[1] - *tStates[2] );
// //       residual[1] += doAct*D*(2.*(*tStates[1]) - *tStates[2] - *tStates[0] );
// //       residual[2] += doAct*D*(2.*(*tStates[2]) - *tStates[0] - *tStates[1] );
// 
//       //NOP!
// //       residual[0] += doAct*D*(2.*(*tStates[0]) - *tStates[1] - *tStates[2] );//K_diss*D*
// //       residual[1] += doAct*D*(2.*(*tStates[1]) - *tStates[2] - *tStates[0] );//K_diss*D*
// //       residual[2] += doAct*D*(2.*(*tStates[2]) - *tStates[0] - *tStates[1] );//K_diss*D*
// 
//       residual[0] += doAct*K_diss*(2.*(*tStates[0]) - *tStates[1] - *tStates[2] );
//       residual[1] += doAct*K_diss*(2.*(*tStates[1]) - *tStates[2] - *tStates[0] );
//       residual[2] += doAct*K_diss*(2.*(*tStates[2]) - *tStates[0] - *tStates[1] );
//     }
//     /* ----JGM-------*/
// }
// 
// //////////////////////////////////////////////////////////////////////////////
