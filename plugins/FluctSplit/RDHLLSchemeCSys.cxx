#include "RDHLLSchemeCSys.hh"
#include "MathTools/MatrixInverter.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FluctSplit/FluctSplitSystem.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/MeshData.hh"

#include "NavierStokes/EulerTerm.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RDHLLSchemeCSys,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitSystemModule>
hllcSchemeSysProvider("SysHLLC");

//////////////////////////////////////////////////////////////////////////////

RDHLLSchemeCSys::RDHLLSchemeCSys(const std::string& name) :
  RDS_SplitterSys(name),
  _sumKplusU(),
  _sumKplus(),
  _invK(),
  _uInflow(),
  _uDiff(),
  _temp(),
  _tempBkp(),
  _tempMat(),
  _tmp(),
  _sumKU(),
  n_IJ(),
  n_IK(),  
  F_I(),
  G_I(),
  F_J(),
  G_J(),
  F_K(),
  G_K(),  
  Fn_IJ(),
  Fn_IK()
{
}

//////////////////////////////////////////////////////////////////////////////

RDHLLSchemeCSys::~RDHLLSchemeCSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void RDHLLSchemeCSys::setup()
{
  RDS_SplitterSys::setup();
  
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
  
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->resizePhysicalData(_pData);  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  n_IJ.resize(dim);
  n_IK.resize(dim);
  
  F_I.resize(_nbEquations);  G_I.resize(_nbEquations);
  
  F_J.resize(_nbEquations);  G_J.resize(_nbEquations);
  
  F_K.resize(_nbEquations);  G_K.resize(_nbEquations);
  
  Fn_IJ.resize(_nbEquations);
  Fn_IK.resize(_nbEquations);  
}

//////////////////////////////////////////////////////////////////////////////

void RDHLLSchemeCSys::distribute(vector<RealVector>& residual)
{
//   std::cout << "void RDHLLSchemeCSys::distribute(vector<RealVector>& residual)\n";
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

  const vector<State*>& tStates = *(this->getMethodData().getDistributionData().tStates);
  
  const CFuint nbStates = tStates.size();  
    

  for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
	
	const CFuint I = iState       % nbStates;
	const CFuint J = (iState + 1) % nbStates ;
	const CFuint K = (iState + 2) % nbStates ;
	
	n_IJ[XX] = (1./6.) *( normals[cellID]->getNodalNormComp(J,XX) - normals[cellID]->getNodalNormComp(I,XX) );
	n_IJ[YY] = (1./6.) *( normals[cellID]->getNodalNormComp(J,YY) - normals[cellID]->getNodalNormComp(I,YY) );

	n_IK[XX] = (1./6.) *( normals[cellID]->getNodalNormComp(K,XX) - normals[cellID]->getNodalNormComp(I,XX) );
	n_IK[YY] = (1./6.) *( normals[cellID]->getNodalNormComp(K,YY) - normals[cellID]->getNodalNormComp(I,YY) );
	
	this->getMethodData().getLinearVar()->computePhysicalData(*tStates[I], _pData);
	const CFreal rho_I = _pData[EulerTerm::RHO];
	const CFreal u_I   = _pData[EulerTerm::VX];
	const CFreal v_I   = _pData[EulerTerm::VY];
	const CFreal p_I   = _pData[EulerTerm::P];
	const CFreal E_I   = _pData[EulerTerm::E];
	const CFreal V_I   = _pData[EulerTerm::V];
	const CFreal A_I   = _pData[EulerTerm::A];

	F_I = 0.;                             G_I = 0.; 
	F_I[rhoID]  = rho_I*u_I;              G_I[rhoID]  = rho_I*v_I;
	F_I[rhoUID] = rho_I*u_I*u_I + p_I;    G_I[rhoUID] = rho_I*u_I*v_I;
	F_I[rhoVID] = rho_I*v_I*u_I;          G_I[rhoVID] = rho_I*v_I*v_I + p_I;
	F_I[rhoEID] = (rho_I*E_I + p_I)*u_I;  G_I[rhoEID] = (rho_I*E_I + p_I)*v_I;

	this->getMethodData().getLinearVar()->computePhysicalData(*tStates[J], _pData);
	const CFreal rho_J = _pData[EulerTerm::RHO];
	const CFreal u_J   = _pData[EulerTerm::VX];
	const CFreal v_J   = _pData[EulerTerm::VY];
	const CFreal p_J   = _pData[EulerTerm::P];
	const CFreal E_J   = _pData[EulerTerm::E];
	const CFreal V_J   = _pData[EulerTerm::V];
	const CFreal A_J   = _pData[EulerTerm::A];

	F_J = 0.;                             G_J = 0.;
	F_J[rhoID]  = rho_J*u_J;              G_J[rhoID]  = rho_J*v_J;
	F_J[rhoUID] = rho_J*u_J*u_J + p_J;    G_J[rhoUID] = rho_J*u_J*v_J ;
	F_J[rhoVID] = rho_J*v_J*u_J;          G_J[rhoVID] = rho_J*v_J*v_J + p_J;
	F_J[rhoEID] = (rho_J*E_J + p_J)*u_J;  G_J[rhoEID] = (rho_J*E_J + p_J)*v_J;

	this->getMethodData().getLinearVar()->computePhysicalData(*tStates[K], _pData);
	const CFreal rho_K = _pData[EulerTerm::RHO];
	const CFreal u_K   = _pData[EulerTerm::VX];
	const CFreal v_K   = _pData[EulerTerm::VY];
	const CFreal p_K   = _pData[EulerTerm::P];
	const CFreal E_K   = _pData[EulerTerm::E];
	const CFreal V_K   = _pData[EulerTerm::V];
	const CFreal A_K   = _pData[EulerTerm::A];

	F_K = 0.;                             G_K = 0.;
	F_K[rhoID]  = rho_K*u_K;              G_K[rhoID]  = rho_K*v_K;
	F_K[rhoUID] = rho_K*u_K*u_K + p_K;    G_K[rhoUID] = rho_K*u_K*v_K;
	F_K[rhoVID] = rho_K*v_K*u_K;          G_K[rhoVID] = rho_K*v_K*v_K + p_K;
	F_K[rhoEID] = (rho_K*E_K + p_K)*u_K;  G_K[rhoEID] = (rho_K*E_K + p_K)*v_K;
  
	/// interaction IJ
	const CFreal vI_n_IJ = u_I*n_IJ[XX] + v_I*n_IJ[YY];
	const CFreal vJ_n_IJ = u_J*n_IJ[XX] + v_J*n_IJ[YY];

	const CFreal SLeft_IJ  = std::min(vI_n_IJ - A_I, vJ_n_IJ - A_J); 
	//vI_n_IJ - A_I; 
	const CFreal SRight_IJ = std::max(vI_n_IJ + A_I, vJ_n_IJ + A_J);
	//vJ_n_IJ + A_J;
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

	/// interaction IK
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
	
    residual[iState] = Fn_IJ + Fn_IK;
  }
}

//////////////////////////////////////////////////////////////////////////////

void RDHLLSchemeCSys::distributePart(vector<RealVector>& residual)
{
  throw Common::NotImplementedException (FromHere(),"void RDHLLSchemeCSys::distributePart()");
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

void RDHLLSchemeCSys::computePicardJacob(vector<RealMatrix*>& jacob)
{
  throw Common::NotImplementedException (FromHere(),"void RDHLLSchemeCSys::computePicardJacob()");
//   // carefull with the signs !!!
//   for (CFuint iState = 0; iState < _nbStatesInCell; ++iState) {
//     const CFuint nStart = iState*_nbStatesInCell;
//     for (CFuint jState = 0; jState < _nbStatesInCell; ++jState) {
//       RealMatrix *const block = jacob[nStart + jState];
// 
//       _tempMat = _invK*(*_kMin[jState]);
//       if (iState == jState) {
// 	for (CFuint iEq = 0; iEq < _nbEquations; ++iEq) {
// 	  _tempMat(iEq,iEq) -= 1.0;
// 	}
//       }
// 
//       _tempMat *= -1.0;
// 
//       (*block) = (*_kPlus[iState])*_tempMat;
//     }
//   }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
