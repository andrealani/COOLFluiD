#include "MHD/MHDProjectionDiffVarSet.hh"
#include "MHD/MHDProjectionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////
       
MHDProjectionDiffVarSet::MHDProjectionDiffVarSet(const std::string& name, 
				       SafePtr<PhysicalModelImpl> model) :
  DiffusiveVarSet(name, model),
  _model(model->getDiffusiveTerm().d_castTo<MHDProjectionDiffTerm>()),
  _twoThird(2./3.),
  _dynViscCoeff(0.0),
  _thermCondCoeff(0.0),
  _useBackUpValues(false),
  _setBackUpValues(false),
  _gradState(),
  _normal(),
  _tau(),
  _uID(0),
  _vID(0),
  _wID(0),
  _TID(0)
{
}
  
//////////////////////////////////////////////////////////////////////////////
      
MHDProjectionDiffVarSet::~MHDProjectionDiffVarSet()
{
}
    
//////////////////////////////////////////////////////////////////////////////
   
void MHDProjectionDiffVarSet::setup()
{
  DiffusiveVarSet::setup();
  _gradState.resize(PhysicalModelStack::getActive()->getNbEq());
 
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _normal.resize(dim);
  
  _tau.resize(dim,dim);
  
  _uID = 1;
  _vID = 2;
  _wID = 3;
  _TID = 7;
}

//////////////////////////////////////////////////////////////////////////////

void MHDProjectionDiffVarSet::computeTransportProperties(const RealVector& state,
							 const vector<RealVector*>& gradients,
							 const RealVector& normal)
{
  /*RealVector& nsData = getModel().getPhysicalData();
  
  // adimensional dynamical viscosity
  if (_useBackUpValues || _freezeDiffCoeff) {
    nsData[MHDProjectionDiffTerm::MU] = _dynViscCoeff;
    nsData[MHDProjectionDiffTerm::LAMBDA] = _thermCondCoeff;
  }
  else {
    // adimensional dynamical viscosity
    nsData[MHDProjectionDiffTerm::MU] = getDynViscosity(state, gradients);
    
    // adimensional thermal conductivity
    nsData[MHDProjectionDiffTerm::LAMBDA] = getThermConductivity(state, nsData[MHDProjectionDiffTerm::MU]);
    
    if (_setBackUpValues) {
      _dynViscCoeff = nsData[MHDProjectionDiffTerm::MU];
      _thermCondCoeff = nsData[MHDProjectionDiffTerm::LAMBDA];
    }
    }*/
}

//////////////////////////////////////////////////////////////////////////////

void MHDProjectionDiffVarSet::computeStressTensor(const RealVector& state,
						  const vector<RealVector*>& gradients,
						  const CFreal& radius) 
{
  /*const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  CFreal divTerm = 0.0;
  if (dim == DIM_2D && radius > MathConsts::CFrealEps()) {
    // if the face is a boundary face, the radius could be 0
    // check against eps instead of 0. for safety
    divTerm = _gradState[_vID]/radius;
  }
  else if (dim == DIM_3D) {
    const RealVector& gradW = *gradients[_wID];
    divTerm = gradW[ZZ];
  }
  
  const RealVector& gradU = *gradients[_uID];
  const RealVector& gradV = *gradients[_vID];
  const CFreal twoThirdDivV = _twoThird*(gradU[XX] + gradV[YY] + divTerm);
  const RealVector& nsData = getModel().getPhysicalData();
  const CFreal coeffTauMu = nsData[MHDProjectionDiffTerm::MU];
  _tau(XX, XX) = coeffTauMu*(2.*gradU[XX] - twoThirdDivV);
  _tau(XX, YY) = _tau(YY, XX) = coeffTauMu*(gradU[YY] + gradV[XX]);
  _tau(YY, YY) = coeffTauMu*(2.*gradV[YY] - twoThirdDivV);
  
  if (dim == DIM_3D) {
    const RealVector& gradW = *gradients[_wID];
    _tau(XX,ZZ) = _tau(ZZ,XX) = coeffTauMu*(gradU[ZZ] + gradW[XX]);
    _tau(YY,ZZ) = _tau(ZZ,YY) = coeffTauMu*(gradV[ZZ] + gradW[YY]);
    _tau(ZZ,ZZ) = coeffTauMu*(2.*gradW[ZZ] - twoThirdDivV);
    }*/
}
      
//////////////////////////////////////////////////////////////////////////////

CFreal MHDProjectionDiffVarSet::getHeatFlux(const RealVector& state,
					    const vector<RealVector*>& gradients,
					    const RealVector& normal)
{
  /// to be modified 
  /*const CFreal lambda = getThermConductivity(state, mu);
  const RealVector& gradT = *gradients[_TID];
  return (-lambda*(MathFunctions::innerProd(gradT,normal)));
  */
  return 0.;
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
