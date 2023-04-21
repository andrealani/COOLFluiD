#include "MHD/MHD.hh"
#include "MHD/MHD3DProjectionDiffPrim.hh"
#include "Environment/ObjectProvider.hh"
#include "MHD/MHDProjectionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHD3DProjectionDiffPrim, DiffusiveVarSet, MHDModule, 2> 
mhd3DProjectionDiffPrimProvider("MHD3DProjectionDiffPrim");
      
//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionDiffPrim::MHD3DProjectionDiffPrim(const std::string& name,
						 Common::SafePtr<Framework::PhysicalModelImpl> model) :
  MHD3DProjectionDiffVarSet(name, model),
  _eulerModel(model->getConvectiveTerm().d_castTo<MHDProjectionTerm>()),
  _tempX()
{
  vector<std::string> names(9);
  names[0] = "rho";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "Bx";
  names[5] = "By";
  names[6] = "Bz";
  names[7] = "p";
  names[8] = "phi";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionDiffPrim::~MHD3DProjectionDiffPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffPrim::setGradientVars(const vector<RealVector*>& states,
					      RealMatrix& values,
					      const CFuint stateSize)
{
  /// PETER
  /*
  Here you have to compute multiple states [rho u v w Bx By Bz T phi] from the given states [rho u v w Bx By Bz p phi]
states and values are matrices storing multiple states
  */


CFreal mu = 1.27;       // Mean molecular weight
CFreal mH = 1.67e-27;   // Mass hydrogen
CFreal kB = 1.38e-23;

// First copy over all states
for (CFuint i = 0; i < 9 ; ++i) {
    for (CFuint j = 0; j < 6 ; ++j) {

      values(i,j) = (*states[j])[i];

      //std::cout << "i = " << i << "\n";
      //std::cout << "j = " << j << "\n";
    }
}

// Now overwrite values(7,j):
for (CFuint j = 0; j < 6; ++j) {
   //values(7,j) = (*states[j])[7]*mu*mH/(2*(*states[j])[0]*kB;
   // Adimensional:
   values(7,j) = (*states[j])[7]/(*states[j])[0]; // T = P[adim]/rho[adim]
}







/*
for (CFuint i = 0; i < nbValues; ++i) {
    for (CFuint j = 0; j < stateSize; ++j) {
      
      if (j==7) {
       values(i,j) = (*states[j])[i]*mu*mH/(2*(*states[0])[i]*kB);
      } else {
        values(i,j) = (*states[j])[i];
      }

    }
}
*/

  /// examples for Navier-Stokes
  // from [p u v w T] to [p u v w T] 
  /*for (CFuint i = 0; i < nbValues; ++i) {
    for (CFuint j = 0; j < stateSize; ++j) {
      values(i,j) = (*states[j])[i];
    }
    }*/

  // example
  // from [rho rhoU rhoV rhoW rhoE] to [p u v w T] 
  /*const CFreal R = _eulerModel->getR();
  const CFreal ovCv = (_eulerModel->getGamma() - 1.)/R;
  
  for (CFuint i = 0; i < stateSize; ++i) {
    const RealVector& state = *states[i];
    const CFreal ovRho = 1./state[0]; 
    values(1,i) = state[1]*ovRho;
    values(2,i) = state[2]*ovRho;
    values(3,i) = state[3]*ovRho;
    
    const CFreal V2 = values(1,i)*values(1,i) + 
      values(2,i)*values(2,i) + values(3,i)*values(3,i);
    values(4,i) = (state[4]*ovRho - 0.5*V2)*ovCv;
    values(0,i) = R*state[0]*values(4,i);
    }*/
  
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffPrim::setGradientVarGradients(const vector<RealVector*>& states,
						      const vector< vector<RealVector*> >& stateGradients,
						      vector< vector<RealVector*> >& gradVarGradients,
						      const CFuint stateSize)
{
 throw Common::NotImplementedException
   (FromHere(), "MHD3DProjectionDiffPrim::setGradientVarGradients()");
}
      
//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffPrim::setStateGradients(const vector<RealVector*>& states,
						const vector< vector<RealVector*> >& gradVarGradients,
						vector< vector<RealVector*> >& stateGradients,
						const CFuint stateSize)
{
  /// PETER
  throw Common::NotImplementedException
   (FromHere(), "MHD3DProjectionDiffPrim::setStateGradients()");
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHD3DProjectionDiffPrim::getDynViscosity(const RealVector& state,
						const vector<RealVector*>& gradients)
{
  throw Common::NotImplementedException
    (FromHere(), "MHD3DProjectionDiffPrim::getDynViscosity()");
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHD3DProjectionDiffPrim::getDensity(const RealVector& state)
{  
  return state[0];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffPrim::setGradientState(const RealVector& state)
{
  /// PETER
  cf_assert(_gradState.size() == state.size());

  /// examples for Navier-Stokes
  /// from [p u v w T] to [p u v w T] 
  // _gradState = state;

  /// from [rho rhoU rhoV rhoW rhoE] to [p u v w T] 
  /*
    const CFreal R = _eulerModel->getR();
  const CFreal cv = R/(_eulerModel->getGamma() - 1.);

  // _gradState = [p u v w T]
  _gradState[1] = state[1]/state[0];
  _gradState[2] = state[2]/state[0];
  _gradState[3] = state[3]/state[0];
  const CFreal V2 = _gradState[1]*_gradState[1] +
    _gradState[2]*_gradState[2] +
    _gradState[3]*_gradState[3];

  _gradState[4] = (state[4] - 0.5*state[0]*V2)/(state[0]*cv);
  _gradState[0] = R*state[0]*_gradState[4];
  */

CFreal mu = 1.27;       // Mean molecular weight
CFreal mH = 1.67e-27;   // Mass hydrogen
CFreal kB = 1.38e-23;

 // Again, copy over the full state vector:
 _gradState = state;
 // And then overwrite _gratState[7] = T
 //_gradState[7] = state[7]*mu*mH/(2.0*state[0]*kB);
 // Adimensional:
 _gradState[7] = state[7]/state[0];

}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffPrim::computeFluxJacobian(const RealVector& state,
						  const RealVector& gradientJacob,
						  const RealVector& normal,
						  const CFreal& radius,
						  RealMatrix& fluxJacob)
{
  throw Common::NotImplementedException
    (FromHere(), "MHD3DProjectionDiffPrim::computeFluxJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffPrim::setComposition(const RealVector& state,
					     const bool isPerturb,
					     const CFuint iVar)
{
}

//////////////////////////////////////////////////////////////////////////////

} // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
