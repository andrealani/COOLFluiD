#include "MultiFluidMHD/MultiFluidMHD.hh"
#include "Euler3DMFMHDConsToRhoiViTiInRhoiViTi.hh"
#include "MultiFluidMHDModel.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "EulerMFMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler3DMFMHDConsToRhoiViTiInRhoiViTi, VarSetMatrixTransformer, 
			    MultiFluidMHDModule, 1> 
euler3DConsToRhovtInRhovtProvider("Euler3DMFMHDConsToRhoiViTiInRhoiViTi");

//////////////////////////////////////////////////////////////////////////////

Euler3DMFMHDConsToRhoiViTiInRhoiViTi::Euler3DMFMHDConsToRhoiViTiInRhoiViTi(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerMFMHDTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DMFMHDConsToRhoiViTiInRhoiViTi::~Euler3DMFMHDConsToRhoiViTiInRhoiViTi()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DMFMHDConsToRhoiViTiInRhoiViTi::setMatrix(const RealVector& state)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DMFMHDConsToRhoiViTiInRhoiViTi::setMatrix()");
//   const CFuint nbSpecies    = _model->getNbScalarVars(0);
//   const CFuint nbMomentum   = _model->getNbScalarVars(1);
//   const CFuint nbEnergyEqs  = _model->getNbScalarVars(2);
//   const CFuint endEM = 8;
//   
//   // Electro magnetic field needs no transformation
//   CFuint endEM = 8;
//   for (CFuint i = 0; i < endEM; ++i){
//     _transMatrix(i,i) = 1.0;
//   }
//   
//   // Transformation from conservative to Rhoi Vi Ti variables in 3D
//   // First transformation of the densities
//   for (CFuint ie = 0; ie < nbSpecies; ++ie) {
//     _transMatrix(ie,ie) = 1.0;
//   }
//   
//   //Transformation taking as example Euler3DConsToRhovtInRhovt
//   for (CFuint ie = 0; ie < nbSpecies; ++ie) {
//     
//     const CFreal R = _model->getR();
//     const CFreal rho = state[0];
//     const CFreal invRho = 1/rho;
//     const CFreal u = state[1];
//     const CFreal v = state[2];
//     const CFreal T = state[3];
//     const CFreal V2 = u*u + v*v;
//     const CFreal gammaMinus1 = _model->getGamma() - 1.;
//     const CFreal kg = gammaMinus1/(R*rho);
//   
// 
//     _transMatrix(1,0) = -u*invRho;
//     _transMatrix(1,1) = invRho;
//   
//     _transMatrix(2,0) = -v*invRho;
//     _transMatrix(2,2) = invRho;
//   
//     _transMatrix(3,0) = kg*(0.5*V2 - R*T/gammaMinus1);
//     _transMatrix(3,1) = -kg*u;
//     _transMatrix(3,2) = -kg*v;
//     _transMatrix(3,3) = kg; 

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
