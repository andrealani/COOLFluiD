#include "NavierStokes/NavierStokes.hh"
#include "Euler2DCharToConsInRef.hh"
#include "EulerPhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DCharToConsInRef, VarSetMatrixTransformer, NavierStokesModule,1> euler2DCharToConsInRefProvider("Euler2DCharToConsInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DCharToConsInRef::Euler2DCharToConsInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DCharToConsInRef::~Euler2DCharToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DCharToConsInRef::setMatrixFromRef()
{
  cf_assert(_model.isNotNull());

  const RealVector& linearData = _model->getPhysicalData();

  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU = linearData[EulerTerm::VX];
  const CFreal avV = linearData[EulerTerm::VY];
  const CFreal speed = sqrt(avU*avU + avV*avV);
  const CFreal cosD = avU/speed;
  const CFreal sinD = avV/speed;
  const CFreal avA   = linearData[EulerTerm::A];
  const CFreal M     = speed/avA;
  const CFreal M2    = M*M;
  const CFreal invM  = 1./M;
  const CFreal eps   = 0.05;
  const CFreal beta  = sqrt(max(eps*eps,std::abs(M2 - 1.)));
  const CFreal chi   = beta/max(M,1.);
  const CFreal rhoDivA = avRho/avA;
  const CFreal invA2   = 1./(avA*avA);
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal k1  = 0.5*beta/(chi*M2);
  const CFreal k2  = k1*rhoDivA;
  const CFreal k3  = 0.5*avRho/(chi*M);
  const CFreal k4  = rhoDivA*speed*speed*0.5 + avRho*avA/gammaMinus1;
  const CFreal k5  = -avU*sinD + avV*cosD;
  const CFreal k6  = avU*cosD + avV*sinD;

  _transMatrix(0,0) = k2;
  _transMatrix(0,1) = k2;
  _transMatrix(0,2) = rhoDivA/M2;
  _transMatrix(0,3) = -invA2;

  _transMatrix(1,0) = k2*avU - k3*sinD;
  _transMatrix(1,1) = k2*avU + k3*sinD;
  _transMatrix(1,2) = rhoDivA*avU/M2 + invM*avRho*cosD;
  _transMatrix(1,3) = -invA2*avU;

  _transMatrix(2,0) = k2*avV + k3*cosD;
  _transMatrix(2,1) = k2*avV - k3*cosD;
  _transMatrix(2,2) = rhoDivA*avV/M2 + invM*avRho*sinD;
  _transMatrix(2,3) = -invA2*avV;

  _transMatrix(3,0) = k1*k4 + k3*k5;
  _transMatrix(3,1) = k1*k4 - k3*k5;
  _transMatrix(3,2) = k4/M2 + k6*avRho*invM;
  _transMatrix(3,3) = -0.5*M2;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
