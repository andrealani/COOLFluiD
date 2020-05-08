#include "KOmega.hh"
#include "Euler2DKLogOmegaConsToPuvtInPuvt.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DKLogOmegaConsToPuvtInPuvt, VarSetMatrixTransformer,
	       KOmegaModule, 1>
euler2DKLogOmegaConsToPuvtInPuvt("Euler2DKLogOmegaConsToPuvtInPuvt");

//////////////////////////////////////////////////////////////////////////////

Euler2DKLogOmegaConsToPuvtInPuvt::Euler2DKLogOmegaConsToPuvtInPuvt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerKLogOmegaTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DKLogOmegaConsToPuvtInPuvt::~Euler2DKLogOmegaConsToPuvtInPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DKLogOmegaConsToPuvtInPuvt::setMatrix(const RealVector& state)
{
  const CFreal R = _model->getR();
  const CFreal p = state[0];
  const CFreal invRho = R*state[3]/p;
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal gammaMinus1 = _model->getGamma() - 1.;

  _transMatrix(0,0) = (u*u + v*v)*0.5*gammaMinus1;
  _transMatrix(0,1) = -gammaMinus1*u;
  _transMatrix(0,2) = -gammaMinus1*v;
  _transMatrix(0,3) = gammaMinus1;

  _transMatrix(1,0) = -u*invRho;
  _transMatrix(1,1) = invRho;

  _transMatrix(2,0) = -v*invRho;
  _transMatrix(2,2) = invRho;

  _transMatrix(3,0) = (_transMatrix(0,0) -
		       p*invRho)*invRho/R;
  _transMatrix(3,1) = _transMatrix(0,1)*invRho/R;
  _transMatrix(3,2) = _transMatrix(0,2)*invRho/R;
  _transMatrix(3,3) = _transMatrix(0,3)*invRho/R;

  _transMatrix(4,4) = invRho;

  _transMatrix(5,5) = invRho;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
