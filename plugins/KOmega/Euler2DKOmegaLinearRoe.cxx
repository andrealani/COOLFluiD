#include "KOmega.hh"
#include "Euler2DKOmegaLinearRoe.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DKOmegaLinearRoe, JacobianLinearizer, KOmegaModule, 1>
euler2DKOmegaLinearRoeProvider("Euler2DKOmegaLinearRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DKOmegaLinearRoe::Euler2DKOmegaLinearRoe(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->getConvectiveTerm().
	d_castTo<EulerKOmegaTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DKOmegaLinearRoe::~Euler2DKOmegaLinearRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DKOmegaLinearRoe::linearize(const vector<State*>& statesInCell)
{
  RealVector& linearData = _model->getPhysicalData();

  const CFuint nbStates = statesInCell.size();
  const CFuint iK = _model->getFirstScalarVar(0);
  // linearizer does this (numerics dependent)
  _sumZ = 0.;
  
  const CFuint nbeqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    State& stateRoe = *statesInCell[iState];
    for (CFuint iEq = 0; iEq < nbeqs; ++iEq) {
      _sumZ[iEq] += stateRoe[iEq];
    }
  }
  cf_assert(_avZ.size() == _sumZ.size());
  _avZ = _sumZ/static_cast<CFreal>(nbStates);

  linearData[EulerTerm::RHO] = _avZ[0]*_avZ[0];
  linearData[EulerTerm::VX]  = _avZ[1]/_avZ[0];
  linearData[EulerTerm::VY]  = _avZ[2]/_avZ[0];
  linearData[EulerTerm::H]   = _avZ[3]/_avZ[0];
  linearData[iK] = _avZ[4]/_avZ[0];
  linearData[iK+1] = _avZ[5]/_avZ[0];
  const CFreal speed2 = (linearData[EulerTerm::VX]*
				  linearData[EulerTerm::VX] +
				  linearData[EulerTerm::VY]*
				  linearData[EulerTerm::VY]);

  linearData[EulerTerm::V] = std::sqrt(speed2);
   CFreal a2 = (_model->getGamma()-1.)*(linearData[EulerTerm::H] - 0.5*speed2 - linearData[iK] );

  /// @todo can this make a difference in speed? test !!
  if (a2 < 0.0) {
    std::string msg =
      std::string("Euler2DKOmegaLinearRoe::linearize() : a2 < 0 => a2 = ")
      + Common::StringOps::to_str(a2);
    throw Common::BadValueException (FromHere(),msg);
  }

  linearData[EulerTerm::A] = sqrt(a2);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
