#include "GReKO.hh"
#include "Euler2DGReKOLinearRoe.hh"
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

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DGReKOLinearRoe, JacobianLinearizer, GReKOModule, 1>
euler2DGReKOLinearRoeProvider("Euler2DGReKOLinearRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DGReKOLinearRoe::Euler2DGReKOLinearRoe(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->getConvectiveTerm().
	d_castTo<EulerGReKOTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DGReKOLinearRoe::~Euler2DGReKOLinearRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOLinearRoe::linearize(const vector<State*>& statesInCell)
{
  RealVector& linearData = _model->getPhysicalData();

  const CFuint nbStates = statesInCell.size();

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
  linearData[EulerTerm::VX] = _avZ[1]/_avZ[0];
  linearData[EulerTerm::VY] = _avZ[2]/_avZ[0];
  linearData[EulerTerm::H]   = _avZ[3]/_avZ[0];
  const CFreal speed2 = (linearData[EulerTerm::VX]*
				  linearData[EulerTerm::VX] +
				  linearData[EulerTerm::VY]*
				  linearData[EulerTerm::VY]);

  linearData[EulerTerm::V] = std::sqrt(speed2);
  const CFreal a2 = (_model->getGamma()-1.)*(linearData[EulerTerm::H] - 0.5*speed2);

  /// @todo can this make a difference in speed? test !!
  if (a2 < 0.0) {
    std::string msg =
      std::string("Euler2DGReKOLinearRoe::linearize() : a2 < 0 => a2 = ")
      + Common::StringOps::to_str(a2);
    throw Common::BadValueException (FromHere(),msg);
  }

  linearData[EulerTerm::A] = sqrt(a2);

  const CFuint iK = _model->getFirstScalarVar(0);
  linearData[iK] = _avZ[4]/_avZ[0];
  linearData[iK+1] = _avZ[5]/_avZ[0];
  linearData[iK+2] = _avZ[6]/_avZ[0];
  linearData[iK+3] = _avZ[7]/_avZ[0];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
