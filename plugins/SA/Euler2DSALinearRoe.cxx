#include "SA.hh"
#include "Euler2DSALinearRoe.hh"
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

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DSALinearRoe, JacobianLinearizer, SAModule, 1>
euler2DSALinearRoeProvider("Euler2DSALinearRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DSALinearRoe::Euler2DSALinearRoe(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
 _model(model->getImplementor()->getConvectiveTerm().
	d_castTo<EulerSATerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DSALinearRoe::~Euler2DSALinearRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSALinearRoe::linearize(const vector<State*>& statesInCell)
{
  RealVector& linearData = _model->getPhysicalData();

  const CFuint nbStates = statesInCell.size();
  // linearizer does this (numerics dependent)

  _sumZ = 0.;
  for (CFuint iEq = 0; iEq < PhysicalModelStack::getActive()->getNbEq(); ++iEq) {
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      _sumZ[iEq] += (*statesInCell[iState])[iEq];
    }
  }

  cf_assert(_avZ.size() == _sumZ.size());
  _avZ = _sumZ/static_cast<CFreal>(nbStates);

  const CFreal overSqrtRho = 1./_avZ[0];


  linearData[EulerTerm::RHO] = _avZ[0]*_avZ[0];
  linearData[EulerTerm::VX] = _avZ[1]*overSqrtRho;
  linearData[EulerTerm::VY] = _avZ[2]*overSqrtRho;
  linearData[EulerTerm::H]   = _avZ[3]*overSqrtRho;
  const CFreal speed2 = (linearData[EulerTerm::VX]*
				  linearData[EulerTerm::VX] +
				  linearData[EulerTerm::VY]*
				  linearData[EulerTerm::VY]);

  linearData[EulerTerm::V] = std::sqrt(speed2);
  const CFreal a2 = (_model->getGamma()-1.)*(linearData[EulerTerm::H] - 0.5*speed2);

  /// @todo can this make a difference in speed? test !!
  if (a2 < 0.0) {
    std::string msg =
      std::string("Euler2DSAConsLinearRoe::linearize() : a2 < 0 => a2 = ") + Common::StringOps::to_str(a2);
    throw Common::BadValueException (FromHere(),msg);
  }

  linearData[EulerTerm::A] = sqrt(a2);

  const CFuint iK = _model->getFirstScalarVar(0);
  linearData[iK] = _avZ[4]*overSqrtRho;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
