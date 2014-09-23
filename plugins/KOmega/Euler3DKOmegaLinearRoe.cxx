#include "KOmega.hh"
#include "Euler3DKOmegaLinearRoe.hh"
#include "Environment/ObjectProvider.hh"

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

Environment::ObjectProvider<Euler3DKOmegaLinearRoe, JacobianLinearizer, KOmegaModule, 1>
euler3DKOmegaLinearRoeProvider("Euler3DKOmegaLinearRoe");

//////////////////////////////////////////////////////////////////////////////

Euler3DKOmegaLinearRoe::Euler3DKOmegaLinearRoe(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->getConvectiveTerm().
	d_castTo<EulerKOmegaTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DKOmegaLinearRoe::~Euler3DKOmegaLinearRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DKOmegaLinearRoe::linearize(const vector<State*>& statesInCell)
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

  linearData[EulerTerm::RHO] = _avZ[0]*_avZ[0];
  linearData[EulerTerm::VX]  = _avZ[1]/_avZ[0];
  linearData[EulerTerm::VY]  = _avZ[2]/_avZ[0];
  linearData[EulerTerm::VZ]  = _avZ[3]/_avZ[0];
  linearData[EulerTerm::H]   = _avZ[4]/_avZ[0];
  const CFreal speed2 = (linearData[EulerTerm::VX]*
                         linearData[EulerTerm::VX] +
                         linearData[EulerTerm::VY]*
                         linearData[EulerTerm::VY] +
                         linearData[EulerTerm::VZ]*
                         linearData[EulerTerm::VZ]);

  linearData[EulerTerm::V] = std::sqrt(speed2);
  const CFreal a2 = (_model->getGamma()-1.)*(linearData[EulerTerm::H] - 0.5*speed2);

  /// @todo can this make a difference in speed? test !!
  if (a2 < 0.0) {
    std::string msg =
      std::string("Euler3DKOmegaLinearRoe::linearize() : a2 < 0 => a2 = ")
      + Common::StringOps::to_str(a2);
    throw Common::BadValueException (FromHere(),msg);
  }

  linearData[EulerTerm::A] = sqrt(a2);

  const CFuint iK  = _model->getFirstScalarVar(0);
  linearData[iK]   = _avZ[5]/_avZ[0];
  linearData[iK+1] = _avZ[6]/_avZ[0];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
