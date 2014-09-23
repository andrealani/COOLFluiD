#include "NavierStokes/NavierStokes.hh"
#include "Euler3DLinearRoe.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/State.hh"
#include "Common/BadValueException.hh"
#include "EulerPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler3DLinearRoe, JacobianLinearizer, NavierStokesModule, 1> euler3DLinearRoeScalarProvider("Euler3DLinearRoe");

//////////////////////////////////////////////////////////////////////////////

Euler3DLinearRoe::Euler3DLinearRoe(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
 _model(model->getImplementor()->
	 getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DLinearRoe::~Euler3DLinearRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLinearRoe::linearize(const vector<State*>& statesInCell)
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
  _avZ = _sumZ/static_cast<CFreal>(nbStates);
  //

  linearData[EulerTerm::RHO] = _avZ[0]*_avZ[0];
  linearData[EulerTerm::VX] = _sumZ[1]/_sumZ[0];
  linearData[EulerTerm::VY] = _sumZ[2]/_sumZ[0];
  linearData[EulerTerm::VZ] = _sumZ[3]/_sumZ[0];
  linearData[EulerTerm::H]  = _sumZ[4]/_sumZ[0];
  linearData[EulerTerm::V]  = sqrt(linearData[EulerTerm::VX]*
				  linearData[EulerTerm::VX] +
				   linearData[EulerTerm::VY]*
				   linearData[EulerTerm::VY] +
				   linearData[EulerTerm::VZ]*
				   linearData[EulerTerm::VZ]);

  const CFreal speed2 = linearData[EulerTerm::V]*linearData[EulerTerm::V];
  const CFreal a2 = (_model->getGamma() - 1.)*(linearData[EulerTerm::H] - 0.5*speed2);

  if (a2 < 0.0) {
    std::string msg = std::string
      ("Euler3DConsLinearRoe::linearizeT() : a2 < 0 => a2 = ") + Common::StringOps::to_str(a2);
    throw Common::BadValueException (FromHere(),msg);
  }

  linearData[EulerTerm::A] = sqrt(a2);
  linearData[EulerTerm::T] = a2/(_model->getGamma()*_model->getR());
  linearData[EulerTerm::P] = linearData[EulerTerm::RHO]*_model->getR()*linearData[EulerTerm::T];
  linearData[EulerTerm::E] = linearData[EulerTerm::H] -
    linearData[EulerTerm::P]/linearData[EulerTerm::RHO];
  linearData[EulerTerm::GAMMA] = _model->getGamma();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
