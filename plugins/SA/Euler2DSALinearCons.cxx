#include "SA.hh"
#include "Euler2DSALinearCons.hh"
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

Environment::ObjectProvider<Euler2DSALinearCons, JacobianLinearizer, SAModule, 1>
euler2DSALinearConsProvider("Euler2DSALinearCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DSALinearCons::Euler2DSALinearCons(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
 _model(model->getImplementor()->getConvectiveTerm().
	d_castTo<EulerSATerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DSALinearCons::~Euler2DSALinearCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSALinearCons::linearize(const vector<State*>& statesInCell)
{
 RealVector& linearData = _model->getPhysicalData();

  const CFuint nbStates = statesInCell.size();
  // linearizer does this (numerics dependent)
  _sumZ = 0.;
  const CFuint nbeqs = m_physmodel->getNbEq();
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    for (CFuint iEq = 0; iEq < nbeqs ; ++iEq) {
      _sumZ[iEq] += (*statesInCell[iState])[iEq];
    }
  }
  _avZ = _sumZ/static_cast<CFreal>(nbStates);
  //////

  linearData[EulerTerm::RHO] = _avZ[0];
  linearData[EulerTerm::VX] = _avZ[1]/_avZ[0];
  linearData[EulerTerm::VY] = _avZ[2]/_avZ[0];
  linearData[EulerTerm::V] = sqrt(linearData[EulerTerm::VX]*
				  linearData[EulerTerm::VX] +
				  linearData[EulerTerm::VY]*
				  linearData[EulerTerm::VY]);

  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal uuvv = linearData[EulerTerm::V]*linearData[EulerTerm::V];
  linearData[EulerTerm::H] = gamma*_avZ[3]/linearData[EulerTerm::RHO] -
    0.5*gammaMinus1*uuvv;

  const CFreal a2 = gammaMinus1*(linearData[EulerTerm::H] - 0.5*uuvv);
  linearData[EulerTerm::A] = sqrt(a2);
  linearData[EulerTerm::T] = a2/(_model->getGamma()*_model->getR());
  linearData[EulerTerm::P] = linearData[EulerTerm::RHO]*_model->getR()*linearData[EulerTerm::T];
  linearData[EulerTerm::E] = linearData[EulerTerm::H] -
    linearData[EulerTerm::P]/linearData[EulerTerm::RHO];
  linearData[EulerTerm::GAMMA] = _model->getGamma();
  
  if (a2 < 0.0) {
    std::string msg =
	std::string("Euler2DSAConsLinearCons::linearize() : a2 < 0 => a2 = ") + Common::StringOps::to_str(a2);
    throw Common::BadValueException (FromHere(),msg);
  }
  
  const CFuint iK = _model->getFirstScalarVar(0);
  linearData[iK] = linearData[EulerTerm::RHO]*_avZ[4];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
