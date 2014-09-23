#include "KOmega.hh"
#include "Euler2DKOmegaLinearCons.hh"
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

Environment::ObjectProvider<Euler2DKOmegaLinearCons, JacobianLinearizer, KOmegaModule, 1>
euler2DKOmegaLinearConsProvider("Euler2DKOmegaLinearCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DKOmegaLinearCons::Euler2DKOmegaLinearCons(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->getConvectiveTerm().
	d_castTo<EulerKOmegaTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DKOmegaLinearCons::~Euler2DKOmegaLinearCons()
{
  
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DKOmegaLinearCons::linearize(const vector<State*>& statesInCell)
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

  linearData[EulerTerm::RHO] = _avZ[0];
  linearData[EulerTerm::VX] = _avZ[1]/_avZ[0];
  linearData[EulerTerm::VY] = _avZ[2]/_avZ[0];
  linearData[EulerTerm::V] = sqrt(linearData[EulerTerm::VX]*
				  linearData[EulerTerm::VX] +
				  linearData[EulerTerm::VY]*
				  linearData[EulerTerm::VY]);
  
  linearData[iK] = _avZ[4];
  linearData[iK+1] = _avZ[5];

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
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
