#include "NavierStokes/NavierStokes.hh"
#include "Euler2DLinearRoe.hh"
#include "Environment/ObjectProvider.hh"
#include "EulerPhysicalModel.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DLinearRoe,
			    JacobianLinearizer,
			    NavierStokesModule, 1>
euler2DLinearRoeProvider("Euler2DLinearRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DLinearRoe::Euler2DLinearRoe(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
 _model(model->getImplementor()->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DLinearRoe::~Euler2DLinearRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLinearRoe::linearize(const vector<State*>& statesInCell)
{
  RealVector& linearData = _model->getPhysicalData();
  const CFuint nbStates = statesInCell.size();
  // linearizer does this (numerics dependent)
  
  _sumZ = 0.;
  const CFreal invNbStates = 1./static_cast<CFreal>(nbStates);
  const CFuint nbeqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    State& stateRoe = *statesInCell[iState];
    for (CFuint iEq = 0; iEq < nbeqs; ++iEq) {
      _sumZ[iEq] += stateRoe[iEq];
    }
  }
  
  cf_assert(_avZ.size() == _sumZ.size());
  _avZ = _sumZ*invNbStates;
  //////
  
  const CFreal ovAvZ0 = 1./_avZ[0];
  linearData[EulerTerm::RHO] = _avZ[0]*_avZ[0];
  linearData[EulerTerm::VX] = _avZ[1]*ovAvZ0;
  linearData[EulerTerm::VY] = _avZ[2]*ovAvZ0;
  const CFreal speed2 = (linearData[EulerTerm::VX]*
			 linearData[EulerTerm::VX] +
			 linearData[EulerTerm::VY]*
			 linearData[EulerTerm::VY]);
  
  linearData[EulerTerm::V] = std::sqrt(speed2);
  linearData[EulerTerm::H]  = _avZ[3]*ovAvZ0;
  
  const CFreal gamma = _model->getGamma();
  const CFreal a2 = (gamma-1.)*(linearData[EulerTerm::H] - 0.5*speed2);
  
  if (a2 < 0.0)
  {
    CFLogWarn("Linearized state gives a2 < 0. LinearData [" << linearData << "]\n");
    
    std::string msg = "Euler2DConsLinearRoe::linearize() : a2 < 0 => a2 = " + Common::StringOps::to_str(a2);
    throw Common::BadValueException (FromHere(),msg);
  }
  
  linearData[EulerTerm::A] = sqrt(a2);
  
  const CFreal Rgas = _model->getR();
  linearData[EulerTerm::T] = a2/(gamma*Rgas);
  const CFreal RT = Rgas*linearData[EulerTerm::T];
  linearData[EulerTerm::P] = linearData[EulerTerm::RHO]*RT;
  linearData[EulerTerm::E] = linearData[EulerTerm::H] - RT;
  linearData[EulerTerm::GAMMA] = gamma;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
