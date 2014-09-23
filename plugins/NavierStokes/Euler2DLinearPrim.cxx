#include "NavierStokes/NavierStokes.hh"
#include "Euler2DLinearPrim.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/State.hh"
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

Environment::ObjectProvider<Euler2DLinearPrim, JacobianLinearizer, NavierStokesModule, 1> 
euler2DLinearPrimProvider("Euler2DLinearPrim");

//////////////////////////////////////////////////////////////////////////////

Euler2DLinearPrim::Euler2DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
 _model(model->getImplementor()->
	 getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DLinearPrim::~Euler2DLinearPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLinearPrim::linearize(const vector<State*>& statesInCell)
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
  //////
  
  linearData[EulerTerm::RHO] = _avZ[0];
  linearData[EulerTerm::VX]  = _avZ[1];
  linearData[EulerTerm::VY]  = _avZ[2];
  
  const CFreal V2 = linearData[EulerTerm::VX]*linearData[EulerTerm::VX] +
    linearData[EulerTerm::VY]*linearData[EulerTerm::VY];
  linearData[EulerTerm::V] = sqrt(V2);
  
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal R =_model->getR();
  
  linearData[EulerTerm::P] = _avZ[3];
  const CFreal p = _model->getPressureFromState(linearData[EulerTerm::P]);
  linearData[EulerTerm::T] = p/(linearData[EulerTerm::RHO]*R); 
  
  const CFreal pOvRho = p/linearData[EulerTerm::RHO];
  linearData[EulerTerm::A] = sqrt(gamma*R*linearData[EulerTerm::T]);
  linearData[EulerTerm::E] = pOvRho/gammaMinus1 + 0.5*V2;
  linearData[EulerTerm::H] = linearData[EulerTerm::E] + pOvRho;
  linearData[EulerTerm::GAMMA] = gamma;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace NavierStokes
} // namespace Physics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
