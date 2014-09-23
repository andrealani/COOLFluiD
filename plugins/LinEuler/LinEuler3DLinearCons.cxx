#include <numeric>

#include "LinEuler/LinearizedEuler.hh"
#include "LinEuler3DLinearCons.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/State.hh"
#include "LinEulerPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LinEuler3DLinearCons, JacobianLinearizer, LinearizedEulerModule, 1> lineuler3DLinearConsProvider("LinEuler3DLinearCons");

//////////////////////////////////////////////////////////////////////////////

LinEuler3DLinearCons::LinEuler3DLinearCons(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
 _model(model->getImplementor()->
	 getConvectiveTerm().d_castTo<LinEulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////
LinEuler3DLinearCons::~LinEuler3DLinearCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinEuler3DLinearCons::linearize(const vector<State*>& statesInCell)
{
  ///@note Mean flow variables, they are constant in time and averaged to the cell center

  RealVector& linearData = _model->getPhysicalData();
  const CFuint nbStates = statesInCell.size();

  _sumZ = 0.;

  for (CFuint iEq = 0; iEq < m_physmodel->getNbEq(); ++iEq) {
    for (CFuint iState = 0; iState < nbStates; ++iState) {
       CFuint IDstate = statesInCell[iState]->getLocalID();
       RealVector meanflow_state = _model->getMeanFlowState(IDstate);
      _sumZ[iEq] += meanflow_state[iEq];
    }
  }
  _avZ = _sumZ/static_cast<CFreal>(nbStates);

  linearData[LinEulerTerm::GAMMA]  = _model->getgamma();
  linearData[LinEulerTerm::rho0]   = _avZ[0];
  linearData[LinEulerTerm::U0]     = _avZ[1];
  linearData[LinEulerTerm::V0]     = _avZ[2];
  linearData[LinEulerTerm::W0]     = _avZ[3];
  linearData[LinEulerTerm::P0]     = _avZ[4];
  linearData[LinEulerTerm::c]      = sqrt(linearData[LinEulerTerm::GAMMA]*linearData[LinEulerTerm::P0]/linearData[LinEulerTerm::rho0]);

}

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearizedEuler
} // namespace Physics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
