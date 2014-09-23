#include "LTE.hh"
#include "Euler2DLinearPrvtLTE.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/State.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DLinearPrvtLTE, JacobianLinearizer, LTEModule, 1> euler2DLinearPrvtLTEProvider("Euler2DLinearPrvtLTE");

//////////////////////////////////////////////////////////////////////////////

Euler2DLinearPrvtLTE::Euler2DLinearPrvtLTE(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->
	 getConvectiveTerm().d_castTo<EulerTerm>()),
  _dhe(3)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DLinearPrvtLTE::~Euler2DLinearPrvtLTE()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLinearPrvtLTE::linearize(const vector<State*>& statesInCell)
{
  // statesInCell is an array of state vectors [p rhou rhov T]

  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  RealVector& linearData = _model->getPhysicalData();
  const RealVector& refData = _model->getReferencePhysicalData();

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
  
  // set the composition
  CFreal pdim = _model->getPressureFromState(_avZ[0])*refData[EulerTerm::P];
  CFreal Tdim = _avZ[3]*_model->getTempRef();
  library->setComposition(Tdim,pdim);
  // set the density, enthalpy and internal energy
  library->setDensityEnthalpyEnergy(Tdim, pdim, _dhe);
  
  linearData[EulerTerm::RHO] = _dhe[0]/refData[EulerTerm::RHO];
  linearData[EulerTerm::P] = _avZ[0];
  linearData[EulerTerm::VX] = _avZ[1]/linearData[EulerTerm::RHO];
  linearData[EulerTerm::VY] = _avZ[2]/linearData[EulerTerm::RHO];
  linearData[EulerTerm::V] = sqrt(linearData[EulerTerm::VX]*
				  linearData[EulerTerm::VX] +
				  linearData[EulerTerm::VY]*
				  linearData[EulerTerm::VY]);

  const CFreal halfV2 = 0.5*linearData[EulerTerm::V]*linearData[EulerTerm::V];
  linearData[EulerTerm::H] = (_dhe[1] + halfV2)/refData[EulerTerm::H];
  linearData[EulerTerm::E] = (_dhe[2] + halfV2)/refData[EulerTerm::E];
  
  library->gammaAndSoundSpeed(Tdim, pdim, _dhe[0], 
			      linearData[EulerTerm::GAMMA], 
			      linearData[EulerTerm::A]);
  linearData[EulerTerm::T] = _avZ[3];
  
  // compute and store also derivatives in linearData
  // dpdT,drdT,dhdT,dedT
  // dpdP,drdP,dhdP,dedP
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LTE
} // namespace Physics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
