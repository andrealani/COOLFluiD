#include "NEQ/NEQ.hh"
#include "Euler2DNEQLinearRoe.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/State.hh"
#include "Framework/MultiScalarTerm.hh"
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

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQLinearRoe, 
			    JacobianLinearizer, NEQModule, 1> 
euler2DNEQLinearRoeProvider("Euler2DNEQLinearRoe");

Environment::ObjectProvider<Euler2DNEQLinearRoe, 
			    JacobianLinearizer, NEQModule, 1> 
euler2DNEQLinearRoeVinokurProvider("Euler2DNEQLinearRoeVinokur");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQLinearRoe::Euler2DNEQLinearRoe(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->
	 getConvectiveTerm().d_castTo< MultiScalarTerm<EulerTerm> >()),
  _ye(),
  _dhe(),
  _tvDim(),
  _evDim()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQLinearRoe::~Euler2DNEQLinearRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQLinearRoe::linearize(const vector<State*>& statesInCell)
{
  // this is an inexact linearization that assumes conservation
  // for all components of the flux except for those related to the momentum
  // equations, where pressure is not assumed quadratic in the Z variables and 
  // computed starting from a simple arithmetic average of temperature
  RealVector& linearData = _model->getPhysicalData();
  cf_assert(_upStates != CFNULL);
  const vector<State*>& upStates = *_upStates;
  
  const CFuint nbStates = statesInCell.size();
  cf_assert(upStates.size() == nbStates);
  
  // linearizer does this (numerics dependent)
  _sumZ = 0.;
  for (CFuint iEq = 0; iEq < PhysicalModelStack::getActive()->getNbEq(); ++iEq) {
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      _sumZ[iEq] += (*statesInCell[iState])[iEq];
    }
  }
  _avZ = _sumZ/static_cast<CFreal>(nbStates);
  //////
  
  // AL: I'm assuming to store T in the update variables
  // (this should be the case  in this kind of application
  // anyway ...)
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFuint TID = nbSpecies + 2;
  CFreal T = 0.0;
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    T += (*upStates[iState])[TID];
  }
  T /= static_cast<CFreal>(nbStates);
  
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  // Set the mixture density (sum of the partial densities)
  CFreal sqRho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    sqRho += _avZ[ie];
  }
  
  const CFreal rho = sqRho*sqRho;
  const CFreal ovSqRho = 1./sqRho;
  _ye.resize(nbSpecies);
  
  // set the species mass fractions
  const CFuint firstSpecies = _model->getFirstScalarVar(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ye[ie] = _avZ[ie]*ovSqRho; 
    linearData[firstSpecies + ie] = _ye[ie];
  }
  
  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!! 
  library->setSpeciesFractions(_ye);
  
  const RealVector& refData =  _model->getReferencePhysicalData();
  CFreal rhoDim = rho*refData[EulerTerm::RHO];
  CFreal Tdim = T*refData[EulerTerm::T];
  
  // this pressure is inconsistent (tvDim required ...)
  CFreal pdim = library->pressure(rhoDim, Tdim, CFNULL);
  CFreal p = pdim/refData[EulerTerm::P];
  
  linearData[EulerTerm::P] = p;
  linearData[EulerTerm::VX] = _avZ[nbSpecies]*ovSqRho;
  linearData[EulerTerm::VY] = _avZ[nbSpecies+1]*ovSqRho;
  
  const CFreal V2 = linearData[EulerTerm::VX]*linearData[EulerTerm::VX] +
    linearData[EulerTerm::VY]*linearData[EulerTerm::VY];
  linearData[EulerTerm::V] = sqrt(V2);
  
  
  /// here I should pass an average of the vibrational temperatures
  library->frozenGammaAndSoundSpeed(Tdim, pdim, rhoDim, 
				    linearData[EulerTerm::GAMMA], 
				    linearData[EulerTerm::A],
				    CFNULL); // this sound speed is inconsistent
  
  
  linearData[EulerTerm::T] = T;
  linearData[EulerTerm::RHO] = rho;
  linearData[EulerTerm::H] = _avZ[TID]/(sqRho*refData[EulerTerm::H]);
  linearData[EulerTerm::E] = linearData[EulerTerm::H] - p/linearData[EulerTerm::RHO];
  
  const CFuint nbTv = _model->getNbScalarVars(1);
  const CFuint firstTv = _model->getFirstScalarVar(1);
  cf_assert(nbTv == 1);
  for (CFuint i = 0; i < nbTv; ++i) {
    linearData[firstTv+i] = _avZ[TID +1 +i]/(sqRho*refData[EulerTerm::H]);
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace NEQ
} // namespace Physics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
