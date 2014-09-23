#include "NEQ/NEQ.hh"
#include "Euler2DNEQLinearRhoivtCNEQ.hh"
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

Environment::ObjectProvider<Euler2DNEQLinearRhoivtCNEQ, JacobianLinearizer, NEQModule, 1>
euler2DNEQLinearRhoivtCNEQProvider("Euler2DNEQLinearRhoivtCNEQ");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQLinearRhoivtCNEQ::Euler2DNEQLinearRhoivtCNEQ(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->
	 getConvectiveTerm().d_castTo< MultiScalarTerm<EulerTerm> >()),
  _ye(),
  _dhe()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQLinearRhoivtCNEQ::~Euler2DNEQLinearRhoivtCNEQ()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQLinearRhoivtCNEQ::linearize(const vector<State*>& statesInCell)
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

  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  // Set the mixture density (sum of the partial densities)
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += _avZ[ie];
  }
  
  // Set the species
  const CFreal ovRho = 1./rho;
  _ye.resize(nbSpecies);
  // set the species mass fractions
  const CFuint firstSpecies = _model->getFirstScalarVar(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ye[ie] = _avZ[ie]*ovRho; 
    linearData[firstSpecies + ie] = _ye[ie];
  }
  
  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!! 
  library->setSpeciesFractions(_ye);
  
  const RealVector& refData =  _model->getReferencePhysicalData();
  CFreal rhoDim = rho*refData[EulerTerm::RHO];
  CFreal T = _avZ[nbSpecies + 2];
  CFreal Tdim = T*refData[EulerTerm::T];

  _dhe.resize(3);

  CFreal pdim = library->pressure(rhoDim, Tdim, CFNULL);
  //cout <<"LinearRhoivtCNEQ!" << endl;
  CFreal p = pdim/refData[EulerTerm::P];
  
  linearData[EulerTerm::P] = p;
  linearData[EulerTerm::VX] = _avZ[nbSpecies];
  linearData[EulerTerm::VY] = _avZ[nbSpecies+1];
  
  const CFreal V2 = linearData[EulerTerm::VX]*linearData[EulerTerm::VX] +
    linearData[EulerTerm::VY]*linearData[EulerTerm::VY];
  linearData[EulerTerm::V] = sqrt(V2);
  
  // unused //  const CFreal ovHref = 1./refData[EulerTerm::H];
  // vibrational temperatures
  //library->setDensityEnthalpyEnergy(Tdim, _tvDim, pdim,_dhe, true);
  library->setDensityEnthalpyEnergy(Tdim, pdim, _dhe);


  
  linearData[EulerTerm::T] = _avZ[nbSpecies+2];
  linearData[EulerTerm::RHO] = rho;

  library->frozenGammaAndSoundSpeed(Tdim, pdim, _dhe[0], linearData[EulerTerm::GAMMA], linearData[EulerTerm::A], CFNULL);//->??
//  library->gammaAndSoundSpeed(Tdim, pdim, _dhe[0],
//                      linearData[EulerTerm::GAMMA],
//                      linearData[EulerTerm::A]
//                      );


  linearData[EulerTerm::H] = _dhe[1]/refData[EulerTerm::H] + 0.5*V2;
  linearData[EulerTerm::E] = _dhe[2]/refData[EulerTerm::H] + 0.5*V2;
  
}

//////////////////////////////////////////////////////////////////////////////

} // namespace NEQ
} // namespace Physics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
