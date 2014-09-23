#include "NEQ/NEQ.hh"
#include "Euler2DNEQLinearConsCNEQ.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"
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

Environment::ObjectProvider<Euler2DNEQLinearConsCNEQ, JacobianLinearizer, NEQModule, 1>
euler2DNEQLinearConsCNEQProvider("Euler2DNEQLinearConsCNEQ");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQLinearConsCNEQ::Euler2DNEQLinearConsCNEQ(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->
	 getConvectiveTerm().d_castTo< MultiScalarTerm<EulerTerm> >()),
  _ye(),
  _dhe(),  
  _mmasses()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQLinearConsCNEQ::~Euler2DNEQLinearConsCNEQ()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQLinearConsCNEQ::linearize(const vector<State*>& statesInCell)
{
  RealVector& linearData = _model->getPhysicalData();

  const CFuint nbStates = statesInCell.size();
  
  // linearizer does this (numerics dependent)
  //Getting Conservatives variable Average
  _sumZ = 0.;
  for (CFuint iEq = 0; iEq < PhysicalModelStack::getActive()->getNbEq(); ++iEq) {
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      _sumZ[iEq] += (*statesInCell[iState])[iEq];
    }
  }
  _avZ = _sumZ/static_cast<CFreal>(nbStates);

  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  SafePtr<PhysicalChemicalLibrary::ExtraData> eData = library->getExtraData();
   
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  _mmasses.resize(nbSpecies);
  library->getMolarMasses(_mmasses);
  SafePtr<RealVector> fcoeff = library->getAtomicityCoeff();

  // Set the mixture density (sum of the partial densities)
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += _avZ[ie];
  }
  
  linearData[EulerTerm::RHO] = rho;
  
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
  
  // U  and V Velocity Average
  linearData[EulerTerm::VX] = _avZ[nbSpecies]*ovRho;
  linearData[EulerTerm::VY] = _avZ[nbSpecies+1]*ovRho;
  
  const CFreal V2 = linearData[EulerTerm::VX]*linearData[EulerTerm::VX] +
    linearData[EulerTerm::VY]*linearData[EulerTerm::VY];
  linearData[EulerTerm::V] = sqrt(V2);
      
  //Temperature Average
  linearData[EulerTerm::E]  = _avZ[nbSpecies+2]*ovRho;
  const CFreal Rgas = library->getRgas();
  CFreal denom = 0.;
  CFreal form  = 0.;
  CFreal riovermi = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal yOvM = _ye[i]/_mmasses[i];
    riovermi += _avZ[i]/_mmasses[i];
    denom += yOvM*(Rgas*(*fcoeff)[i]);
    form += _ye[i]*eData->enthalpyForm[i];
  }
  
  // here could be changed for ionized case ....
  linearData[EulerTerm::T] = (linearData[EulerTerm::E] - form - 0.5*V2)/denom;
  
  //Enthalpy Average
  const CFreal P = linearData[EulerTerm::T]*Rgas*riovermi;
  linearData[EulerTerm::P]= P; 
  linearData[EulerTerm::H]= linearData[EulerTerm::E] + P*ovRho;
  
  //Speed of Sound  Average
  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  
  const CFuint start = (library->presenceElectron()) ? 1 : 0;
  for (CFuint i = start; i < nbSpecies; ++i) {
    const CFreal sigmai = linearData[firstSpecies + i]/_mmasses[i];
    numBeta += sigmai;
    denBeta += sigmai*(*fcoeff)[i];
  }
  
  const CFreal beta = numBeta/denBeta;
  const CFreal RT = Rgas*linearData[EulerTerm::T];
  eData->dEdT = 0.0;

  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal yOvM = _ye[i]/_mmasses[i]; 
    eData->energyTr[i] = (*fcoeff)[i]*RT/_mmasses[i] + eData->enthalpyForm[i];
    eData->dEdT += yOvM*(*fcoeff)[i]*Rgas;
  }
  
 linearData[EulerTerm::A] = std::sqrt((1. + beta)*linearData[EulerTerm::P]/linearData[EulerTerm::RHO]);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace NEQ
} // namespace Physics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
