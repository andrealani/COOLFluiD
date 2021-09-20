#include "NEQ.hh"
#include "Euler3DNEQRhoivtToCons.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/EulerPhysicalModel.hh"
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

Environment::ObjectProvider<Euler3DNEQRhoivtToCons, VarSetTransformer, NEQModule, 1>
euler3DNEQRhoivtToConsProvider("Euler3DNEQRhoivtToCons");

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQRhoivtToCons::Euler3DNEQRhoivtToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerNEQTerm>()),
  _dhe(3),
  _ye()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQRhoivtToCons::~Euler3DNEQRhoivtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivtToCons::transform(const State& state, State& result)
{
  //cout<<"Euler3DNEQRhoivtToCons\n";
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const CFuint nbSpecies = _model->getNbScalarVars(0);
//   const EquationSubSysDescriptor& eqSS =
//     PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
//   const CFuint iEqSS = eqSS.getEqSS();
//   const CFuint nbEqSS = eqSS.getTotalNbEqSS();
  
  // partial densities
  _ye.resize(nbSpecies);
  CFreal rho = 0.0;
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = state[ie];
    rho += state[ie];
    _ye[ie] = state[ie];
  }
  _ye /= rho;
  
  const RealVector& refData =  _model->getReferencePhysicalData();
  CFreal rhoDim = rho*refData[EulerTerm::RHO];
  CFreal T = state[nbSpecies + 3]; // state[nbSpecies + 2]; // VS: this was 2 for 3D hence CNEQ3D didn't work, corrected to 3 for 3-D
  CFreal Tdim = T*refData[EulerTerm::T];
  CFreal* rhoi = &const_cast<State&>(state)[0];
  library->setState(rhoi, &Tdim);
  
  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!! 
  library->setSpeciesFractions(_ye);
  
  CFreal pdim = library->pressure(rhoDim, Tdim, CFNULL);
  const CFreal u = state[nbSpecies];
  const CFreal v = state[nbSpecies + 1];
  const CFreal w = state[nbSpecies + 2];
  const CFreal V2 = u*u + v*v + w*w;
  library->setDensityEnthalpyEnergy(Tdim, pdim,_dhe);
  
  result[nbSpecies] = rho*u;
  result[nbSpecies + 1] = rho*v;
  result[nbSpecies + 2] = rho*w;
  result[nbSpecies + 3] = rho*(_dhe[2]/refData[EulerTerm::H] + 0.5*V2);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQRhoivtToCons::transformFromRef(const RealVector& data, State& result) 
{
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFuint firstSpecies = _model->getFirstScalarVar(0);
  const CFreal rho = data[EulerTerm::RHO];
//   const EquationSubSysDescriptor& eqSS =
//     PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
//   const CFuint iEqSS = eqSS.getEqSS();
//   const CFuint nbEqSS = eqSS.getTotalNbEqSS();
  
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = rho*data[firstSpecies + ie];
  }
  
  result[nbSpecies] = rho*data[EulerTerm::VX];
  result[nbSpecies + 1] = rho*data[EulerTerm::VY];
  result[nbSpecies + 2] = rho*data[EulerTerm::VZ];
  result[nbSpecies + 3] = rho*data[EulerTerm::E];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
